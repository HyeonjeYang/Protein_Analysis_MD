"""Config normalization, preset resolution, and lock-file compilation."""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import yaml

from idrptm.presets import (
    analysis_preset,
    merge_overrides,
    report_preset,
    simulation_preset,
)
from idrptm.provenance import build_trajectory_folder_name
from idrptm.schema import WorkflowConfig
from idrptm.storage import estimate_project_storage
from idrptm.units import CANONICAL_UNITS

LOCK_FILE_NAME = "project.lock.yaml"
RESOLVED_FILE_NAME = "config_resolved.json"


@dataclass(frozen=True)
class NormalizedConfig:
    """A normalized user-facing config before preset expansion."""

    source_path: Path | None
    project_name: str
    project_dir: Path
    input: dict[str, Any]
    protocol: dict[str, Any]
    analysis: dict[str, Any]
    report: dict[str, Any]
    execution: dict[str, Any]
    trajectory: dict[str, Any]
    visualization: bool
    sweep: dict[str, Any]
    legacy_workflow: dict[str, Any] | None = None


@dataclass(frozen=True)
class ResolvedConfig:
    """A preset-resolved config ready to be converted to WorkflowConfig YAML."""

    normalized: NormalizedConfig
    simulation: dict[str, Any]
    environment: dict[str, Any]
    analysis: dict[str, Any]
    report: dict[str, Any]
    execution: dict[str, Any]
    workflow: dict[str, Any]
    warnings: tuple[str, ...] = ()


@dataclass(frozen=True)
class LockedConfig:
    """Compiled project lock paths and metadata."""

    project_dir: Path
    lock_yaml: Path
    resolved_json: Path
    workflow: WorkflowConfig
    resolved: dict[str, Any]
    warnings: tuple[str, ...] = ()


def resolve_config_target(target: str | Path) -> Path:
    """Resolve a config path or compiled project directory to a workflow YAML file."""

    path = Path(target)
    if path.is_dir():
        lock_path = path / LOCK_FILE_NAME
        if not lock_path.exists():
            raise FileNotFoundError(
                f"{path} is a directory but does not contain {LOCK_FILE_NAME}. "
                "Run 'pamd compile CONFIG.yaml' first."
            )
        return lock_path
    return path


def load_user_config(path: str | Path) -> dict[str, Any]:
    """Load a user YAML config."""

    config_path = Path(path)
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    if not isinstance(data, dict):
        raise ValueError(f"Config {config_path} must contain a YAML mapping.")
    return data


def normalize_config(
    user_config: dict[str, Any],
    *,
    source_path: str | Path | None = None,
) -> NormalizedConfig:
    """Normalize old explicit configs and new short configs into one shape."""

    path = Path(source_path) if source_path is not None else None
    if _is_legacy_workflow(user_config):
        project_name = _legacy_project_name(user_config)
        project_dir = Path(user_config.get("runner", {}).get("work_dir", f"runs/{project_name}"))
        trajectory = _legacy_trajectory_options(user_config)
        return NormalizedConfig(
            source_path=path,
            project_name=project_name,
            project_dir=project_dir,
            input=_legacy_input(user_config),
            protocol=_legacy_protocol(user_config),
            analysis={"preset": "custom", "overrides": user_config.get("analysis", {})},
            report={"preset": "standard"},
            execution=dict(user_config.get("execution") or {}),
            trajectory=trajectory,
            visualization=bool(user_config.get("visualization", True)),
            sweep=user_config.get("sweep", {}),
            legacy_workflow=dict(user_config),
        )

    project = user_config.get("project", {})
    if isinstance(project, str):
        project_name = project
        project_payload: dict[str, Any] = {}
    else:
        project_name = str(project.get("name") or "protein_analysis_md_project")
        project_payload = dict(project)
    input_block = dict(user_config.get("input") or {})
    if "protein" not in input_block and "protein" in user_config:
        input_block["protein"] = user_config["protein"]
    if path is not None:
        input_block = _resolve_input_paths(input_block, path.parent)
    protocol = dict(user_config.get("protocol") or {})
    protocol.setdefault("preset", "smoke_single_chain")
    if path is not None:
        protocol = _resolve_protocol_paths(protocol, path.parent)
    trajectory = _trajectory_options(
        project_name=project_name,
        project_payload=project_payload,
        user_config=user_config,
        input_block=input_block,
        protocol=protocol,
    )
    explicit_outdir = (
        user_config.get("outdir")
        if isinstance(project, str)
        else project_payload.get("outdir") or user_config.get("outdir")
    )
    project_dir = Path(explicit_outdir) if explicit_outdir else Path(trajectory["work_dir"])
    return NormalizedConfig(
        source_path=path,
        project_name=project_name,
        project_dir=project_dir,
        input=input_block,
        protocol=protocol,
        analysis=dict(user_config.get("analysis") or {"preset": "standard_idr"}),
        report=dict(user_config.get("report") or {"preset": "standard"}),
        execution=dict(user_config.get("execution") or {}),
        trajectory=trajectory,
        visualization=_visualization_enabled(user_config),
        sweep=dict(user_config.get("sweep") or {}),
    )


def _resolve_input_paths(input_block: dict[str, Any], base_dir: Path) -> dict[str, Any]:
    resolved = dict(input_block)
    protein = dict(resolved.get("protein") or {})
    if protein.get("fasta"):
        fasta = Path(str(protein["fasta"]))
        if not fasta.is_absolute():
            protein["fasta"] = str((base_dir / fasta).resolve())
    resolved["protein"] = protein
    ptm = dict(resolved.get("ptm") or {})
    if ptm.get("file"):
        ptm_file = Path(str(ptm["file"]))
        if not ptm_file.is_absolute():
            ptm["file"] = str((base_dir / ptm_file).resolve())
    if ptm:
        resolved["ptm"] = ptm
    return resolved


def _resolve_protocol_paths(protocol: dict[str, Any], base_dir: Path) -> dict[str, Any]:
    resolved = dict(protocol)
    simulation = dict(resolved.get("simulation") or {})
    if simulation.get("residue_parameters"):
        residue_parameters = Path(str(simulation["residue_parameters"]))
        if not residue_parameters.is_absolute():
            simulation["residue_parameters"] = str((base_dir / residue_parameters).resolve())
    if simulation:
        resolved["simulation"] = simulation
    return resolved


def resolve_presets(normalized_config: NormalizedConfig) -> ResolvedConfig:
    """Resolve preset names and user overrides into complete settings."""

    warnings: list[str] = []
    if normalized_config.legacy_workflow is not None:
        workflow = _normalize_legacy_workflow(normalized_config.legacy_workflow)
        protocol = normalized_config.protocol
        return ResolvedConfig(
            normalized=normalized_config,
            simulation=protocol.get("simulation", {}),
            environment=protocol.get("environment", {}),
            analysis=workflow.get("analysis", {}),
            report=normalized_config.report,
            execution=workflow.get("execution", {}),
            workflow=workflow,
            warnings=tuple(warnings),
        )

    protocol = normalized_config.protocol
    simulation_name = str(protocol.get("preset") or "smoke_single_chain")
    simulation = simulation_preset(simulation_name)
    simulation = merge_overrides(simulation, protocol.get("simulation"))
    environment = merge_overrides(
        {
            "pH": 7.4,
            "ionic_M": 0.15,
            "temperature_K": 298.0,
        },
        protocol.get("environment"),
    )
    analysis_name = str(normalized_config.analysis.get("preset") or "standard_idr")
    analysis = analysis_preset(analysis_name)
    analysis = merge_overrides(analysis, normalized_config.analysis.get("overrides"))
    report_name = str(normalized_config.report.get("preset") or "standard")
    report = report_preset(report_name)
    report = merge_overrides(
        report,
        {key: value for key, value in normalized_config.report.items() if key != "preset"},
    )
    workflow = _workflow_from_normalized(
        normalized=normalized_config,
        simulation=simulation,
        environment=environment,
        analysis=analysis,
        report=report,
    )
    return ResolvedConfig(
        normalized=normalized_config,
        simulation=simulation,
        environment=environment,
        analysis=analysis,
        report=report,
        execution=normalized_config.execution,
        workflow=workflow,
        warnings=tuple(warnings),
    )


def compile_config(resolved_config: ResolvedConfig) -> LockedConfig:
    """Compile a resolved config into a validated locked WorkflowConfig object."""

    workflow = WorkflowConfig.model_validate(resolved_config.workflow)
    storage_payload: dict[str, Any]
    warnings = list(resolved_config.warnings)
    try:
        if not _storage_can_resolve_without_prompt(workflow):
            raise ValueError(
                "sequence source requires ambiguous UniProt selection; storage estimate "
                "will be available after using an accession or cached FASTA"
            )
        storage = estimate_project_storage(workflow)
        storage_payload = storage.to_dict()
    except Exception as exc:
        storage_payload = {"warning": f"Storage estimate skipped: {exc}"}
        warnings.append(storage_payload["warning"])

    workflow_payload = workflow.model_dump(mode="json")
    workflow_payload["sequence"] = None
    workflow_payload["protein"] = None
    workflow_payload["compiled"] = {
        "source_config": str(resolved_config.normalized.source_path)
        if resolved_config.normalized.source_path
        else None,
        "simulation_preset": resolved_config.normalized.protocol.get("preset", "custom"),
        "analysis_preset": resolved_config.normalized.analysis.get("preset", "custom"),
        "report_preset": resolved_config.normalized.report.get("preset", "standard"),
        "simulation": resolved_config.simulation,
        "environment": resolved_config.environment,
        "report": resolved_config.report,
        "execution": resolved_config.execution,
        "trajectory": resolved_config.normalized.trajectory,
        "visualization": resolved_config.normalized.visualization,
        "sweep": resolved_config.normalized.sweep,
        "warnings": warnings,
    }
    workflow_payload["units"] = CANONICAL_UNITS
    workflow_payload["storage_estimate"] = storage_payload
    locked_workflow = WorkflowConfig.model_validate(workflow_payload)

    resolved = {
        "project": {
            "name": resolved_config.normalized.project_name,
            "outdir": str(resolved_config.normalized.project_dir),
        },
        "workflow": workflow_payload,
        "simulation": resolved_config.simulation,
        "environment": resolved_config.environment,
        "analysis": resolved_config.analysis,
        "report": resolved_config.report,
        "execution": resolved_config.execution,
        "trajectory": resolved_config.normalized.trajectory,
        "visualization": resolved_config.normalized.visualization,
        "sweep": resolved_config.normalized.sweep,
        "units": CANONICAL_UNITS,
        "storage_estimate": storage_payload,
        "warnings": warnings,
    }
    return LockedConfig(
        project_dir=resolved_config.normalized.project_dir,
        lock_yaml=resolved_config.normalized.project_dir / LOCK_FILE_NAME,
        resolved_json=resolved_config.normalized.project_dir / RESOLVED_FILE_NAME,
        workflow=locked_workflow,
        resolved=resolved,
        warnings=tuple(warnings),
    )


def compile_config_file(
    config_path: str | Path,
    *,
    force: bool = False,
) -> LockedConfig:
    """Compile a user config file and write project lock artifacts."""

    source = Path(config_path)
    normalized = normalize_config(load_user_config(source), source_path=source)
    resolved = resolve_presets(normalized)
    locked = compile_config(resolved)
    write_locked_config(locked, force=force)
    return locked


def write_locked_config(locked: LockedConfig, *, force: bool = False) -> None:
    """Write lock YAML, resolved JSON, storage estimate, and manifest preview."""

    locked.project_dir.mkdir(parents=True, exist_ok=True)
    if locked.lock_yaml.exists() and not force:
        raise FileExistsError(f"{locked.lock_yaml} already exists; pass --force to overwrite.")
    workflow_payload = locked.workflow.model_dump(mode="json")
    locked.lock_yaml.write_text(
        yaml.safe_dump(workflow_payload, sort_keys=False),
        encoding="utf-8",
    )
    locked.resolved_json.write_text(
        json.dumps(locked.resolved, indent=2, default=str) + "\n",
        encoding="utf-8",
    )
    storage = locked.resolved.get("storage_estimate", {})
    if isinstance(storage, dict) and "runs" in storage:
        write_storage_estimate_from_dict(storage, locked.project_dir / "storage_estimate.json")
        _write_manifest_preview(storage, locked.project_dir / "manifest_preview.csv")


def write_storage_estimate_from_dict(payload: dict[str, Any], path: Path) -> None:
    """Write a storage estimate payload without rebuilding design variants."""

    path.write_text(json.dumps(payload, indent=2, default=str) + "\n", encoding="utf-8")


def _write_manifest_preview(storage_payload: dict[str, Any], path: Path) -> None:
    rows = storage_payload.get("runs", [])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["run_id", "n_beads", "n_frames", "dcd_mb"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "run_id": row.get("run_id"),
                    "n_beads": row.get("n_beads"),
                    "n_frames": row.get("n_frames"),
                    "dcd_mb": row.get("dcd_mb"),
                }
            )


def _is_legacy_workflow(data: dict[str, Any]) -> bool:
    return any(key in data for key in ("sequence", "protein", "proteins", "calvados", "runner"))


def _legacy_project_name(data: dict[str, Any]) -> str:
    project = data.get("project", "protein_analysis_md_project")
    if isinstance(project, dict):
        return str(project.get("name") or "protein_analysis_md_project")
    return str(project)


def _legacy_input(data: dict[str, Any]) -> dict[str, Any]:
    if "protein" in data:
        return {"protein": data["protein"], "ptm": data["protein"].get("ptm", data.get("ptm", {}))}
    if "sequence" in data:
        return {"protein": data["sequence"], "ptm": data.get("ptm", {})}
    if "proteins" in data:
        return {"proteins": data["proteins"]}
    return {}


def _legacy_trajectory_options(data: dict[str, Any]) -> dict[str, Any]:
    runner = dict(data.get("runner") or {})
    return {
        "traj_name": runner.get("traj_name"),
        "traj_flag": runner.get("traj_flag"),
        "trajectory_root": str(runner.get("trajectory_root", "runs")),
        "include_timestamp": bool(runner.get("include_timestamp", False)),
        "timestamp_format": runner.get("timestamp_format", "%Y%m%d_%H%M%S"),
        "work_dir": str(runner.get("work_dir", f"runs/{_legacy_project_name(data)}")),
    }


def _trajectory_options(
    *,
    project_name: str,
    project_payload: dict[str, Any],
    user_config: dict[str, Any],
    input_block: dict[str, Any],
    protocol: dict[str, Any],
) -> dict[str, Any]:
    trajectory = dict(user_config.get("trajectory") or {})
    traj_name = (
        project_payload.get("traj_name")
        or project_payload.get("trajectory_name")
        or trajectory.get("traj_name")
        or trajectory.get("trajectory_name")
    )
    traj_flag = (
        project_payload.get("traj_flag")
        or project_payload.get("flag")
        or trajectory.get("traj_flag")
        or trajectory.get("flag")
    )
    trajectory_root = Path(
        project_payload.get("trajectory_root")
        or trajectory.get("trajectory_root")
        or user_config.get("trajectory_root")
        or "runs"
    )
    include_timestamp = bool(
        project_payload.get(
            "include_timestamp",
            trajectory.get("include_timestamp", True),
        )
    )
    timestamp_format = str(
        project_payload.get(
            "timestamp_format",
            trajectory.get("timestamp_format", "%Y%m%d_%H%M%S"),
        )
    )
    simulation = dict(protocol.get("simulation") or {})
    production = dict(simulation.get("production") or {})
    ptm = dict(input_block.get("ptm") or {})
    cleavage = dict(input_block.get("cleavage") or {})
    folder = build_trajectory_folder_name(
        project_name=project_name,
        protein_hint=_protein_hint(input_block),
        preset=str(protocol.get("preset") or "custom"),
        total_time_ns=production.get("total_time_ns"),
        replicates=simulation.get("replicates"),
        ptm_mode=ptm.get("mode"),
        cleavage_mode=cleavage.get("mode"),
        traj_name=str(traj_name) if traj_name else None,
        traj_flag=str(traj_flag) if traj_flag else None,
        include_timestamp=include_timestamp,
        timestamp_format=timestamp_format,
    )
    return {
        "traj_name": traj_name,
        "traj_flag": traj_flag,
        "trajectory_root": str(trajectory_root),
        "include_timestamp": include_timestamp,
        "timestamp_format": timestamp_format,
        "work_dir": str(trajectory_root / folder),
    }


def _visualization_enabled(user_config: dict[str, Any]) -> bool:
    if "visualization" in user_config:
        return bool(user_config["visualization"])
    report = user_config.get("report")
    if isinstance(report, dict) and "visualization" in report:
        return bool(report["visualization"])
    return True


def _protein_hint(input_block: dict[str, Any]) -> str | None:
    protein = dict(input_block.get("protein") or {})
    if protein:
        return str(
            protein.get("accession")
            or protein.get("name")
            or protein.get("query")
            or ""
        )
    proteins = input_block.get("proteins")
    if isinstance(proteins, list) and proteins:
        names = [
            str(item.get("accession") or item.get("name") or item.get("query") or "")
            for item in proteins
            if isinstance(item, dict)
        ]
        return "_".join(name for name in names if name) or None
    return None


def _legacy_protocol(data: dict[str, Any]) -> dict[str, Any]:
    calvados = data.get("calvados", {})
    return {
        "preset": "custom",
        "environment": {
            "pH": calvados.get("ph"),
            "ionic_M": calvados.get("ionic_m"),
            "temperature_K": calvados.get("temperature_k"),
        },
        "simulation": calvados.get("simulation", {}),
    }


def _normalize_legacy_workflow(data: dict[str, Any]) -> dict[str, Any]:
    workflow = dict(data)
    project = workflow.get("project")
    if isinstance(project, dict):
        workflow["project"] = project.get("name") or "protein_analysis_md_project"
        workflow.setdefault("runner", {})
        if project.get("outdir") and "work_dir" not in workflow["runner"]:
            workflow["runner"]["work_dir"] = project["outdir"]
    calvados = dict(workflow.get("calvados") or {})
    if "platform" in calvados:
        simulation = dict(calvados.get("simulation") or {})
        simulation.setdefault("platform", calvados.pop("platform"))
        calvados["simulation"] = simulation
    workflow["calvados"] = calvados
    workflow.setdefault("execution", {})
    workflow.setdefault("compiled", {})
    workflow.setdefault("units", CANONICAL_UNITS)
    workflow.setdefault("storage_estimate", {})
    return workflow


def _workflow_from_normalized(
    *,
    normalized: NormalizedConfig,
    simulation: dict[str, Any],
    environment: dict[str, Any],
    analysis: dict[str, Any],
    report: dict[str, Any],
) -> dict[str, Any]:
    protein = _protein_from_input(normalized.input)
    production = dict(simulation.get("production") or {})
    workflow = {
        "project": normalized.project_name,
        "replicates": int(simulation.get("replicates", 1)),
        "visualization": normalized.visualization,
        "protein": protein,
        "calvados": {
            "model": simulation.get("model", "CALVADOS2"),
            "residue_parameters": simulation.get("residue_parameters"),
            "box_nm": simulation.get("box_nm", [34.4, 34.4, 34.4]),
            "temperature_k": environment.get("temperature_K", 298.0),
            "ph": environment.get("pH", 7.4),
            "ionic_m": environment.get("ionic_M", 0.15),
            "topol": simulation.get("topol", "center"),
            "simulation": {
                "integrator": "calvados_default",
                "dt_ps": simulation.get("dt_ps", 0.01),
                "total_time_ns": production.get("total_time_ns", 0.5),
                "frame_interval_ns": production.get("frame_interval_ns", 0.1),
                "runtime_hours": simulation.get("runtime_hours", 0),
                "platform": simulation.get("platform", "CPU"),
                "restart": simulation.get("restart", "checkpoint"),
                "checkpoint_file": simulation.get("checkpoint_file", "restart.chk"),
                "random_seed": simulation.get("random_seed"),
            },
        },
        "runner": {
            "backend": simulation.get("backend", "local"),
            "work_dir": str(normalized.project_dir),
            "dry_run": True,
            "traj_name": normalized.trajectory.get("traj_name"),
            "traj_flag": normalized.trajectory.get("traj_flag"),
            "trajectory_root": normalized.trajectory.get("trajectory_root", "runs"),
            "include_timestamp": normalized.trajectory.get("include_timestamp", True),
            "timestamp_format": normalized.trajectory.get("timestamp_format", "%Y%m%d_%H%M%S"),
            "progress": bool(simulation.get("progress", True)),
            "progress_interval_s": float(simulation.get("progress_interval_s", 5.0)),
        },
        "execution": normalized.execution,
        "analysis": _analysis_to_workflow(analysis),
        "compiled": {
            "report": report,
        },
        "units": CANONICAL_UNITS,
        "storage_estimate": {},
    }
    return workflow


def _protein_from_input(input_block: dict[str, Any]) -> dict[str, Any]:
    protein = dict(input_block.get("protein") or {})
    if not protein:
        raise ValueError("New config format requires input.protein.")
    protein["ptm"] = _ptm_from_input(input_block.get("ptm") or {})
    cleavage = input_block.get("cleavage") or {"mode": "none"}
    protein["cleavage_sets"] = [_cleavage_from_input(cleavage)]
    return protein


def _ptm_from_input(ptm: dict[str, Any]) -> dict[str, Any]:
    mode = ptm.get("mode", "none")
    if mode in {"none", "uniprot_suggestions"}:
        return {"mode": "wt", "include_wt": True, "sites": []}
    if mode == "combinatorial":
        return {
            "mode": "all_sites" if ptm.get("include_all_sites", True) else "single_site_scan",
            "include_wt": True,
            "sites": [_ptm_site(site) for site in ptm.get("sites", [])],
        }
    if mode == "file":
        from idrptm.ptm import parse_ptm_file

        path = Path(ptm["file"])
        sites = [
            {
                "position": request.biological_position,
                "residue": request.expected_residue,
                "ptm": request.ptm.name,
            }
            for request in parse_ptm_file(path)
        ]
        return {"mode": "explicit", "include_wt": True, "sites": sites}
    states = ptm.get("states", [])
    for state in states:
        modifications = state.get("modifications", [])
        if modifications:
            return {
                "mode": "explicit",
                "include_wt": any(not item.get("modifications") for item in states),
                "sites": [_ptm_site(site) for site in modifications],
            }
    return {"mode": "wt", "include_wt": True, "sites": []}


def _ptm_site(site: dict[str, Any]) -> dict[str, Any]:
    return {
        "position": site.get("position", site.get("site")),
        "residue": site["residue"],
        "ptm": site["ptm"],
    }


def _cleavage_from_input(cleavage: dict[str, Any]) -> dict[str, Any]:
    mode = cleavage.get("mode", "none")
    payload = {
        "name": cleavage.get("name", "intact" if mode == "none" else mode),
        "mode": mode,
        "include_intact": mode == "none",
        "individual_fragments": False if mode == "none" else True,
        "fragment_mixture": False if mode == "none" else True,
    }
    for key in (
        "enzyme",
        "protease",
        "cut_after_sites",
        "manual_cuts",
        "n_cuts",
        "seed",
        "min_fragment_length",
        "output",
        "candidate_sites",
        "rate_model",
        "global_rate_per_ns",
        "site_rate_per_ns",
        "max_time_ns",
        "max_cuts",
        "terminus",
        "step_size",
        "max_removed",
    ):
        if key in cleavage:
            payload[key] = cleavage[key]
    return payload


def _analysis_to_workflow(analysis: dict[str, Any]) -> dict[str, Any]:
    supported = {
        "rg",
        "ree",
        "contacts",
        "ps",
        "scaling",
        "lifetime",
        "contact_lifetime",
        "msd",
    }
    observables = [item for item in analysis.get("observables", []) if item in supported]
    if not observables:
        observables = ["rg", "ree", "contacts", "ps", "scaling"]
    payload: dict[str, Any] = {"observables": observables}
    for key in (
        "contact_cutoff_nm",
        "min_sequence_separation",
        "max_lag",
        "fit_min_s",
        "fit_max_s",
        "fit_to",
        "smoothing",
        "decomposition",
        "free_energy",
    ):
        if key in analysis:
            payload[key] = analysis[key]
    return payload


def _storage_can_resolve_without_prompt(workflow: WorkflowConfig) -> bool:
    for protein in workflow.proteins:
        if protein.source == "uniprot" and not protein.accession:
            return False
    return True
