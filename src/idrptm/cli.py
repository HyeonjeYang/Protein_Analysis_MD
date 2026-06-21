"""Command-line interface for protein_analysis_md."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Annotated, Literal, cast

import typer
import yaml

from idrptm import __version__
from idrptm.runner import RunPhase

app = typer.Typer(
    name="pamd",
    help="Prepare, run, and analyze CALVADOS-backed protein/IDR MD workflows.",
    no_args_is_help=True,
)


ConfigOption = Annotated[
    Path | None,
    typer.Option("--config", "-c", help="Workflow YAML configuration file."),
]
TrajectoryReader = Literal["mdtraj", "mdanalysis"]


@app.callback()
def callback(
    version: Annotated[
        bool,
        typer.Option("--version", help="Show the protein_analysis_md version and exit."),
    ] = False,
) -> None:
    """Top-level CLI callback."""

    if version:
        typer.echo(f"protein_analysis_md {__version__}")
        raise typer.Exit()


@app.command("init")
def init_command(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Directory where a starter project would be created."),
    ] = Path("."),
    force: Annotated[
        bool,
        typer.Option("--force", help="Allow overwriting existing starter files."),
    ] = False,
) -> None:
    """Placeholder for creating a starter workflow directory."""

    typer.echo(
        f"Would initialize a protein_analysis_md project at {output} "
        f"(force={force})."
    )


@app.command("search-uniprot")
def search_uniprot_command(
    query: Annotated[str, typer.Argument(help="UniProt search query.")],
    reviewed: Annotated[
        bool,
        typer.Option("--reviewed/--all", help="Restrict to Swiss-Prot reviewed entries."),
    ] = True,
    organism: Annotated[
        str | None,
        typer.Option("--organism", help="Organism scientific name."),
    ] = None,
    size: Annotated[int, typer.Option("--size", help="Number of candidates to show.")] = 10,
) -> None:
    """Search UniProt/Swiss-Prot and print ranked candidates."""

    from idrptm.uniprot import format_candidates, search_uniprot

    try:
        candidates = search_uniprot(query, reviewed=reviewed, organism=organism, size=size)
    except Exception as exc:
        typer.echo(f"UniProt search failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(format_candidates(candidates))


@app.command("fetch-sequence")
def fetch_sequence_command(
    query: Annotated[str, typer.Argument(help="UniProt query or accession.")],
    reviewed: Annotated[
        bool,
        typer.Option("--reviewed/--all", help="Restrict to Swiss-Prot reviewed entries."),
    ] = True,
    organism: Annotated[
        str | None,
        typer.Option("--organism", help="Organism scientific name."),
    ] = None,
    accession: Annotated[
        str | None,
        typer.Option("--accession", help="Explicit UniProt accession."),
    ] = None,
    interactive: Annotated[
        bool,
        typer.Option("--interactive", help="Prompt to select among ambiguous candidates."),
    ] = False,
    yes: Annotated[
        bool,
        typer.Option("--yes", "-y", help="Accept a strong unique candidate non-interactively."),
    ] = False,
    refresh: Annotated[
        bool,
        typer.Option("--refresh", help="Refresh cached sequence metadata."),
    ] = False,
    cache_dir: Annotated[
        Path,
        typer.Option("--cache-dir", help="Sequence cache directory."),
    ] = Path("data/sequences"),
) -> None:
    """Fetch one UniProt/Swiss-Prot sequence into the local cache."""

    from idrptm.uniprot import fetch_sequence

    try:
        result = fetch_sequence(
            query,
            accession=accession,
            reviewed=reviewed,
            organism=organism,
            interactive=interactive,
            yes=yes,
            refresh=refresh,
            cache_dir=cache_dir,
        )
    except Exception as exc:
        typer.echo(f"UniProt fetch failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Fetched {result.candidate.accession} {result.candidate.entry_name}")
    typer.echo(f"FASTA: {result.fasta_path}")
    typer.echo(f"Metadata: {result.metadata_path}")


@app.command("estimate-size")
def estimate_size_command(
    config: Annotated[Path, typer.Argument(help="Workflow YAML configuration file.")],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for storage_estimate.json."),
    ] = None,
) -> None:
    """Estimate trajectory storage for a workflow before running simulations."""

    from idrptm.configuration import resolve_config_target
    from idrptm.storage import estimate_from_config_file, format_storage_table

    try:
        estimate = estimate_from_config_file(resolve_config_target(config), output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Storage estimate failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(format_storage_table(estimate))


@app.command("env-check")
def env_check_command(
    json_output: Annotated[
        bool,
        typer.Option("--json", help="Also write machine-readable environment_check.json."),
    ] = False,
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="JSON output path when --json is used."),
    ] = Path("environment_check.json"),
) -> None:
    """Show whether commands appear to run locally or on a Remote-SSH host."""

    from protein_analysis_md.environment_check import (
        collect_environment_info,
        format_environment_report,
        write_environment_json,
    )

    info = collect_environment_info()
    typer.echo(format_environment_report(info))
    if json_output:
        path = write_environment_json(info, output)
        typer.echo(f"Wrote JSON: {path}")


@app.command("compile")
def compile_command(
    config: Annotated[Path, typer.Argument(help="Short or explicit workflow YAML.")],
    force: Annotated[
        bool,
        typer.Option("--force", help="Overwrite existing project.lock.yaml."),
    ] = False,
) -> None:
    """Compile a concise config into locked project configuration files."""

    from idrptm.configuration import compile_config_file

    try:
        locked = compile_config_file(config, force=force)
    except Exception as exc:
        typer.echo(f"Compile failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    storage = locked.resolved.get("storage_estimate", {})
    total_gb = storage.get("total_dcd_gb") if isinstance(storage, dict) else None
    typer.echo(f"Project: {locked.workflow.project}")
    typer.echo(f"Project dir: {locked.project_dir}")
    typer.echo(f"Lock: {locked.lock_yaml}")
    typer.echo(f"Resolved JSON: {locked.resolved_json}")
    typer.echo(f"Simulation preset: {locked.workflow.compiled.get('simulation_preset', 'custom')}")
    typer.echo(f"Analysis preset: {locked.workflow.compiled.get('analysis_preset', 'custom')}")
    if total_gb is not None:
        typer.echo(f"Estimated DCD storage: {float(total_gb):.4f} GB")
    for warning in locked.warnings:
        typer.echo(f"Warning: {warning}")


@app.command("design")
def design_command(
    config: Annotated[
        Path,
        typer.Argument(help="Workflow YAML configuration file."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for design outputs."),
    ] = None,
) -> None:
    """Generate WT/PTM design manifest, FASTA files, and metadata stubs."""

    from idrptm.configuration import resolve_config_target
    from idrptm.design import design_from_config_file

    try:
        result = design_from_config_file(resolve_config_target(config), output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Design failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote manifest: {result.manifest_path}")
    typer.echo(f"Wrote {len(result.fasta_paths)} variant FASTA file(s).")
    typer.echo(f"Wrote {len(result.metadata_paths)} per-run metadata stub(s).")


@app.command("prepare")
def prepare_command(
    config: Annotated[
        Path,
        typer.Argument(help="Workflow YAML configuration file."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="CALVADOS run-directory root."),
    ] = None,
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run", help="Print planned run directories without writing files."),
    ] = False,
) -> None:
    """Prepare CALVADOS run directories from a workflow config."""

    from idrptm.calvados_adapter import prepare_from_config_file
    from idrptm.configuration import resolve_config_target

    try:
        result = prepare_from_config_file(
            resolve_config_target(config),
            output_dir=output_dir,
            dry_run=dry_run,
        )
    except Exception as exc:
        typer.echo(f"Prepare failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    if result.dry_run:
        typer.echo(f"Dry run: would prepare {len(result.run_directories)} run directories.")
        for run_dir in result.run_directories:
            typer.echo(f"Would prepare: {run_dir.path}")
        return

    typer.echo(f"Prepared {len(result.run_directories)} CALVADOS run directories.")
    typer.echo(f"Manifest: {result.manifest_path}")


@app.command("run")
def run_command(
    target: Annotated[
        Path,
        typer.Argument(help="Prepared run directory, or project directory with --all-runs."),
    ] = Path("."),
    phase: Annotated[
        str,
        typer.Option("--phase", help="Run phase: equilibration, production, or all."),
    ] = "all",
    all_runs: Annotated[
        bool,
        typer.Option("--all-runs", help="Run every prepared run directory in a project."),
    ] = False,
    python_executable: Annotated[
        str,
        typer.Option("--python", help="Python executable used to invoke generated run scripts."),
    ] = "python",
    dry_run: Annotated[
        bool,
        typer.Option(
            "--dry-run/--execute",
            help="Print planned commands instead of executing them.",
        ),
    ] = True,
) -> None:
    """Run prepared CALVADOS directories locally through generated run scripts."""

    from idrptm.runner import execute_local_runs, plan_local_runs, write_planned_status
    from protein_analysis_md.environment_check import collect_environment_info

    try:
        if phase not in {"equilibration", "production", "all"}:
            raise ValueError("Phase must be 'equilibration', 'production', or 'all'.")
        plans = plan_local_runs(
            target,
            phase=cast(RunPhase, phase),
            all_runs=all_runs,
            python_executable=python_executable,
        )
        env_info = collect_environment_info(cwd=target)
        _print_run_environment_summary(target, plans, env_info)
        _enforce_execution_policy(target, env_info)
        if dry_run:
            typer.echo(f"Dry run: planned {len(plans)} local run(s).")
            for plan in plans:
                write_planned_status(plan)
                typer.echo(f"{plan.run_dir}: {' '.join(plan.command)}")
            return
        _confirm_local_nontrivial_run(target, env_info)
        results = execute_local_runs(plans)
    except Exception as exc:
        typer.echo(f"Run failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    failed = [result for result in results if result.status == "failed"]
    for result in results:
        typer.echo(f"{result.run_dir}: {result.status} ({result.status_json})")
    if failed:
        raise typer.Exit(1)


def _print_run_environment_summary(
    target: Path,
    plans: tuple[object, ...],
    env_info: dict[str, object],
) -> None:
    estimate = _trajectory_estimate_for_target(target)
    configured_platform = _configured_platform_for_plans(plans)
    openmm_platforms = env_info.get("openmm_platforms") or []
    typer.echo("Execution environment:")
    typer.echo(f"  hostname: {env_info.get('hostname')}")
    typer.echo(f"  cwd: {env_info.get('cwd')}")
    typer.echo(f"  interpretation: {env_info.get('interpretation')}")
    typer.echo(f"  configured OpenMM platform: {configured_platform or 'unknown'}")
    typer.echo(
        "  OpenMM platforms available: "
        + (", ".join(str(item) for item in openmm_platforms) if openmm_platforms else "unknown")
    )
    typer.echo(f"  estimated trajectory size: {_format_size_estimate(estimate)}")


def _enforce_execution_policy(target: Path, env_info: dict[str, object]) -> None:
    policy = _execution_policy_for_target(target)
    require_remote = bool(policy.get("require_remote_for_md", False))
    expected = policy.get("expected_hostname_contains")
    hostname_text = " ".join(
        str(env_info.get(key) or "") for key in ("hostname", "fqdn")
    )
    interpretation = str(env_info.get("interpretation") or "")
    appears_remote = "remote SSH host" in interpretation
    if require_remote and not appears_remote:
        raise ValueError(
            "execution.require_remote_for_md is true, but env-check did not detect "
            "SSH/remote execution."
        )
    if require_remote and expected and str(expected) not in hostname_text:
        raise ValueError(
            "execution.expected_hostname_contains does not match this host: "
            f"{expected!r}"
        )


def _confirm_local_nontrivial_run(target: Path, env_info: dict[str, object]) -> None:
    interpretation = str(env_info.get("interpretation") or "")
    if "running locally" not in interpretation:
        return
    estimate = _trajectory_estimate_for_target(target)
    total_gb = float(estimate.get("total_dcd_gb") or 0.0) if estimate else 0.0
    if total_gb < 1.0:
        return
    message = (
        "This appears to be running on the local machine. Continue? "
        f"Estimated DCD output is {total_gb:.3f} GB."
    )
    if not typer.confirm(message, default=False):
        raise typer.Exit(1)


def _execution_policy_for_target(target: Path) -> dict[str, object]:
    project_root = _project_root_for_target(target)
    if project_root is None:
        return {}
    lock = _read_yaml(project_root / "project.lock.yaml")
    if isinstance(lock.get("execution"), dict):
        return dict(lock["execution"])
    resolved = _read_json(project_root / "config_resolved.json")
    workflow = resolved.get("workflow") if isinstance(resolved, dict) else {}
    if isinstance(workflow, dict) and isinstance(workflow.get("execution"), dict):
        return dict(workflow["execution"])
    return {}


def _trajectory_estimate_for_target(target: Path) -> dict[str, object]:
    project_root = _project_root_for_target(target) or target
    storage = _read_json(project_root / "storage_estimate.json")
    return storage if isinstance(storage, dict) else {}


def _configured_platform_for_plans(plans: tuple[object, ...]) -> str | None:
    for plan in plans:
        run_dir = getattr(plan, "run_dir", None)
        if run_dir is None:
            continue
        config = _read_yaml(Path(run_dir) / "config.yaml")
        platform_name = config.get("platform")
        if platform_name:
            return str(platform_name)
    return None


def _project_root_for_target(target: Path) -> Path | None:
    path = target.resolve()
    candidates = [path, *path.parents]
    for candidate in candidates:
        if (candidate / "project.lock.yaml").is_file() or (candidate / "manifest.csv").is_file():
            return candidate
    return None


def _format_size_estimate(estimate: dict[str, object]) -> str:
    if not estimate:
        return "unknown"
    total_gb = estimate.get("total_dcd_gb")
    warning = estimate.get("warning_level")
    if total_gb is None:
        return "unknown"
    return f"{float(total_gb):.4f} GB DCD ({warning or 'unknown'} warning level)"


def _read_yaml(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _read_json(path: Path) -> dict[str, object]:
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


@app.command("run-recipe")
def run_recipe_command(
    recipe: Annotated[Path, typer.Argument(help="Python recipe file.")],
) -> None:
    """Execute a Python recipe; compile an ``experiment`` variable if present."""

    import runpy

    try:
        namespace = runpy.run_path(str(recipe))
        experiment = namespace.get("experiment") or namespace.get("exp")
        if experiment is not None and hasattr(experiment, "compile"):
            locked = experiment.compile(force=True)
            typer.echo(f"Compiled recipe to: {locked.lock_yaml}")
        else:
            typer.echo(f"Executed recipe: {recipe}")
    except Exception as exc:
        typer.echo(f"Recipe failed: {exc}", err=True)
        raise typer.Exit(1) from exc


@app.command("wizard")
def wizard_command() -> None:
    """Interactively generate a concise YAML config without running MD."""

    try:
        project_name = typer.prompt("Project name", default="my_project")
        source = typer.prompt("Protein source (direct/fasta/uniprot)", default="direct")
        protein: dict[str, object] = {"source": source, "charge_termini": "both"}
        if source == "direct":
            protein["name"] = typer.prompt("Protein name", default="ProteinA")
            protein["sequence"] = typer.prompt("Sequence")
        elif source == "fasta":
            protein["name"] = typer.prompt("Protein name", default="ProteinA")
            protein["fasta"] = typer.prompt("FASTA path")
        elif source == "uniprot":
            protein["query"] = typer.prompt("UniProt query")
            protein["organism"] = typer.prompt("Organism", default="Homo sapiens")
            protein["reviewed_only"] = typer.confirm("Reviewed Swiss-Prot only?", default=True)
            protein["interactive_select"] = True
        else:
            raise ValueError("Protein source must be direct, fasta, or uniprot.")
        experiment_type = typer.prompt(
            "Experiment type (smoke/production/ptm/cleavage/phase)",
            default="smoke",
        )
        preset = {
            "smoke": "smoke_single_chain",
            "production": "production_single_chain",
            "ptm": "short_single_chain",
            "cleavage": "cleavage_smoke",
            "phase": "phase_smoke",
        }.get(experiment_type, "smoke_single_chain")
        output = Path(typer.prompt("Output config path", default=f"configs/{project_name}.yaml"))
        config = {
            "project": {"name": project_name, "outdir": f"runs/{project_name}"},
            "input": {
                "protein": protein,
                "ptm": {"mode": "none"},
                "cleavage": {"mode": "none"},
            },
            "protocol": {"preset": preset},
            "analysis": {"preset": "standard_idr"},
            "report": {"preset": "standard"},
        }
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
    except Exception as exc:
        typer.echo(f"Wizard failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Wrote config: {output}")


@app.command("hpc-script")
def hpc_script_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Prepared project directory containing manifest.csv."),
    ],
    scheduler: Annotated[
        str,
        typer.Option("--scheduler", help="Scheduler script type. Currently only slurm."),
    ] = "slurm",
    phase: Annotated[
        str,
        typer.Option("--phase", help="Run phase: equilibration, production, or all."),
    ] = "all",
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Directory for generated scheduler files."),
    ] = None,
    python_executable: Annotated[
        str,
        typer.Option("--python", help="Python executable used by array tasks."),
    ] = "python",
    submit: Annotated[
        bool,
        typer.Option("--submit", help="Submit the generated script. Not implemented yet."),
    ] = False,
) -> None:
    """Generate an HPC scheduler script without submitting it."""

    from idrptm.hpc import write_slurm_array_script

    try:
        if scheduler != "slurm":
            raise ValueError("Only --scheduler slurm is supported.")
        if phase not in {"equilibration", "production", "all"}:
            raise ValueError("Phase must be 'equilibration', 'production', or 'all'.")
        if submit:
            raise ValueError("Submission is intentionally not automatic yet; run sbatch manually.")
        script = write_slurm_array_script(
            project_dir,
            output_dir=output_dir,
            phase=cast(RunPhase, phase),
            python_executable=python_executable,
        )
    except Exception as exc:
        typer.echo(f"HPC script generation failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Wrote SLURM array script: {script}")


@app.command("analyze")
def analyze_command(
    run_dir: Annotated[
        Path,
        typer.Argument(
            help="Prepared CALVADOS run directory containing top.pdb and trajectory.dcd."
        ),
    ],
    config: ConfigOption = None,
    topology: Annotated[
        Path | None,
        typer.Option("--topology", help="Topology PDB path. Defaults to RUN_DIR/top.pdb."),
    ] = None,
    trajectory: Annotated[
        Path | None,
        typer.Option(
            "--trajectory",
            "-t",
            help="Trajectory DCD path. Defaults to RUN_DIR/trajectory.dcd.",
        ),
    ] = None,
    output_dir: Annotated[
        Path | None,
        typer.Option(
            "--output-dir",
            "-o",
            help="Analysis output directory. Defaults to RUN_DIR/analysis.",
        ),
    ] = None,
    trajectory_reader: Annotated[
        str,
        typer.Option(
            "--trajectory-reader",
            help="Trajectory reader backend: mdtraj or mdanalysis.",
        ),
    ] = "mdtraj",
    force: Annotated[
        bool,
        typer.Option("--force", help="Recompute analysis even if cache is valid."),
    ] = False,
) -> None:
    """Analyze a prepared CALVADOS trajectory."""

    from idrptm.analysis.pipeline import analyze_run_directory

    try:
        if trajectory_reader not in {"mdtraj", "mdanalysis"}:
            raise ValueError("Trajectory reader must be 'mdtraj' or 'mdanalysis'.")
        result = analyze_run_directory(
            run_dir,
            config_path=config,
            topology=topology,
            trajectory=trajectory,
            trajectory_reader=cast(TrajectoryReader, trajectory_reader),
            output_dir=output_dir,
            force=force,
        )
    except Exception as exc:
        typer.echo(f"Analyze failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote analysis outputs to: {result.output_dir}")
    typer.echo(f"Summary: {result.summary_json}")


@app.command("status")
def status_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
) -> None:
    """Summarize project run status."""

    from idrptm.project import format_project_status, summarize_project_status

    typer.echo(format_project_status(summarize_project_status(project_dir)))


@app.command("resume")
def resume_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
    phase: Annotated[str, typer.Option("--phase", help="Run phase.")] = "all",
    force: Annotated[bool, typer.Option("--force", help="Rerun completed runs too.")] = False,
    python_executable: Annotated[
        str,
        typer.Option("--python", help="Python executable used to invoke generated run scripts."),
    ] = "python",
    dry_run: Annotated[
        bool,
        typer.Option("--dry-run/--execute", help="Plan resume commands or execute them."),
    ] = True,
) -> None:
    """Resume failed or incomplete runs."""

    from idrptm.project import resume_project

    try:
        if phase not in {"equilibration", "production", "all"}:
            raise ValueError("Phase must be 'equilibration', 'production', or 'all'.")
        results = resume_project(
            project_dir,
            phase=cast(RunPhase, phase),
            force=force,
            dry_run=dry_run,
            python_executable=python_executable,
        )
    except Exception as exc:
        typer.echo(f"Resume failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    for result in results:
        command = getattr(result, "command", None)
        run_dir = getattr(result, "run_dir", None)
        status = getattr(result, "status", "planned")
        typer.echo(f"{run_dir}: {status} {command or ''}")


@app.command("clean")
def clean_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
    yes: Annotated[bool, typer.Option("--yes", "-y", help="Confirm safe cleanup.")] = False,
) -> None:
    """Safely remove cache/temp files, never raw trajectories."""

    from idrptm.project import clean_project

    removed = clean_project(project_dir, yes=yes)
    if not yes:
        typer.echo("Dry run: pass --yes to remove cache/temp files.")
        return
    typer.echo(f"Removed {len(removed)} file(s).")


@app.command("list")
def list_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
) -> None:
    """List runs from the optional SQLite registry."""

    from idrptm.registry import list_runs

    for row in list_runs(project_dir):
        typer.echo(
            ",".join(
                str(row.get(key, ""))
                for key in ("run_id", "ptm_state", "status", "run_dir")
            )
        )


@app.command("summary")
def summary_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
) -> None:
    """Show registry run counts by status."""

    from idrptm.registry import summarize_registry

    for status, count in sorted(summarize_registry(project_dir).items()):
        typer.echo(f"{status},{count}")


@app.command("query")
def query_command(
    project_dir: Annotated[Path, typer.Argument(help="Compiled/prepared project directory.")],
    expression: Annotated[str, typer.Argument(help="Pandas-style query expression.")],
) -> None:
    """Query the run registry with a simple expression."""

    from idrptm.registry import query_runs

    for row in query_runs(project_dir, expression):
        typer.echo(
            ",".join(
                str(row.get(key, ""))
                for key in ("run_id", "ptm_state", "status", "run_dir")
            )
        )


@app.command("compare")
def compare_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Project directory containing manifest.csv and analyzed runs."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Comparison output directory."),
    ] = None,
) -> None:
    """Compare each PTM condition against the detected WT condition."""

    from idrptm.analysis.compare import compare_project

    try:
        result = compare_project(project_dir, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Compare failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"WT condition: {result.wt_condition}")
    typer.echo(f"Wrote comparison outputs to: {result.output_dir}")
    typer.echo(f"Summary: {result.outputs['summary_csv']}")


@app.command("report")
def report_command(
    project_dir: Annotated[
        Path,
        typer.Argument(help="Project directory containing manifest.csv and analyzed runs."),
    ],
    output_dir: Annotated[
        Path | None,
        typer.Option("--output-dir", "-o", help="Report output directory."),
    ] = None,
) -> None:
    """Generate WT-vs-PTM comparison figures and a Markdown report."""

    from idrptm.plotting.report import generate_report

    try:
        result = generate_report(project_dir, output_dir=output_dir)
    except Exception as exc:
        typer.echo(f"Report failed: {exc}", err=True)
        raise typer.Exit(1) from exc

    typer.echo(f"Wrote report: {result.report_path}")
    typer.echo(f"Wrote {len(result.figure_paths)} figure file(s).")


@app.command("repo-check")
def repo_check_command(
    output: Annotated[
        Path,
        typer.Option("--output", "-o", help="Markdown output path."),
    ] = Path("REPO_CHECK.md"),
) -> None:
    """Run practical public-repository readiness checks."""

    from protein_analysis_md.repo_check import run_repo_check

    try:
        result = run_repo_check(output_path=output)
    except Exception as exc:
        typer.echo(f"Repo check failed: {exc}", err=True)
        raise typer.Exit(1) from exc
    typer.echo(f"Wrote repo check: {result.output_path}")
    readme = result.findings["readme"]
    typer.echo(f"README lines: {readme['line_count']}")
    generated = len(result.findings["generated_files"])
    tracked_generated = len(result.findings["tracked_generated"])
    typer.echo(f"Generated/output file warnings: {generated}")
    typer.echo(f"Tracked generated/output file warnings: {tracked_generated}")


def main() -> None:
    """Console-script entry point."""

    app()


if __name__ == "__main__":
    main()
