"""Prepare CALVADOS-compatible run directories without mutating upstream files."""

from __future__ import annotations

import csv
import json
import shutil
from dataclasses import dataclass
from pathlib import Path

import yaml

from idrptm.design import design_variants, write_design_outputs
from idrptm.residue_params import write_residue_parameters
from idrptm.schema import WorkflowConfig, load_config


@dataclass(frozen=True)
class CalvadosRunDirectory:
    """Paths that identify a prepared CALVADOS run directory."""

    path: Path
    input_fasta: Path
    residues_csv: Path
    config_yaml: Path
    components_yaml: Path
    run_script: Path
    metadata_json: Path


@dataclass(frozen=True)
class PreparationResult:
    """Result of preparing or planning CALVADOS run directories."""

    output_dir: Path
    manifest_path: Path
    run_directories: tuple[CalvadosRunDirectory, ...]
    dry_run: bool


def prepare_calvados_run_directory(output_dir: str | Path) -> CalvadosRunDirectory:
    """Return canonical paths for one CALVADOS run directory."""

    run_dir = Path(output_dir)
    return CalvadosRunDirectory(
        path=run_dir,
        input_fasta=run_dir / "input.fasta",
        residues_csv=run_dir / "residues.csv",
        config_yaml=run_dir / "config.yaml",
        components_yaml=run_dir / "components.yaml",
        run_script=run_dir / "run.py",
        metadata_json=run_dir / "metadata.json",
    )


def prepare_from_config_file(
    config_path: str | Path,
    output_dir: str | Path | None = None,
    dry_run: bool = False,
) -> PreparationResult:
    """Load a workflow config and prepare CALVADOS run directories."""

    config = load_config(config_path)
    return prepare_from_config(config=config, output_dir=output_dir, dry_run=dry_run)


def prepare_from_config(
    config: WorkflowConfig,
    output_dir: str | Path | None = None,
    dry_run: bool = False,
) -> PreparationResult:
    """Prepare one CALVADOS run directory per manifest row."""

    root = Path(output_dir) if output_dir is not None else config.runner.work_dir
    manifest_path = root / "manifest.csv"
    if dry_run:
        run_directories = tuple(
            prepare_calvados_run_directory(root / "runs" / variant.variant_id)
            for variant in design_variants(config)
        )
        return PreparationResult(
            output_dir=root,
            manifest_path=manifest_path,
            run_directories=run_directories,
            dry_run=True,
        )

    design_output = write_design_outputs(config, output_dir=root)
    manifest_rows = _read_manifest(design_output.manifest_path)
    prepared: list[CalvadosRunDirectory] = []
    for row in manifest_rows:
        run_dir = (root / row["metadata_path"]).parent
        paths = prepare_calvados_run_directory(run_dir)
        paths.path.mkdir(parents=True, exist_ok=True)

        source_fasta = root / row["fasta_path"]
        shutil.copyfile(source_fasta, paths.input_fasta)
        residue_result = write_residue_parameters(
            source_csv=config.calvados.residue_parameters,
            output_csv=paths.residues_csv,
            ph=config.calvados.ph,
        )

        config_yaml = _calvados_config(config, row)
        components_yaml = _calvados_components(config, row)
        paths.config_yaml.write_text(
            yaml.safe_dump(config_yaml, sort_keys=False),
            encoding="utf-8",
        )
        paths.components_yaml.write_text(
            yaml.safe_dump(components_yaml, sort_keys=False),
            encoding="utf-8",
        )
        paths.run_script.write_text(_run_script_text(), encoding="utf-8")
        metadata = _metadata_json(config, row, paths, residue_result.metadata())
        paths.metadata_json.write_text(
            json.dumps(metadata, indent=2) + "\n",
            encoding="utf-8",
        )
        prepared.append(paths)

    return PreparationResult(
        output_dir=root,
        manifest_path=design_output.manifest_path,
        run_directories=tuple(prepared),
        dry_run=False,
    )


def _read_manifest(manifest_path: Path) -> list[dict[str, str]]:
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _calvados_config(config: WorkflowConfig, manifest_row: dict[str, str]) -> dict[str, object]:
    simulation = config.calvados.simulation
    variant_id = manifest_row["variant_id"]
    random_seed = _manifest_random_seed(manifest_row, simulation.random_seed)
    topol = config.calvados.topol
    if _is_multi_component_row(manifest_row) and topol == "center":
        topol = manifest_row.get("placement") or "grid"
    return {
        "sysname": variant_id,
        "box": config.calvados.box_nm,
        "temp": config.calvados.temperature_k,
        "ionic": config.calvados.ionic_m,
        "pH": config.calvados.ph,
        "topol": topol,
        "fixed_lambda": 0,
        "eps_lj": 0.2,
        "cutoff_lj": 2.0,
        "cutoff_yu": 4.0,
        "wfreq": simulation.save_every_steps,
        "steps": simulation.total_steps,
        "runtime": simulation.runtime_hours,
        "platform": simulation.platform,
        "threads": 1,
        "restart": simulation.restart,
        "frestart": simulation.checkpoint_file,
        "verbose": config.calvados.verbose,
        "slab_eq": False,
        "bilayer_eq": False,
        "pressure_coupling": False,
        "box_eq": False,
        "pressure": [0, 0, 0],
        "boxscaling_xyz": [True, True, True],
        "k_eq": 0.02,
        "steps_eq": 1000,
        "ext_force": False,
        "ext_force_expr": "step(d2-18)*d2; d2=periodicdistance(x, y, z, 0, 0, z)^2",
        "friction_coeff": 0.01,
        "slab_width": 100,
        "slab_outer": 40,
        "random_number_seed": random_seed,
        "report_potential_energy": False,
        "logfreq": 1000000,
        "gpu_id": 0,
        "custom_restraints": False,
        "custom_restraint_type": "harmonic",
        "fcustom_restraints": "custom_restraints.txt",
    }


def _calvados_components(config: WorkflowConfig, manifest_row: dict[str, str]) -> dict[str, object]:
    components = _manifest_components(manifest_row, config)
    return {
        "defaults": {
            "molecule_type": config.calvados.molecule_type,
            "nmol": config.calvados.nmol,
            "restraint": False,
            "charge_termini": config.calvados.charge_termini,
            "fresidues": "residues.csv",
            "ffasta": "input.fasta",
            "alpha": 0,
            "kb": 8033.0,
            "ext_restraint": True,
            "cutoff_restr": 0.9,
            "pdb_folder": "pdbs",
            "restraint_type": "harmonic",
            "k_harmonic": 700.0,
            "fdomains": "domains.yaml",
            "k_go": 15.0,
            "use_com": True,
            "periodic": False,
            "colabfold": 0,
            "bfac_shift": 0.8,
            "bfac_width": 50.0,
            "pae_shift": 0.3,
            "pae_width": 15.0,
            "rna_kb1": 8033.0,
            "rna_kb2": 8033.0,
            "rna_ka": 7.24,
            "rna_pa": 3.14,
            "rna_nb_sigma": 0.4,
            "rna_nb_scale": 15,
            "rna_nb_cutoff": 0.6,
            "n_ends": 1,
            "ptm_name": "protein_analysis_md_ptm",
            "ptm_locations": [],
        },
        "system": {
            component["component_name"]: {
                "molecule_type": component.get("molecule_type") or config.calvados.molecule_type,
                "nmol": int(component.get("copies", config.calvados.nmol)),
            }
            for component in components
        },
    }


def _manifest_components(
    manifest_row: dict[str, str],
    config: WorkflowConfig,
) -> list[dict[str, object]]:
    components_json = manifest_row.get("components_json")
    if components_json:
        parsed = json.loads(components_json)
        if not isinstance(parsed, list):
            raise ValueError("Manifest components_json must decode to a list.")
        return [dict(component) for component in parsed]
    return [
        {
            "component_name": manifest_row["variant_id"],
            "protein_name": manifest_row["base_sequence_name"],
            "protein_variant_id": manifest_row["variant_id"],
            "ptm_state": manifest_row["ptm_state"],
            "copies": 1,
            "molecule_type": config.calvados.molecule_type,
            "original_sequence": manifest_row["original_sequence"],
            "simulation_sequence": manifest_row["simulation_sequence"],
        }
    ]


def _is_multi_component_row(manifest_row: dict[str, str]) -> bool:
    return (
        manifest_row.get("is_multi_component") == "1"
        or int(manifest_row.get("component_count", "1")) > 1
        or int(manifest_row.get("total_molecule_copies", "1")) > 1
    )


def _run_script_text() -> str:
    return '''"""Run a prepared CALVADOS simulation directory."""

from argparse import ArgumentParser

from calvados import sim


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("--path", nargs="?", default=".", const=".", type=str)
    parser.add_argument("--config", nargs="?", default="config.yaml", const="config.yaml", type=str)
    parser.add_argument(
        "--components",
        nargs="?",
        default="components.yaml",
        const="components.yaml",
        type=str,
    )
    args = parser.parse_args()
    sim.run(path=args.path, fconfig=args.config, fcomponents=args.components)
'''


def _metadata_json(
    config: WorkflowConfig,
    manifest_row: dict[str, str],
    paths: CalvadosRunDirectory,
    residue_metadata: dict[str, object],
) -> dict[str, object]:
    return {
        "schema_version": 1,
        "project": config.project,
        "variant_id": manifest_row["variant_id"],
        "base_variant_id": manifest_row.get("base_variant_id") or manifest_row["variant_id"],
        "replicate": _manifest_int(manifest_row.get("replicate"), 1),
        "replicate_id": manifest_row.get("replicate_id") or "rep001",
        "random_seed": _manifest_random_seed(
            manifest_row,
            config.calvados.simulation.random_seed,
        ),
        "system_name": manifest_row.get("system_name", manifest_row["variant_id"]),
        "placement": manifest_row.get("placement", config.calvados.topol),
        "is_multi_component": _is_multi_component_row(manifest_row),
        "ptm_state": manifest_row["ptm_state"],
        "ptm_count": int(manifest_row["ptm_count"]),
        "ptm_sites_1based": manifest_row["ptm_sites_1based"],
        "ptm_sites_0based": manifest_row["ptm_sites_0based"],
        "original_sequence": manifest_row["original_sequence"],
        "simulation_sequence": manifest_row["simulation_sequence"],
        "components": _manifest_components(manifest_row, config),
        "calvados": {
            "model": config.calvados.model,
            "config_yaml": paths.config_yaml.name,
            "components_yaml": paths.components_yaml.name,
            "run_script": paths.run_script.name,
        },
        "simulation": config.calvados.simulation.metadata()
        | {
            "random_seed": _manifest_random_seed(
                manifest_row,
                config.calvados.simulation.random_seed,
            )
        },
        "residue_parameters": residue_metadata,
        "files": {
            "input_fasta": paths.input_fasta.name,
            "residues_csv": paths.residues_csv.name,
            "metadata_json": paths.metadata_json.name,
        },
        "upstream_policy": {
            "mutates_calvados_source": False,
        },
    }


def _manifest_int(value: str | None, default: int) -> int:
    if value in (None, ""):
        return default
    return int(value)


def _manifest_random_seed(
    manifest_row: dict[str, str],
    default: int | None,
) -> int | None:
    value = manifest_row.get("random_seed")
    if value in (None, ""):
        return default
    return int(value)
