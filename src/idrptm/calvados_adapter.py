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

        config_yaml = _calvados_config(config, row["variant_id"])
        components_yaml = _calvados_components(config, row["variant_id"])
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


def _calvados_config(config: WorkflowConfig, variant_id: str) -> dict[str, object]:
    steps = config.calvados.nframes * config.calvados.wfreq
    return {
        "sysname": variant_id,
        "box": config.calvados.box_nm,
        "temp": config.calvados.temperature_k,
        "ionic": config.calvados.ionic_m,
        "pH": config.calvados.ph,
        "topol": config.calvados.topol,
        "wfreq": config.calvados.wfreq,
        "steps": steps,
        "runtime": config.calvados.runtime,
        "platform": config.calvados.platform,
        "restart": config.calvados.restart,
        "frestart": config.calvados.frestart,
        "verbose": config.calvados.verbose,
    }


def _calvados_components(config: WorkflowConfig, variant_id: str) -> dict[str, object]:
    return {
        "defaults": {
            "molecule_type": config.calvados.molecule_type,
            "nmol": config.calvados.nmol,
            "restraint": False,
            "charge_termini": config.calvados.charge_termini,
            "fresidues": "residues.csv",
            "ffasta": "input.fasta",
        },
        "system": {
            variant_id: {},
        },
    }


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
        "ptm_state": manifest_row["ptm_state"],
        "ptm_count": int(manifest_row["ptm_count"]),
        "ptm_sites_1based": manifest_row["ptm_sites_1based"],
        "ptm_sites_0based": manifest_row["ptm_sites_0based"],
        "original_sequence": manifest_row["original_sequence"],
        "simulation_sequence": manifest_row["simulation_sequence"],
        "calvados": {
            "model": config.calvados.model,
            "config_yaml": paths.config_yaml.name,
            "components_yaml": paths.components_yaml.name,
            "run_script": paths.run_script.name,
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
