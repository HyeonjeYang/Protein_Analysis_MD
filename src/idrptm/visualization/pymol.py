"""Organize prepared/completed runs for later PyMOL visualization."""

from __future__ import annotations

import csv
import json
import os
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from idrptm.provenance import slugify

ExportMode = Literal["symlink", "copy"]


@dataclass(frozen=True)
class PyMOLRunExport:
    """Files written for one run's PyMOL-friendly view."""

    run_id: str
    source_run_dir: Path
    output_dir: Path
    topology: Path | None
    trajectory: Path | None
    script: Path
    metadata: Path | None
    parameters: Path | None
    status: str
    missing: tuple[str, ...]


@dataclass(frozen=True)
class PyMOLExportResult:
    """Result of exporting a project for PyMOL visualization."""

    output_dir: Path
    manifest_csv: Path
    load_all_pml: Path
    readme: Path
    runs: tuple[PyMOLRunExport, ...]


def export_pymol_project(
    project_dir: str | Path,
    *,
    output_dir: str | Path | None = None,
    mode: ExportMode = "symlink",
    include_missing: bool = False,
    force: bool = False,
) -> PyMOLExportResult:
    """Create a PyMOL-friendly directory with linked or copied trajectory assets."""

    root = Path(project_dir).resolve()
    out = Path(output_dir).resolve() if output_dir is not None else root / "pymol"
    out.mkdir(parents=True, exist_ok=True)
    run_exports: list[PyMOLRunExport] = []
    for row in _manifest_rows(root):
        run_id = row.get("variant_id") or row.get("run_id") or "run"
        run_dir = _run_dir_from_row(root, row)
        topology = _find_topology(run_dir)
        trajectory = _find_trajectory(run_dir, row)
        missing = tuple(
            label
            for label, path in (("topology", topology), ("trajectory", trajectory))
            if path is None
        )
        if missing and not include_missing:
            continue
        run_out = out / slugify(run_id)
        run_out.mkdir(parents=True, exist_ok=True)
        exported_topology = (
            _place_asset(topology, run_out / "topology.pdb", mode=mode, force=force)
            if topology is not None
            else None
        )
        exported_trajectory = (
            _place_asset(trajectory, run_out / "trajectory.dcd", mode=mode, force=force)
            if trajectory is not None
            else None
        )
        metadata = _optional_asset(run_dir / "metadata.json", run_out, mode=mode, force=force)
        parameters = _optional_asset(run_dir / "parameters.txt", run_out, mode=mode, force=force)
        script = run_out / "load.pml"
        script.write_text(
            _run_pml(
                object_name=slugify(run_id),
                has_topology=exported_topology is not None,
                has_trajectory=exported_trajectory is not None,
                metadata_path=metadata,
            ),
            encoding="utf-8",
        )
        run_exports.append(
            PyMOLRunExport(
                run_id=run_id,
                source_run_dir=run_dir,
                output_dir=run_out,
                topology=exported_topology,
                trajectory=exported_trajectory,
                script=script,
                metadata=metadata,
                parameters=parameters,
                status="ready" if not missing else "missing_" + "_".join(missing),
                missing=missing,
            )
        )
    manifest = out / "pymol_manifest.csv"
    _write_manifest(manifest, run_exports)
    load_all = out / "load_all.pml"
    load_all.write_text(_load_all_pml(out, run_exports), encoding="utf-8")
    readme = out / "README.md"
    readme.write_text(_readme(out, mode, run_exports), encoding="utf-8")
    return PyMOLExportResult(
        output_dir=out,
        manifest_csv=manifest,
        load_all_pml=load_all,
        readme=readme,
        runs=tuple(run_exports),
    )


def _manifest_rows(root: Path) -> list[dict[str, str]]:
    manifest = root / "manifest.csv"
    if manifest.exists():
        with manifest.open(newline="", encoding="utf-8") as handle:
            return list(csv.DictReader(handle))
    runs_dir = root / "runs"
    if not runs_dir.exists():
        raise FileNotFoundError(f"Could not find manifest.csv or runs/ under {root}.")
    return [
        {
            "variant_id": path.name,
            "metadata_path": (path / "metadata.json").relative_to(root).as_posix(),
        }
        for path in sorted(runs_dir.iterdir())
        if path.is_dir()
    ]


def _run_dir_from_row(root: Path, row: dict[str, str]) -> Path:
    metadata_path = row.get("metadata_path")
    if metadata_path:
        return (root / metadata_path).parent
    run_dir = row.get("run_dir")
    if run_dir:
        return root / run_dir
    return root / "runs" / (row.get("variant_id") or row.get("run_id") or "")


def _find_topology(run_dir: Path) -> Path | None:
    path = run_dir / "top.pdb"
    return path if path.exists() else None


def _find_trajectory(run_dir: Path, row: dict[str, str]) -> Path | None:
    candidates = []
    for key in ("variant_id", "run_id"):
        if row.get(key):
            candidates.append(run_dir / f"{row[key]}.dcd")
    candidates.append(run_dir / "trajectory.dcd")
    candidates.extend(sorted(run_dir.glob("*.dcd")))
    for path in candidates:
        if path.exists():
            return path
    return None


def _place_asset(source: Path, destination: Path, *, mode: ExportMode, force: bool) -> Path:
    if destination.exists() or destination.is_symlink():
        if not force:
            return destination
        destination.unlink()
    if mode == "copy":
        shutil.copy2(source, destination)
    else:
        destination.symlink_to(_relative_target(source, destination.parent))
    return destination


def _optional_asset(
    source: Path,
    destination_dir: Path,
    *,
    mode: ExportMode,
    force: bool,
) -> Path | None:
    if not source.exists():
        return None
    return _place_asset(source, destination_dir / source.name, mode=mode, force=force)


def _relative_target(source: Path, start: Path) -> Path:
    return Path(os.path.relpath(source.resolve(), start.resolve()))


def _run_pml(
    *,
    object_name: str,
    has_topology: bool,
    has_trajectory: bool,
    metadata_path: Path | None,
) -> str:
    lines = [
        "# Generated by protein_analysis_md. Run with: pymol load.pml",
        "reinitialize",
    ]
    lines.extend(
        _object_pml_lines(
            object_name=object_name,
            topology_path=Path("topology.pdb") if has_topology else None,
            trajectory_path=Path("trajectory.dcd") if has_trajectory else None,
            metadata_path=metadata_path,
        )
    )
    lines.extend(
        [
            "set orthoscopic, on",
            "set ray_opaque_background, off",
            "bg_color white",
            "orient",
        ]
    )
    return "\n".join(lines) + "\n"


def _object_pml_lines(
    *,
    object_name: str,
    topology_path: Path | None,
    trajectory_path: Path | None,
    metadata_path: Path | None,
) -> list[str]:
    lines: list[str] = []
    if topology_path is not None:
        lines.append(f"load {topology_path.as_posix()}, {object_name}")
    if trajectory_path is not None:
        lines.append(f"load_traj {trajectory_path.as_posix()}, {object_name}")
    lines.extend(
        [
            f"hide everything, {object_name}",
            f"show spheres, {object_name}",
            f"set sphere_scale, 0.35, {object_name}",
            f"spectrum count, rainbow, {object_name}",
        ]
    )
    ptm_sites = _ptm_sites(metadata_path)
    if ptm_sites:
        selector = "+".join(str(site) for site in ptm_sites)
        selection_name = f"{object_name}_ptm_sites"
        lines.extend(
            [
                f"select {selection_name}, {object_name} and resi {selector}",
                f"color yelloworange, {selection_name}",
                f"show spheres, {selection_name}",
            ]
        )
    return lines


def _ptm_sites(metadata_path: Path | None) -> list[int]:
    if metadata_path is None or not metadata_path.exists():
        return []
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except Exception:
        return []
    raw = metadata.get("ptm_sites_1based", "")
    if isinstance(raw, list):
        return [int(item) for item in raw if str(item).isdigit()]
    if isinstance(raw, str):
        return [int(item) for item in raw.replace(";", ",").split(",") if item.strip().isdigit()]
    return []


def _write_manifest(path: Path, runs: list[PyMOLRunExport]) -> None:
    fieldnames = [
        "run_id",
        "status",
        "source_run_dir",
        "pymol_dir",
        "topology",
        "trajectory",
        "script",
        "metadata",
        "parameters",
        "missing",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
        writer.writeheader()
        for item in runs:
            writer.writerow(
                {
                    "run_id": item.run_id,
                    "status": item.status,
                    "source_run_dir": item.source_run_dir,
                    "pymol_dir": item.output_dir,
                    "topology": item.topology or "",
                    "trajectory": item.trajectory or "",
                    "script": item.script,
                    "metadata": item.metadata or "",
                    "parameters": item.parameters or "",
                    "missing": ",".join(item.missing),
                }
            )


def _load_all_pml(root: Path, runs: list[PyMOLRunExport]) -> str:
    lines = [
        "# Generated by protein_analysis_md. Run with: pymol load_all.pml",
        "reinitialize",
    ]
    for item in runs:
        if item.topology is None:
            continue
        lines.extend(
            _object_pml_lines(
                object_name=slugify(item.run_id),
                topology_path=item.topology.relative_to(root),
                trajectory_path=item.trajectory.relative_to(root) if item.trajectory else None,
                metadata_path=item.metadata,
            )
        )
    lines.extend(
        [
            "set orthoscopic, on",
            "set ray_opaque_background, off",
            "bg_color white",
            "orient",
        ]
    )
    return "\n".join(lines) + "\n"


def _readme(root: Path, mode: ExportMode, runs: list[PyMOLRunExport]) -> str:
    ready = sum(1 for item in runs if item.status == "ready")
    return "\n".join(
        [
            "# PyMOL Export",
            "",
            "Generated by `protein_analysis_md`.",
            "",
            f"- Export directory: `{root}`",
            f"- Asset mode: `{mode}`",
            f"- Ready runs: {ready}/{len(runs)}",
            "",
            "Open one run:",
            "",
            "```bash",
            "pymol <RUN_ID>/load.pml",
            "```",
            "",
            "Open all exported runs:",
            "",
            "```bash",
            "pymol load_all.pml",
            "```",
            "",
            "By default, large trajectory files are symlinked rather than copied.",
            "Use `pamd pymol --mode copy` when you need a portable folder.",
            "",
        ]
    )
