"""Create shareable project bundles."""

from __future__ import annotations

import tarfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class BundleResult:
    """Path and file count for a project bundle."""

    archive_path: Path
    included_files: int
    include_trajectories: bool


def pack_project(
    project_dir: str | Path,
    *,
    output: str | Path | None = None,
    include_trajectories: bool = False,
) -> BundleResult:
    """Create a tar.gz bundle of shareable project outputs."""

    root = Path(project_dir).resolve()
    archive = (
        Path(output).resolve()
        if output is not None
        else root / f"{root.name}_bundle.tar.gz"
    )
    archive.parent.mkdir(parents=True, exist_ok=True)
    files = tuple(_iter_bundle_files(root, include_trajectories=include_trajectories))
    with tarfile.open(archive, "w:gz") as tar:
        for path in files:
            tar.add(path, arcname=Path(root.name) / path.relative_to(root), recursive=False)
    return BundleResult(
        archive_path=archive,
        included_files=len(files),
        include_trajectories=include_trajectories,
    )


def _iter_bundle_files(root: Path, *, include_trajectories: bool) -> list[Path]:
    allowed_names = {
        "project.lock.yaml",
        "config_resolved.json",
        "storage_estimate.json",
        "manifest.csv",
        "manifest_preview.csv",
        "pymol_manifest.csv",
        "load_all.pml",
        "README.md",
        "report.md",
        "comparison_summary.csv",
        "comparison_metadata.json",
        "parameters.txt",
        "metadata.json",
        "metadata.yaml",
        "config.yaml",
        "components.yaml",
        "input.fasta",
        "residues.csv",
        "run_status.json",
        "summary.json",
    }
    allowed_suffixes = {
        ".png",
        ".pdf",
        ".csv",
        ".json",
        ".yaml",
        ".yml",
        ".txt",
        ".md",
        ".pml",
        ".parquet",
        ".npy",
    }
    trajectory_suffixes = {".dcd", ".xtc", ".trr", ".nc", ".h5", ".chk", ".pdb", ".xml"}
    files: list[Path] = []
    for path in root.rglob("*"):
        if not path.is_file() or path == root:
            continue
        if path.name == root.name + "_bundle.tar.gz":
            continue
        if path.suffix in trajectory_suffixes and not include_trajectories:
            continue
        if (
            path.name in allowed_names
            or path.suffix in allowed_suffixes
            or (include_trajectories and path.suffix in trajectory_suffixes)
        ):
            if path.suffix in trajectory_suffixes and not include_trajectories:
                continue
            files.append(path)
    return sorted(dict.fromkeys(files))
