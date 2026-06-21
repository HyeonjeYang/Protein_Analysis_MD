"""Finalize projects after simulations complete."""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from idrptm.schema import load_config


@dataclass(frozen=True)
class FinalizeResult:
    """Summary of finalize outputs."""

    project_dir: Path
    analyzed_runs: tuple[Path, ...]
    skipped_analysis: tuple[Path, ...]
    comparison_dir: Path | None
    report_path: Path | None
    pymol_dir: Path | None
    visualization: bool


def finalize_project(
    project_dir: str | Path,
    *,
    analysis_parallel: int = 1,
    force_analysis: bool = False,
    visualization: bool | None = None,
    pymol: bool = True,
) -> FinalizeResult:
    """Run analysis/comparison and optional visualization outputs for a project."""

    root = Path(project_dir)
    config_path = root / "project.lock.yaml"
    config = load_config(config_path) if config_path.exists() else None
    visualization_enabled = (
        bool(config.visualization) if visualization is None and config is not None else True
    )
    if visualization is not None:
        visualization_enabled = visualization
    run_dirs = _run_dirs(root)
    to_analyze = tuple(
        run_dir
        for run_dir in run_dirs
        if force_analysis or not (run_dir / "analysis" / "summary.json").exists()
    )
    skipped = tuple(run_dir for run_dir in run_dirs if run_dir not in set(to_analyze))
    analyzed = _analyze_runs(
        to_analyze,
        config_path=config_path if config_path.exists() else None,
        max_workers=analysis_parallel,
        force=force_analysis,
    )

    from idrptm.analysis.compare import compare_project

    comparison = compare_project(root)
    report_path = None
    pymol_dir = None
    if visualization_enabled:
        from idrptm.plotting.report import generate_report

        report = generate_report(root)
        report_path = report.report_path
        if pymol:
            from idrptm.visualization.pymol import export_pymol_project

            pymol_result = export_pymol_project(root, include_missing=False, force=True)
            pymol_dir = pymol_result.output_dir

    return FinalizeResult(
        project_dir=root,
        analyzed_runs=tuple(analyzed),
        skipped_analysis=skipped,
        comparison_dir=comparison.output_dir,
        report_path=report_path,
        pymol_dir=pymol_dir,
        visualization=visualization_enabled,
    )


def _analyze_runs(
    run_dirs: tuple[Path, ...],
    *,
    config_path: Path | None,
    max_workers: int,
    force: bool,
) -> list[Path]:
    if not run_dirs:
        return []
    from idrptm.analysis.pipeline import analyze_run_directory

    workers = max(1, int(max_workers))
    analyzed: list[Path] = []
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {
            executor.submit(
                analyze_run_directory,
                run_dir,
                config_path=config_path,
                force=force,
            ): run_dir
            for run_dir in run_dirs
        }
        for future in as_completed(futures):
            run_dir = futures[future]
            future.result()
            analyzed.append(run_dir)
    return sorted(analyzed)


def _run_dirs(root: Path) -> tuple[Path, ...]:
    manifest = root / "manifest.csv"
    if manifest.exists():
        dirs: list[Path] = []
        with manifest.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("metadata_path"):
                    dirs.append((root / row["metadata_path"]).parent)
        return tuple(sorted(dict.fromkeys(dirs)))
    runs_root = root / "runs"
    if runs_root.exists():
        return tuple(sorted(path for path in runs_root.iterdir() if path.is_dir()))
    return ()
