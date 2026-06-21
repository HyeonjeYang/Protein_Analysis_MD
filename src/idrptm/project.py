"""Project status, resume, and safe cleanup helpers."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from idrptm.runner import RunPhase, execute_local_runs, plan_local_runs


@dataclass(frozen=True)
class ProjectStatus:
    """Summary of run states in a project directory."""

    project_dir: Path
    counts: dict[str, int]
    run_statuses: tuple[dict[str, object], ...]


def summarize_project_status(project_dir: str | Path) -> ProjectStatus:
    """Summarize prepared/running/completed/failed/analyzed/reported states."""

    root = Path(project_dir)
    statuses = tuple(_status_for_run(run_dir) for run_dir in _run_dirs(root))
    counts: dict[str, int] = {}
    for status in statuses:
        for state in _states_for_status(status):
            counts[state] = counts.get(state, 0) + 1
    if (root / "report" / "report.md").exists():
        counts["reported"] = counts.get("reported", 0) + 1
    return ProjectStatus(project_dir=root, counts=counts, run_statuses=statuses)


def resume_project(
    project_dir: str | Path,
    *,
    phase: RunPhase = "all",
    force: bool = False,
    dry_run: bool = True,
    python_executable: str = "python",
) -> tuple[object, ...]:
    """Resume failed or incomplete local runs, skipping completed runs by default."""

    root = Path(project_dir)
    candidate_dirs = []
    for status in summarize_project_status(root).run_statuses:
        if force or status.get("status") != "completed":
            candidate_dirs.append(Path(str(status["run_dir"])))
    plans = tuple(
        plan
        for run_dir in candidate_dirs
        for plan in plan_local_runs(
            run_dir,
            phase=phase,
            all_runs=False,
            python_executable=python_executable,
        )
    )
    if dry_run:
        return plans
    return execute_local_runs(plans)


def clean_project(project_dir: str | Path, *, yes: bool = False) -> tuple[Path, ...]:
    """Safely remove derived analysis cache files and temporary scheduler lists."""

    if not yes:
        return ()
    root = Path(project_dir)
    removed: list[Path] = []
    for path in root.rglob("cache_manifest.json"):
        path.unlink(missing_ok=True)
        removed.append(path)
    for path in (root / "run_dirs.txt", root / "run_slurm_array.sh"):
        if path.exists():
            path.unlink()
            removed.append(path)
    return tuple(removed)


def format_project_status(status: ProjectStatus) -> str:
    """Format a project status summary for CLI output."""

    lines = [f"project_dir: {status.project_dir}", "status,count"]
    for key in sorted(status.counts):
        lines.append(f"{key},{status.counts[key]}")
    return "\n".join(lines)


def _run_dirs(root: Path) -> tuple[Path, ...]:
    manifest = root / "manifest.csv"
    if manifest.exists():
        import csv

        dirs = []
        with manifest.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle):
                if row.get("metadata_path"):
                    dirs.append((root / row["metadata_path"]).parent)
        return tuple(sorted(dict.fromkeys(dirs)))
    runs_root = root / "runs"
    if runs_root.exists():
        return tuple(sorted(path for path in runs_root.iterdir() if path.is_dir()))
    return ()


def _status_for_run(run_dir: Path) -> dict[str, object]:
    status_path = run_dir / "run_status.json"
    if status_path.exists():
        payload = json.loads(status_path.read_text(encoding="utf-8"))
    else:
        payload = {"status": "prepared" if (run_dir / "run.py").exists() else "unknown"}
    payload["run_dir"] = str(run_dir)
    payload["has_checkpoint"] = any(
        (run_dir / name).exists() for name in ("restart.chk", "checkpoint.chk", "checkpoint.pdb")
    )
    payload["has_analysis"] = (run_dir / "analysis" / "summary.json").exists()
    return payload


def _states_for_status(status: dict[str, object]) -> tuple[str, ...]:
    states = [str(status.get("status", "unknown"))]
    if status.get("has_checkpoint"):
        states.append("checkpointed")
    if status.get("has_analysis"):
        states.append("analyzed")
    return tuple(states)
