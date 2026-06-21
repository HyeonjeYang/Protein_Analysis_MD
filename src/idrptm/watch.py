"""Human-friendly project progress summaries."""

from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass
from pathlib import Path

import yaml


@dataclass(frozen=True)
class ProjectWatchSnapshot:
    """One watch snapshot for a project."""

    project_dir: Path
    pipeline_status: dict[str, object]
    n_runs: int
    status_counts: dict[str, int]
    ready_trajectories: int
    total_dcd_mb: float
    frame_progress: tuple[int, int] | None


def summarize_watch(project_dir: str | Path) -> ProjectWatchSnapshot:
    """Collect status, DCD size, and frame progress for a project."""

    root = Path(project_dir)
    run_dirs = _run_dirs(root)
    counts: dict[str, int] = {}
    ready = 0
    dcd_bytes = 0
    frames_done = 0
    frames_total = 0
    for run_dir in run_dirs:
        status = _run_status(run_dir)
        counts[status] = counts.get(status, 0) + 1
        dcd = _find_dcd(run_dir)
        if dcd is not None:
            ready += 1
            dcd_bytes += dcd.stat().st_size
        progress = _frame_progress(run_dir)
        if progress is not None:
            done, total = progress
            frames_done += done
            frames_total += total
    frame_progress = (frames_done, frames_total) if frames_total else None
    return ProjectWatchSnapshot(
        project_dir=root,
        pipeline_status=_read_json(root / "pipeline_status.json"),
        n_runs=len(run_dirs),
        status_counts=counts,
        ready_trajectories=ready,
        total_dcd_mb=round(dcd_bytes / 1_000_000, 3),
        frame_progress=frame_progress,
    )


def format_watch(snapshot: ProjectWatchSnapshot) -> str:
    """Format one watch snapshot for terminal display."""

    pipeline = snapshot.pipeline_status
    lines = [
        f"project_dir: {snapshot.project_dir}",
        f"pipeline: {pipeline.get('stage', 'unknown')} / {pipeline.get('status', 'unknown')}",
        f"runs: {snapshot.n_runs}",
        "status_counts: "
        + ", ".join(
            f"{key}={value}" for key, value in sorted(snapshot.status_counts.items())
        ),
        f"trajectories_with_dcd: {snapshot.ready_trajectories}/{snapshot.n_runs}",
        f"total_dcd_mb: {snapshot.total_dcd_mb:.3f}",
    ]
    if snapshot.frame_progress is not None:
        done, total = snapshot.frame_progress
        percent = (100.0 * done / total) if total else 0.0
        lines.append(f"frames_estimate: {done}/{total} ({percent:.1f}%)")
    running = pipeline.get("running")
    pending = pipeline.get("pending")
    completed = pipeline.get("completed")
    if running:
        lines.append(f"running: {', '.join(str(item) for item in running)}")
    if pending is not None:
        lines.append(f"pending: {pending}")
    if isinstance(completed, list):
        lines.append(f"pipeline_completed_items: {len(completed)}")
    lines.append(f"updated_at: {time.strftime('%Y-%m-%d %H:%M:%S')}")
    return "\n".join(lines)


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


def _run_status(run_dir: Path) -> str:
    status = _read_json(run_dir / "run_status.json")
    if status.get("status"):
        return str(status["status"])
    if (run_dir / "run.py").exists():
        return "prepared"
    return "unknown"


def _frame_progress(run_dir: Path) -> tuple[int, int] | None:
    config = _read_yaml(run_dir / "config.yaml")
    steps = _as_int(config.get("steps"))
    wfreq = _as_int(config.get("wfreq"))
    if not steps or not wfreq:
        return None
    total = max(1, steps // wfreq)
    status = _read_json(run_dir / "run_status.json")
    progress = status.get("progress")
    if isinstance(progress, dict) and isinstance(progress.get("frames_written_estimate"), int):
        return min(int(progress["frames_written_estimate"]), total), total
    dcd = _find_dcd(run_dir)
    top = run_dir / "top.pdb"
    n_atoms = _count_atoms(top)
    if dcd is None or n_atoms is None:
        return 0, total
    return min(_estimate_dcd_frames(dcd, n_atoms), total), total


def _find_dcd(run_dir: Path) -> Path | None:
    candidates = sorted(run_dir.glob("*.dcd"), key=lambda path: path.stat().st_mtime)
    return candidates[-1] if candidates else None


def _count_atoms(topology_path: Path) -> int | None:
    if not topology_path.exists():
        return None
    count = 0
    with topology_path.open("r", encoding="utf-8", errors="ignore") as handle:
        for line in handle:
            if line.startswith(("ATOM", "HETATM")):
                count += 1
    return count or None


def _estimate_dcd_frames(dcd_path: Path, n_atoms: int) -> int:
    size = dcd_path.stat().st_size
    if size <= 0 or n_atoms <= 0:
        return 0
    frame_sizes = (12 * n_atoms + 80, 12 * n_atoms + 24)
    estimates = []
    for frame_size in frame_sizes:
        if size > frame_size:
            estimates.append(max(0, (size - 4096) // frame_size))
            estimates.append(max(0, (size - 1024) // frame_size))
            estimates.append(max(0, (size - 512) // frame_size))
    return max(estimates) if estimates else 0


def _read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _read_yaml(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _as_int(value: object) -> int | None:
    try:
        return int(value)
    except (TypeError, ValueError):
        return None
