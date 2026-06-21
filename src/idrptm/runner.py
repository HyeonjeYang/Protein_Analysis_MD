"""Local execution wrappers for prepared CALVADOS run directories."""

from __future__ import annotations

import csv
import json
import subprocess
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import yaml

from idrptm.provenance import execution_environment, write_parameter_txt

RunStatus = Literal["planned", "submitted", "completed", "failed"]
RunPhase = Literal["equilibration", "production", "all"]


@dataclass(frozen=True)
class RunPlan:
    """A planned execution command for a prepared run directory."""

    run_dir: Path
    command: tuple[str, ...]
    phase: RunPhase = "all"
    status: RunStatus = "planned"


@dataclass(frozen=True)
class RunResult:
    """Execution result for one prepared run directory."""

    run_dir: Path
    command: tuple[str, ...]
    phase: RunPhase
    status: RunStatus
    returncode: int | None
    status_json: Path


def plan_local_run(
    run_dir: str | Path,
    *,
    phase: RunPhase = "all",
    python_executable: str = "python",
) -> RunPlan:
    """Create a local execution plan for one prepared CALVADOS run directory."""

    path = Path(run_dir)
    return RunPlan(
        run_dir=path,
        command=_run_command(path, phase=phase, python_executable=python_executable),
        phase=phase,
    )


def plan_local_runs(
    target: str | Path,
    *,
    phase: RunPhase = "all",
    all_runs: bool = False,
    python_executable: str = "python",
) -> tuple[RunPlan, ...]:
    """Create local execution plans for one run directory or a project directory."""

    return tuple(
        plan_local_run(
            run_dir,
            phase=phase,
            python_executable=python_executable,
        )
        for run_dir in discover_run_directories(target, all_runs=all_runs)
    )


def execute_local_runs(
    plans: tuple[RunPlan, ...],
    *,
    progress: bool = True,
    progress_interval_s: float = 5.0,
) -> tuple[RunResult, ...]:
    """Execute planned local runs and write ``run_status.json`` files."""

    results: list[RunResult] = []
    for plan in plans:
        results.append(
            _execute_plan(
                plan,
                progress=progress,
                progress_interval_s=progress_interval_s,
            )
        )
    return tuple(results)


def discover_run_directories(target: str | Path, *, all_runs: bool = False) -> tuple[Path, ...]:
    """Resolve a target path into prepared run directories."""

    path = Path(target)
    if all_runs:
        manifest_path = path / "manifest.csv"
        if manifest_path.exists():
            return _run_dirs_from_manifest(path, manifest_path)
        runs_dir = path / "runs"
        if runs_dir.exists():
            return tuple(sorted(_iter_prepared_run_dirs(runs_dir)))
        raise FileNotFoundError(
            f"Could not find manifest.csv or runs/ under project directory {path}."
        )
    if not path.exists():
        raise FileNotFoundError(f"Run directory does not exist: {path}")
    return (path,)


def write_planned_status(plan: RunPlan) -> Path:
    """Record a dry-run status file for a planned run."""

    status_path = plan.run_dir / "run_status.json"
    _write_status(
        status_path,
        {
            "status": "planned",
            "phase": plan.phase,
            "run_dir": str(plan.run_dir),
            "command": list(plan.command),
            "timestamp": _now_iso(),
        },
    )
    return status_path


def _run_command(
    run_dir: Path,
    *,
    phase: RunPhase,
    python_executable: str,
) -> tuple[str, ...]:
    script = _script_for_phase(run_dir, phase)
    executable_path = Path(python_executable)
    executable = (
        str(executable_path.absolute())
        if executable_path.parent != Path(".")
        else python_executable
    )
    if script.name == "run.py":
        return (
            executable,
            script.name,
            "--path",
            ".",
            "--config",
            "config.yaml",
            "--components",
            "components.yaml",
        )
    return (executable, script.name)


def _script_for_phase(run_dir: Path, phase: RunPhase) -> Path:
    phase_scripts = {
        "equilibration": "run_equilibration.py",
        "production": "run_production.py",
        "all": "run_all.py",
    }
    preferred = run_dir / phase_scripts[phase]
    if preferred.exists():
        return preferred
    legacy_script = run_dir / "run.py"
    if phase == "all" and legacy_script.exists():
        return legacy_script
    raise FileNotFoundError(
        f"No run script for phase {phase!r} in {run_dir}. "
        "Run 'pamd prepare' first, or use --phase all for legacy run.py directories."
    )


def _execute_plan(
    plan: RunPlan,
    *,
    progress: bool,
    progress_interval_s: float,
) -> RunResult:
    status_path = plan.run_dir / "run_status.json"
    started_at = _now_iso()
    started_perf = time.perf_counter()
    stdout_path = plan.run_dir / "execution_stdout.log"
    stderr_path = plan.run_dir / "execution_stderr.log"
    environment = execution_environment()
    _write_status(
        status_path,
        {
            "status": "submitted",
            "phase": plan.phase,
            "run_dir": str(plan.run_dir),
            "command": list(plan.command),
            "started_at": started_at,
            "execution_environment": environment,
        },
    )
    _write_runtime_parameters(plan, status_path)
    try:
        with stdout_path.open("w", encoding="utf-8") as stdout, stderr_path.open(
            "w",
            encoding="utf-8",
        ) as stderr:
            process = subprocess.Popen(
                plan.command,
                cwd=plan.run_dir,
                stdout=stdout,
                stderr=stderr,
                text=True,
            )
            bar = _make_progress_bar(plan, enabled=progress)
            returncode: int | None = None
            try:
                while True:
                    snapshot = _progress_snapshot(plan.run_dir)
                    _update_progress_bar(bar, snapshot)
                    _write_status(
                        status_path,
                        {
                            "status": "submitted",
                            "phase": plan.phase,
                            "run_dir": str(plan.run_dir),
                            "command": list(plan.command),
                            "started_at": started_at,
                            "updated_at": _now_iso(),
                            "elapsed_wall_s": round(time.perf_counter() - started_perf, 3),
                            "progress": snapshot,
                            "execution_environment": environment,
                        },
                    )
                    try:
                        returncode = process.wait(timeout=progress_interval_s)
                        break
                    except subprocess.TimeoutExpired:
                        continue
            finally:
                if bar is not None:
                    bar.close()
            if returncode is None:
                returncode = process.wait()
    except Exception as exc:
        elapsed_wall_s = round(time.perf_counter() - started_perf, 3)
        _write_status(
            status_path,
            {
                "status": "failed",
                "phase": plan.phase,
                "run_dir": str(plan.run_dir),
                "command": list(plan.command),
                "started_at": started_at,
                "ended_at": _now_iso(),
                "elapsed_wall_s": elapsed_wall_s,
                "elapsed_wall_h": round(elapsed_wall_s / 3600.0, 6),
                "error": str(exc),
                "execution_environment": environment,
            },
        )
        _write_runtime_parameters(plan, status_path)
        return RunResult(
            run_dir=plan.run_dir,
            command=plan.command,
            phase=plan.phase,
            status="failed",
            returncode=None,
            status_json=status_path,
        )

    ended_at = _now_iso()
    elapsed_wall_s = round(time.perf_counter() - started_perf, 3)
    status: RunStatus = "completed" if returncode == 0 else "failed"
    _write_status(
        status_path,
        {
            "status": status,
            "phase": plan.phase,
            "run_dir": str(plan.run_dir),
            "command": list(plan.command),
            "started_at": started_at,
            "ended_at": ended_at,
            "elapsed_wall_s": elapsed_wall_s,
            "elapsed_wall_h": round(elapsed_wall_s / 3600.0, 6),
            "returncode": returncode,
            "progress": _progress_snapshot(plan.run_dir),
            "stdout_log": str(stdout_path),
            "stderr_log": str(stderr_path),
            "stdout_tail": _tail(stdout_path),
            "stderr_tail": _tail(stderr_path),
            "execution_environment": environment,
        },
    )
    _write_runtime_parameters(plan, status_path)
    return RunResult(
        run_dir=plan.run_dir,
        command=plan.command,
        phase=plan.phase,
        status=status,
        returncode=returncode,
        status_json=status_path,
    )


def _run_dirs_from_manifest(root: Path, manifest_path: Path) -> tuple[Path, ...]:
    run_dirs: list[Path] = []
    with manifest_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if row.get("run_dir"):
                run_dirs.append(root / row["run_dir"])
            elif row.get("metadata_path"):
                run_dirs.append((root / row["metadata_path"]).parent)
    return tuple(sorted(dict.fromkeys(run_dirs)))


def _iter_prepared_run_dirs(root: Path) -> tuple[Path, ...]:
    scripts = {
        "run.py",
        "run_all.py",
        "run_equilibration.py",
        "run_production.py",
    }
    return tuple(
        path
        for path in root.iterdir()
        if path.is_dir() and any((path / script).exists() for script in scripts)
    )


def _write_runtime_parameters(plan: RunPlan, status_path: Path) -> None:
    payload = {
        "run_dir": str(plan.run_dir),
        "phase": plan.phase,
        "command": list(plan.command),
        "config": _read_yaml(plan.run_dir / "config.yaml"),
        "components": _read_yaml(plan.run_dir / "components.yaml"),
        "metadata": _read_json(plan.run_dir / "metadata.json"),
        "run_status": _read_json(status_path),
    }
    write_parameter_txt(
        plan.run_dir / "parameters.txt",
        payload,
        title="protein_analysis_md runtime parameters",
    )


def _progress_snapshot(run_dir: Path) -> dict[str, object]:
    config = _read_yaml(run_dir / "config.yaml")
    steps = _as_int(config.get("steps"))
    wfreq = _as_int(config.get("wfreq"))
    total_frames = steps // wfreq if steps and wfreq else None
    dcd_path = _find_dcd_path(run_dir, str(config.get("sysname") or ""))
    n_atoms = _count_atoms(run_dir / "top.pdb")
    frames_written = (
        _estimate_dcd_frames(dcd_path, n_atoms)
        if dcd_path is not None and n_atoms is not None
        else None
    )
    if frames_written is not None and total_frames is not None:
        frames_written = min(frames_written, total_frames)
    step_estimate = frames_written * wfreq if frames_written is not None and wfreq else None
    return {
        "dcd_path": str(dcd_path) if dcd_path is not None else None,
        "n_atoms": n_atoms,
        "frames_written_estimate": frames_written,
        "total_frames": total_frames,
        "step_estimate": step_estimate,
        "total_steps": steps,
        "wfreq": wfreq,
    }


def _make_progress_bar(plan: RunPlan, *, enabled: bool):
    if not enabled:
        return None
    snapshot = _progress_snapshot(plan.run_dir)
    total = snapshot.get("total_frames")
    try:
        from tqdm import tqdm
    except Exception:
        return None
    return tqdm(
        total=int(total) if isinstance(total, int) and total > 0 else None,
        unit="frame",
        desc=plan.run_dir.name,
        dynamic_ncols=True,
        leave=True,
    )


def _update_progress_bar(bar: object | None, snapshot: dict[str, object]) -> None:
    if bar is None:
        return
    frames = snapshot.get("frames_written_estimate")
    if isinstance(frames, int):
        bar.n = frames
    postfix = {}
    step = snapshot.get("step_estimate")
    total_steps = snapshot.get("total_steps")
    if isinstance(step, int) and isinstance(total_steps, int):
        postfix["step"] = f"{step}/{total_steps}"
    if postfix and hasattr(bar, "set_postfix"):
        bar.set_postfix(postfix)
    if hasattr(bar, "refresh"):
        bar.refresh()


def _find_dcd_path(run_dir: Path, sysname: str) -> Path | None:
    preferred = []
    if sysname:
        preferred.append(run_dir / f"{sysname}.dcd")
    preferred.append(run_dir / "trajectory.dcd")
    for path in preferred:
        if path.exists():
            return path
    candidates = sorted(run_dir.glob("*.dcd"), key=lambda path: path.stat().st_mtime)
    return candidates[-1] if candidates else None


def _count_atoms(topology_path: Path) -> int | None:
    if not topology_path.exists():
        return None
    count = 0
    try:
        with topology_path.open("r", encoding="utf-8", errors="ignore") as handle:
            for line in handle:
                if line.startswith(("ATOM", "HETATM")):
                    count += 1
    except OSError:
        return None
    return count or None


def _estimate_dcd_frames(dcd_path: Path, n_atoms: int) -> int | None:
    try:
        size = dcd_path.stat().st_size
    except OSError:
        return None
    if size <= 0 or n_atoms <= 0:
        return None
    frame_sizes = (12 * n_atoms + 80, 12 * n_atoms + 24)
    estimates = []
    for frame_size in frame_sizes:
        if size > frame_size:
            estimates.append(max(0, (size - 4096) // frame_size))
            estimates.append(max(0, (size - 1024) // frame_size))
            estimates.append(max(0, (size - 512) // frame_size))
    return max(estimates) if estimates else 0


def _read_yaml(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = yaml.safe_load(path.read_text(encoding="utf-8")) or {}
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _read_json(path: Path) -> dict[str, object]:
    if not path.exists():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def _as_int(value: object) -> int | None:
    if value in (None, ""):
        return None
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


def _tail(path: Path, limit: int = 4000) -> str:
    if not path.exists():
        return ""
    try:
        text = path.read_text(encoding="utf-8", errors="replace")
    except OSError:
        return ""
    return text[-limit:]


def _write_status(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
