"""Local execution wrappers for prepared CALVADOS run directories."""

from __future__ import annotations

import csv
import json
import subprocess
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

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


def execute_local_runs(plans: tuple[RunPlan, ...]) -> tuple[RunResult, ...]:
    """Execute planned local runs and write ``run_status.json`` files."""

    results: list[RunResult] = []
    for plan in plans:
        results.append(_execute_plan(plan))
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


def _execute_plan(plan: RunPlan) -> RunResult:
    status_path = plan.run_dir / "run_status.json"
    started_at = _now_iso()
    _write_status(
        status_path,
        {
            "status": "submitted",
            "phase": plan.phase,
            "run_dir": str(plan.run_dir),
            "command": list(plan.command),
            "started_at": started_at,
        },
    )
    try:
        completed = subprocess.run(
            plan.command,
            cwd=plan.run_dir,
            check=False,
            capture_output=True,
            text=True,
        )
    except Exception as exc:
        _write_status(
            status_path,
            {
                "status": "failed",
                "phase": plan.phase,
                "run_dir": str(plan.run_dir),
                "command": list(plan.command),
                "started_at": started_at,
                "ended_at": _now_iso(),
                "error": str(exc),
            },
        )
        return RunResult(
            run_dir=plan.run_dir,
            command=plan.command,
            phase=plan.phase,
            status="failed",
            returncode=None,
            status_json=status_path,
        )

    status: RunStatus = "completed" if completed.returncode == 0 else "failed"
    _write_status(
        status_path,
        {
            "status": status,
            "phase": plan.phase,
            "run_dir": str(plan.run_dir),
            "command": list(plan.command),
            "started_at": started_at,
            "ended_at": _now_iso(),
            "returncode": completed.returncode,
            "stdout_tail": completed.stdout[-4000:],
            "stderr_tail": completed.stderr[-4000:],
        },
    )
    return RunResult(
        run_dir=plan.run_dir,
        command=plan.command,
        phase=plan.phase,
        status=status,
        returncode=completed.returncode,
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


def _write_status(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _now_iso() -> str:
    return datetime.now(timezone.utc).isoformat()
