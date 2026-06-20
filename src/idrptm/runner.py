"""Local execution placeholders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

RunStatus = Literal["planned", "submitted", "completed", "failed"]


@dataclass(frozen=True)
class RunPlan:
    """A planned execution command for a prepared run directory."""

    run_dir: Path
    command: tuple[str, ...]
    status: RunStatus = "planned"


def plan_local_run(run_dir: str | Path, python_executable: str = "python") -> RunPlan:
    """Create a local execution plan for a future CALVADOS run."""

    path = Path(run_dir)
    return RunPlan(run_dir=path, command=(python_executable, str(path / "run.py")))
