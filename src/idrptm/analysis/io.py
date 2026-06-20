"""Trajectory I/O placeholders."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class TrajectoryInput:
    """Files required to load one trajectory."""

    trajectory: Path
    topology: Path | None = None


def load_trajectory(_: TrajectoryInput) -> object:
    """Load a trajectory in Stage 2."""

    raise NotImplementedError("Trajectory loading will be implemented in Stage 2.")
