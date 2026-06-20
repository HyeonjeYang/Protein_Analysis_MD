"""Trajectory I/O and common trajectory data model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class TrajectoryInput:
    """Files required to load one trajectory."""

    trajectory: Path
    topology: Path | None = None


@dataclass(frozen=True)
class TrajectoryData:
    """Common coordinate container used by the pure analysis core."""

    positions: NDArray[np.float64]
    topology_path: Path
    trajectory_path: Path
    time_ps: NDArray[np.float64] | None = None
    length_unit: str = "nm"

    @property
    def n_frames(self) -> int:
        return int(self.positions.shape[0])

    @property
    def n_residues(self) -> int:
        return int(self.positions.shape[1])


class TrajectoryFileError(FileNotFoundError):
    """Raised when a required trajectory input file is missing."""


def load_calvados_trajectory(
    run_dir: str | Path,
    *,
    topology: str | Path | None = None,
    trajectory: str | Path | None = None,
) -> TrajectoryData:
    """Load CALVADOS ``top.pdb`` and ``trajectory.dcd`` with mdtraj."""

    run_path = Path(run_dir)
    topology_path = Path(topology) if topology is not None else run_path / "top.pdb"
    trajectory_path = Path(trajectory) if trajectory is not None else run_path / "trajectory.dcd"
    _require_file(topology_path, "topology PDB")
    _require_file(trajectory_path, "trajectory DCD")

    try:
        import mdtraj as md
    except ImportError as exc:
        raise ImportError(
            "mdtraj is required to load CALVADOS trajectories. "
            "Install idr-ptm-md with its analysis dependencies."
        ) from exc

    loaded = md.load(str(trajectory_path), top=str(topology_path))
    positions = np.asarray(loaded.xyz, dtype=float)
    if positions.ndim != 3 or positions.shape[-1] != 3:
        raise ValueError(
            f"Loaded trajectory has unexpected coordinate shape {positions.shape}; "
            "expected (n_frames, n_atoms, 3)."
        )
    time = np.asarray(loaded.time, dtype=float) if loaded.time is not None else None
    return TrajectoryData(
        positions=positions,
        topology_path=topology_path,
        trajectory_path=trajectory_path,
        time_ps=time,
    )


def load_trajectory(trajectory_input: TrajectoryInput) -> TrajectoryData:
    """Load a trajectory from explicit files with mdtraj."""

    if trajectory_input.topology is None:
        raise ValueError("A topology PDB is required to load a trajectory.")
    return load_calvados_trajectory(
        trajectory_input.trajectory.parent,
        topology=trajectory_input.topology,
        trajectory=trajectory_input.trajectory,
    )


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise TrajectoryFileError(
            f"Missing {label}: {path}. Expected CALVADOS outputs named "
            "'top.pdb' and 'trajectory.dcd' in the run directory, or pass explicit paths."
        )
