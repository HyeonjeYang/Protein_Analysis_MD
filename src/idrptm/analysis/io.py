"""Trajectory I/O and common trajectory data model."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
from numpy.typing import NDArray

from idrptm.units import CANONICAL_UNITS, angstrom_to_nm


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
    length_unit: str = CANONICAL_UNITS["length"]
    input_position_unit: str = CANONICAL_UNITS["length"]
    canonical_position_unit: str = CANONICAL_UNITS["length"]

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
    engine: Literal["mdtraj", "mdanalysis"] = "mdtraj",
) -> TrajectoryData:
    """Load CALVADOS ``top.pdb`` and ``trajectory.dcd`` with mdtraj or MDAnalysis."""

    run_path = Path(run_dir)
    topology_path = Path(topology) if topology is not None else run_path / "top.pdb"
    trajectory_path = Path(trajectory) if trajectory is not None else run_path / "trajectory.dcd"
    _require_file(topology_path, "topology PDB")
    _require_file(trajectory_path, "trajectory DCD")

    if engine == "mdanalysis":
        return _load_with_mdanalysis(topology_path, trajectory_path)
    if engine != "mdtraj":
        raise ValueError(f"Unsupported trajectory reader: {engine}")
    return _load_with_mdtraj(topology_path, trajectory_path)


def _load_with_mdtraj(topology_path: Path, trajectory_path: Path) -> TrajectoryData:
    """Load a trajectory with mdtraj, whose coordinates are already in nm."""

    try:
        import mdtraj as md
    except ImportError as exc:
        raise ImportError(
            "mdtraj is required to load CALVADOS trajectories. "
            "Install idr-ptm-md with its analysis dependencies."
        ) from exc

    loaded = md.load(str(trajectory_path), top=str(topology_path))
    positions = np.asarray(loaded.xyz, dtype=float)
    _validate_positions(positions)
    time = np.asarray(loaded.time, dtype=float) if loaded.time is not None else None
    return TrajectoryData(
        positions=positions,
        topology_path=topology_path,
        trajectory_path=trajectory_path,
        time_ps=time,
        length_unit=CANONICAL_UNITS["length"],
        input_position_unit=CANONICAL_UNITS["length"],
        canonical_position_unit=CANONICAL_UNITS["length"],
    )


def _load_with_mdanalysis(topology_path: Path, trajectory_path: Path) -> TrajectoryData:
    """Load with MDAnalysis and convert Angstrom coordinates to canonical nm."""

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise ImportError(
            "MDAnalysis is required when engine='mdanalysis'. Install MDAnalysis "
            "or use the default mdtraj trajectory reader."
        ) from exc

    universe = mda.Universe(str(topology_path), str(trajectory_path))
    frames: list[NDArray[np.float64]] = []
    times: list[float] = []
    for timestep in universe.trajectory:
        frames.append(np.asarray(universe.atoms.positions, dtype=float).copy())
        if hasattr(timestep, "time"):
            times.append(float(timestep.time))
    positions_angstrom = np.asarray(frames, dtype=float)
    _validate_positions(positions_angstrom)
    positions_nm = np.asarray(angstrom_to_nm(positions_angstrom), dtype=float)
    time = np.asarray(times, dtype=float) if len(times) == len(frames) else None
    return TrajectoryData(
        positions=positions_nm,
        topology_path=topology_path,
        trajectory_path=trajectory_path,
        time_ps=time,
        length_unit=CANONICAL_UNITS["length"],
        input_position_unit="angstrom",
        canonical_position_unit=CANONICAL_UNITS["length"],
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


def _validate_positions(positions: NDArray[np.float64]) -> None:
    if positions.ndim != 3 or positions.shape[-1] != 3:
        raise ValueError(
            f"Loaded trajectory has unexpected coordinate shape {positions.shape}; "
            "expected (n_frames, n_atoms, 3)."
        )


def _require_file(path: Path, label: str) -> None:
    if not path.is_file():
        raise TrajectoryFileError(
            f"Missing {label}: {path}. Expected CALVADOS outputs named "
            "'top.pdb' and 'trajectory.dcd' in the run directory, or pass explicit paths."
        )
