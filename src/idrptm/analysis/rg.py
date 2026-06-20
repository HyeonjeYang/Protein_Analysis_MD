"""Radius-of-gyration analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from idrptm.analysis._validation import as_masses, as_position_trajectory


def rg_timeseries(positions: ArrayLike, masses: ArrayLike | None = None) -> NDArray[np.float64]:
    """Return radius of gyration for each frame.

    Parameters
    ----------
    positions
        Coordinates with shape ``(n_frames, n_residues, 3)`` or one frame with
        shape ``(n_residues, 3)``.
    masses
        Optional positive residue/bead masses. Uniform masses are used by default.
    """

    trajectory = as_position_trajectory(positions)
    mass_values = as_masses(masses, trajectory.shape[1])
    total_mass = mass_values.sum()
    centers = np.einsum("fij,i->fj", trajectory, mass_values) / total_mass
    centered = trajectory - centers[:, np.newaxis, :]
    squared_distances = np.sum(centered * centered, axis=-1)
    return np.sqrt(np.einsum("fi,i->f", squared_distances, mass_values) / total_mass)


def compute_rg(positions: ArrayLike, masses: ArrayLike | None = None) -> NDArray[np.float64]:
    """Backward-compatible wrapper for :func:`rg_timeseries`."""

    return rg_timeseries(positions, masses=masses)
