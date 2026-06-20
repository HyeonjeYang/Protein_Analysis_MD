"""Shared validation helpers for pure coordinate analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def as_position_trajectory(positions: ArrayLike) -> NDArray[np.float64]:
    """Return positions as a ``(n_frames, n_residues, 3)`` float array."""

    array = np.asarray(positions, dtype=float)
    if array.ndim == 2:
        array = array[np.newaxis, :, :]
    if array.ndim != 3 or array.shape[-1] != 3:
        raise ValueError("positions must have shape (n_frames, n_residues, 3) or (n_residues, 3).")
    if array.shape[0] < 1:
        raise ValueError("positions must contain at least one frame.")
    if array.shape[1] < 1:
        raise ValueError("positions must contain at least one residue/bead.")
    return array


def as_masses(masses: ArrayLike | None, n_residues: int) -> NDArray[np.float64]:
    """Return residue masses as a normalized one-dimensional float array."""

    if masses is None:
        values = np.ones(n_residues, dtype=float)
    else:
        values = np.asarray(masses, dtype=float)
    if values.shape != (n_residues,):
        raise ValueError(f"masses must have shape ({n_residues},).")
    if np.any(values <= 0):
        raise ValueError("masses must be positive.")
    return values


def pairwise_distances(frame: NDArray[np.float64]) -> NDArray[np.float64]:
    """Return an all-pairs Euclidean distance matrix for one frame."""

    delta = frame[:, np.newaxis, :] - frame[np.newaxis, :, :]
    return np.linalg.norm(delta, axis=-1)


def validate_lag_count(n_frames: int, max_lag: int | None) -> int:
    """Validate and normalize a maximum lag value."""

    if max_lag is None:
        return n_frames - 1
    if max_lag < 0:
        raise ValueError("max_lag must be non-negative.")
    return min(max_lag, n_frames - 1)
