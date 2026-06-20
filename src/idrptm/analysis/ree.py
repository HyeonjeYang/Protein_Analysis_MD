"""End-to-end distance analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from idrptm.analysis._validation import as_position_trajectory


def ree_timeseries(
    positions: ArrayLike,
    start_index: int = 0,
    end_index: int = -1,
) -> NDArray[np.float64]:
    """Return end-to-end distance for each frame."""

    trajectory = as_position_trajectory(positions)
    n_residues = trajectory.shape[1]
    start = _normalize_index(start_index, n_residues)
    end = _normalize_index(end_index, n_residues)
    if start == end:
        raise ValueError("start_index and end_index must refer to different residues.")
    delta = trajectory[:, end, :] - trajectory[:, start, :]
    return np.linalg.norm(delta, axis=-1)


def compute_ree(
    positions: ArrayLike,
    start_index: int = 0,
    end_index: int = -1,
) -> NDArray[np.float64]:
    """Backward-compatible wrapper for :func:`ree_timeseries`."""

    return ree_timeseries(positions, start_index=start_index, end_index=end_index)


def _normalize_index(index: int, length: int) -> int:
    normalized = index + length if index < 0 else index
    if normalized < 0 or normalized >= length:
        raise ValueError(f"Residue index {index} is out of bounds for length {length}.")
    return normalized
