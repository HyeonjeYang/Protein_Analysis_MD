"""Contact-map analysis."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances


def contact_map_from_positions(
    positions: ArrayLike,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> NDArray[np.float64]:
    """Return contact probability matrix from positions.

    A pair is in contact when its distance is less than or equal to ``cutoff``.
    The returned matrix contains contact frequencies over frames and is symmetric.
    Pairs with sequence separation smaller than ``min_sequence_separation`` are
    set to zero. The default excludes only self-contacts.
    """

    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    if min_sequence_separation < 1:
        raise ValueError("min_sequence_separation must be at least 1.")

    trajectory = as_position_trajectory(positions)
    n_residues = trajectory.shape[1]
    contacts = np.zeros((n_residues, n_residues), dtype=float)
    valid_pairs = _valid_pair_mask(n_residues, min_sequence_separation)

    for frame in trajectory:
        frame_contacts = pairwise_distances(frame) <= cutoff
        contacts += frame_contacts & valid_pairs

    return contacts / trajectory.shape[0]


def compute_contact_map(
    positions: ArrayLike,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> NDArray[np.float64]:
    """Backward-compatible wrapper for :func:`contact_map_from_positions`."""

    return contact_map_from_positions(
        positions,
        cutoff=cutoff,
        min_sequence_separation=min_sequence_separation,
    )


def _valid_pair_mask(n_residues: int, min_sequence_separation: int) -> NDArray[np.bool_]:
    indices = np.arange(n_residues)
    return np.abs(indices[:, np.newaxis] - indices[np.newaxis, :]) >= min_sequence_separation
