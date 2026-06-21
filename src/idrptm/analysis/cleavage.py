"""Analysis helpers for pre-cleaved sequence and fragment-mixture simulations."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances
from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.multichain import (
    cluster_size_timeseries,
    inter_chain_contact_map,
    per_chain_ree,
    per_chain_rg,
)


def fragment_resolved_rg(positions: ArrayLike, fragment_ids: ArrayLike) -> pd.DataFrame:
    """Return Rg time series per fragment."""

    return per_chain_rg(positions, fragment_ids).rename(columns={"chain_id": "fragment_id"})


def fragment_resolved_ree(positions: ArrayLike, fragment_ids: ArrayLike) -> pd.DataFrame:
    """Return end-to-end distance time series per fragment."""

    return per_chain_ree(positions, fragment_ids).rename(columns={"chain_id": "fragment_id"})


def inter_fragment_contact_map(
    positions: ArrayLike,
    fragment_ids: ArrayLike,
    cutoff: float,
) -> NDArray[np.float64]:
    """Return residue contact probabilities between fragments only."""

    return inter_chain_contact_map(positions, fragment_ids, cutoff=cutoff)


def fragment_cluster_size(
    positions: ArrayLike,
    fragment_ids: ArrayLike,
    cutoff: float,
) -> pd.DataFrame:
    """Return fragment cluster-size time series."""

    return cluster_size_timeseries(positions, fragment_ids, cutoff=cutoff)


def intact_vs_cleaved_delta_contact_map(
    intact_positions: ArrayLike,
    cleaved_positions: ArrayLike,
    *,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> NDArray[np.float64]:
    """Return cleaved minus intact contact probability maps."""

    intact = contact_map_from_positions(
        intact_positions,
        cutoff=cutoff,
        min_sequence_separation=min_sequence_separation,
    )
    cleaved = contact_map_from_positions(
        cleaved_positions,
        cutoff=cutoff,
        min_sequence_separation=min_sequence_separation,
    )
    if intact.shape != cleaved.shape:
        raise ValueError("Intact and cleaved contact maps must have the same shape.")
    return cleaved - intact


def original_sequence_coordinate_contact_map(
    positions: ArrayLike,
    original_sequence_indices: ArrayLike,
    *,
    n_original_residues: int,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> NDArray[np.float64]:
    """Project fragment coordinates onto original-sequence contact-map coordinates."""

    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    if n_original_residues < 1:
        raise ValueError("n_original_residues must be positive.")
    trajectory = as_position_trajectory(positions)
    original = np.asarray(original_sequence_indices, dtype=int)
    if original.shape != (trajectory.shape[1],):
        raise ValueError(f"original_sequence_indices must have shape ({trajectory.shape[1]},).")
    if np.any(original < 1) or np.any(original > n_original_residues):
        raise ValueError("original_sequence_indices must be one-based positions in range.")

    contact_map = np.zeros((n_original_residues, n_original_residues), dtype=float)
    for frame in trajectory:
        contacts = pairwise_distances(frame) <= cutoff
        for i in range(trajectory.shape[1]):
            orig_i = original[i] - 1
            for j in range(i + min_sequence_separation, trajectory.shape[1]):
                orig_j = original[j] - 1
                if abs(orig_i - orig_j) < min_sequence_separation:
                    continue
                if contacts[i, j]:
                    contact_map[orig_i, orig_j] += 1.0
                    contact_map[orig_j, orig_i] += 1.0
    return contact_map / trajectory.shape[0]
