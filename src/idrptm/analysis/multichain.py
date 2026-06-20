"""Multi-chain and multi-protein analysis utilities."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike, NDArray

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances
from idrptm.analysis.contacts import contact_map_from_positions
from idrptm.analysis.lifetime import contact_lifetime
from idrptm.analysis.msd import com_msd
from idrptm.analysis.ree import ree_timeseries
from idrptm.analysis.rg import rg_timeseries


def per_chain_rg(positions: ArrayLike, chain_ids: ArrayLike) -> pd.DataFrame:
    """Return Rg time series for each chain."""

    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    rows: list[pd.DataFrame] = []
    for chain in _ordered_unique(chains):
        table = pd.DataFrame(
            {
                "frame": np.arange(trajectory.shape[0], dtype=int),
                "chain_id": chain,
                "rg": rg_timeseries(trajectory[:, chains == chain, :]),
            }
        )
        rows.append(table)
    return pd.concat(rows, ignore_index=True)


def per_chain_ree(positions: ArrayLike, chain_ids: ArrayLike) -> pd.DataFrame:
    """Return end-to-end distance time series for each chain."""

    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    rows: list[pd.DataFrame] = []
    for chain in _ordered_unique(chains):
        subset = trajectory[:, chains == chain, :]
        values = (
            ree_timeseries(subset)
            if subset.shape[1] > 1
            else np.full(trajectory.shape[0], np.nan)
        )
        rows.append(
            pd.DataFrame(
                {
                    "frame": np.arange(trajectory.shape[0], dtype=int),
                    "chain_id": chain,
                    "ree": values,
                }
            )
        )
    return pd.concat(rows, ignore_index=True)


def intra_chain_contact_map(
    positions: ArrayLike,
    chain_ids: ArrayLike,
    cutoff: float,
    min_sequence_separation: int = 1,
) -> dict[str, NDArray[np.float64]]:
    """Return one residue contact map per chain."""

    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    maps: dict[str, NDArray[np.float64]] = {}
    for chain in _ordered_unique(chains):
        maps[chain] = contact_map_from_positions(
            trajectory[:, chains == chain, :],
            cutoff=cutoff,
            min_sequence_separation=min_sequence_separation,
        )
    return maps


def inter_chain_contact_map(
    positions: ArrayLike,
    chain_ids: ArrayLike,
    cutoff: float,
) -> NDArray[np.float64]:
    """Return a residue-residue contact probability map for inter-chain pairs only."""

    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    contacts = np.zeros((trajectory.shape[1], trajectory.shape[1]), dtype=float)
    inter_mask = chains[:, np.newaxis] != chains[np.newaxis, :]
    for frame in trajectory:
        contacts += (pairwise_distances(frame) <= cutoff) & inter_mask
    return contacts / trajectory.shape[0]


def com_distance_timeseries(
    positions: ArrayLike,
    chain_ids: ArrayLike,
) -> pd.DataFrame:
    """Return center-of-mass distances for each chain pair and frame."""

    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    centers = _chain_centers(trajectory, chains)
    rows: list[dict[str, float | int | str]] = []
    chain_names = list(centers)
    for frame in range(trajectory.shape[0]):
        for i, chain_i in enumerate(chain_names):
            for chain_j in chain_names[i + 1 :]:
                distance = np.linalg.norm(centers[chain_i][frame] - centers[chain_j][frame])
                rows.append(
                    {
                        "frame": frame,
                        "chain_i": chain_i,
                        "chain_j": chain_j,
                        "distance": float(distance),
                    }
                )
    return pd.DataFrame(rows, columns=["frame", "chain_i", "chain_j", "distance"])


def chain_com_msd(
    positions: ArrayLike,
    chain_ids: ArrayLike,
    max_lag: int | None = None,
) -> pd.DataFrame:
    """Return center-of-mass MSD for each chain."""

    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    rows: list[pd.DataFrame] = []
    for chain in _ordered_unique(chains):
        table = com_msd(trajectory[:, chains == chain, :], max_lag=max_lag)
        table["chain_id"] = chain
        rows.append(table[["chain_id", "lag", "msd", "n_origins"]])
    return pd.concat(rows, ignore_index=True)


def cluster_size_timeseries(
    positions: ArrayLike,
    chain_ids: ArrayLike,
    cutoff: float,
) -> pd.DataFrame:
    """Return chain-level cluster sizes using inter-chain contacts as graph edges."""

    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    chain_names = _ordered_unique(chains)
    rows: list[dict[str, float | int]] = []
    for frame_index, frame in enumerate(trajectory):
        parent = {chain: chain for chain in chain_names}
        for index_i, chain_i in enumerate(chain_names):
            for chain_j in chain_names[index_i + 1 :]:
                if _chains_contact(frame, chains, chain_i, chain_j, cutoff):
                    _union(parent, chain_i, chain_j)
        sizes: dict[str, int] = {}
        for chain in chain_names:
            root = _find(parent, chain)
            sizes[root] = sizes.get(root, 0) + 1
        values = np.asarray(list(sizes.values()), dtype=float)
        rows.append(
            {
                "frame": frame_index,
                "n_clusters": int(values.size),
                "largest_cluster_size": int(values.max()),
                "mean_cluster_size": float(values.mean()),
            }
        )
    return pd.DataFrame(
        rows,
        columns=["frame", "n_clusters", "largest_cluster_size", "mean_cluster_size"],
    )


def inter_protein_contact_lifetime(
    positions: ArrayLike,
    chain_ids: ArrayLike,
    cutoff: float,
    max_lag: int | None = None,
) -> pd.DataFrame:
    """Return contact lifetime for chain-pair inter-protein contacts."""

    if cutoff <= 0:
        raise ValueError("cutoff must be positive.")
    trajectory = as_position_trajectory(positions)
    chains = _chain_ids(chain_ids, trajectory.shape[1])
    chain_names = _ordered_unique(chains)
    pairs = [
        (chain_i, chain_j)
        for i, chain_i in enumerate(chain_names)
        for chain_j in chain_names[i + 1 :]
    ]
    states = np.zeros((trajectory.shape[0], len(pairs)), dtype=bool)
    for frame_index, frame in enumerate(trajectory):
        for pair_index, (chain_i, chain_j) in enumerate(pairs):
            states[frame_index, pair_index] = _chains_contact(
                frame,
                chains,
                chain_i,
                chain_j,
                cutoff,
            )
    return contact_lifetime(states, max_lag=max_lag)


def _chain_ids(chain_ids: ArrayLike, n_residues: int) -> NDArray[np.str_]:
    chains = np.asarray(chain_ids, dtype=str)
    if chains.shape != (n_residues,):
        raise ValueError(f"chain_ids must have shape ({n_residues},).")
    if len(set(chains.tolist())) < 1:
        raise ValueError("chain_ids must contain at least one chain.")
    return chains


def _ordered_unique(values: NDArray[np.str_]) -> list[str]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values.tolist():
        if value not in seen:
            seen.add(value)
            ordered.append(value)
    return ordered


def _chain_centers(
    trajectory: NDArray[np.float64],
    chains: NDArray[np.str_],
) -> dict[str, NDArray[np.float64]]:
    return {
        chain: trajectory[:, chains == chain, :].mean(axis=1)
        for chain in _ordered_unique(chains)
    }


def _chains_contact(
    frame: NDArray[np.float64],
    chains: NDArray[np.str_],
    chain_i: str,
    chain_j: str,
    cutoff: float,
) -> bool:
    first = frame[chains == chain_i]
    second = frame[chains == chain_j]
    distances = np.linalg.norm(first[:, np.newaxis, :] - second[np.newaxis, :, :], axis=-1)
    return bool(np.any(distances <= cutoff))


def _find(parent: dict[str, str], item: str) -> str:
    while parent[item] != item:
        parent[item] = parent[parent[item]]
        item = parent[item]
    return item


def _union(parent: dict[str, str], first: str, second: str) -> None:
    root_first = _find(parent, first)
    root_second = _find(parent, second)
    if root_first != root_second:
        parent[root_second] = root_first
