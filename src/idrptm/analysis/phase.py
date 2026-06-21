"""Synthetic phase-separation analysis helpers."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from idrptm.analysis._validation import as_position_trajectory


def density_profile_z(
    positions: ArrayLike,
    *,
    box_z_nm: float,
    bins: int = 50,
) -> pd.DataFrame:
    """Compute a simple bead density profile along z."""

    trajectory = as_position_trajectory(positions)
    z = trajectory[:, :, 2].ravel()
    counts, edges = np.histogram(z, bins=bins, range=(0.0, box_z_nm))
    centers = 0.5 * (edges[:-1] + edges[1:])
    return pd.DataFrame({"z_nm": centers, "count": counts})


def largest_cluster_fraction(cluster_sizes: ArrayLike) -> float:
    """Return largest cluster divided by total molecules."""

    values = np.asarray(cluster_sizes, dtype=float)
    if values.size == 0 or values.sum() == 0:
        return 0.0
    return float(values.max() / values.sum())
