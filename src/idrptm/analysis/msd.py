"""Mean-squared-displacement analysis."""

from __future__ import annotations

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

from idrptm.analysis._validation import as_masses, as_position_trajectory, validate_lag_count


def com_msd(
    positions: ArrayLike,
    masses: ArrayLike | None = None,
    max_lag: int | None = None,
) -> pd.DataFrame:
    """Return center-of-mass mean-squared displacement by lag."""

    trajectory = as_position_trajectory(positions)
    max_lag_value = validate_lag_count(trajectory.shape[0], max_lag)
    mass_values = as_masses(masses, trajectory.shape[1])
    centers = np.einsum("fij,i->fj", trajectory, mass_values) / mass_values.sum()

    rows: list[dict[str, float | int]] = []
    for lag in range(max_lag_value + 1):
        if lag == 0:
            squared = np.zeros(trajectory.shape[0], dtype=float)
        else:
            displacements = centers[lag:] - centers[:-lag]
            squared = np.sum(displacements * displacements, axis=-1)
        rows.append(
            {
                "lag": lag,
                "msd": float(squared.mean()),
                "n_origins": int(squared.size),
            }
        )
    return pd.DataFrame(rows, columns=["lag", "msd", "n_origins"])


def compute_msd(
    positions: ArrayLike,
    masses: ArrayLike | None = None,
    max_lag: int | None = None,
) -> pd.DataFrame:
    """Backward-compatible wrapper for :func:`com_msd`."""

    return com_msd(positions, masses=masses, max_lag=max_lag)
