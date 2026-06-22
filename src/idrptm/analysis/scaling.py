"""Polymer-scaling analysis."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy import stats

from idrptm.analysis._validation import as_position_trajectory, pairwise_distances


@dataclass(frozen=True)
class FloryFit:
    """Log-log fit result for ``R(s) = prefactor * s**nu``."""

    nu: float
    intercept: float
    prefactor: float
    r_value: float
    p_value: float
    stderr: float
    n_points: int


def internal_distance_scaling(positions: ArrayLike) -> pd.DataFrame:
    """Return internal-distance scaling metrics by sequence separation."""

    trajectory = as_position_trajectory(positions)
    n_residues = trajectory.shape[1]
    accumulators = {separation: [] for separation in range(1, n_residues)}

    for frame in trajectory:
        distances = pairwise_distances(frame)
        for separation in range(1, n_residues):
            accumulators[separation].extend(np.diag(distances, k=separation).tolist())

    rows = [
        {
            "s": separation,
            "distance": float(np.mean(values)),
            "mean_distance_nm": float(np.mean(values)),
            "rms_distance_nm": float(np.sqrt(np.mean(np.square(values)))),
            "mean_square_distance_nm2": float(np.mean(np.square(values))),
            "std_distance_nm": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
            "n_pairs": len(values),
        }
        for separation, values in accumulators.items()
    ]
    return pd.DataFrame(
        rows,
        columns=[
            "s",
            "distance",
            "mean_distance_nm",
            "rms_distance_nm",
            "mean_square_distance_nm2",
            "std_distance_nm",
            "n_pairs",
        ],
    )


def fit_flory_exponent(
    scaling: pd.DataFrame | None = None,
    *,
    s: ArrayLike | None = None,
    distances: ArrayLike | None = None,
    min_s: int | None = None,
    max_s: int | None = None,
) -> FloryFit:
    """Fit the Flory exponent ``nu`` from internal-distance scaling data."""

    if scaling is not None:
        distance_column = "distance" if "distance" in scaling else "mean_distance_nm"
        if "s" not in scaling or distance_column not in scaling:
            raise ValueError("scaling must contain 's' and a distance column.")
        s_values = scaling["s"].to_numpy(dtype=float)
        distance_values = scaling[distance_column].to_numpy(dtype=float)
    else:
        if s is None or distances is None:
            raise ValueError("Provide either scaling or both s and distances.")
        s_values = np.asarray(s, dtype=float)
        distance_values = np.asarray(distances, dtype=float)

    if s_values.shape != distance_values.shape:
        raise ValueError("s and distances must have the same shape.")

    mask = np.isfinite(s_values) & np.isfinite(distance_values)
    mask &= (s_values > 0) & (distance_values > 0)
    if min_s is not None:
        mask &= s_values >= min_s
    if max_s is not None:
        mask &= s_values <= max_s

    fit_s = s_values[mask]
    fit_distances = distance_values[mask]
    if fit_s.size < 2:
        raise ValueError("At least two positive finite points are required for fitting.")

    result = stats.linregress(np.log(fit_s), np.log(fit_distances))
    return FloryFit(
        nu=float(result.slope),
        intercept=float(result.intercept),
        prefactor=float(np.exp(result.intercept)),
        r_value=float(result.rvalue),
        p_value=float(result.pvalue),
        stderr=float(result.stderr),
        n_points=int(fit_s.size),
    )


def fit_scaling_exponent(*args: object, **kwargs: object) -> FloryFit:
    """Backward-compatible wrapper for :func:`fit_flory_exponent`."""

    return fit_flory_exponent(*args, **kwargs)
