"""Lightweight smoothing and coarse-graining utilities for analysis curves."""

from __future__ import annotations

from typing import Literal

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy import ndimage, signal, stats


def logspace_smooth_1d(
    x: ArrayLike,
    y: ArrayLike,
    window_log10: float = 0.2,
    min_points: int = 5,
    robust: bool = True,
    center: str = "geometric",
) -> pd.DataFrame:
    """Smooth a one-dimensional curve in log10(x) space without dropping raw rows."""

    x_values, y_values = _as_xy(x, y)
    _validate_positive_x(x_values)
    if window_log10 <= 0:
        raise ValueError("window_log10 must be positive.")
    if min_points <= 0:
        raise ValueError("min_points must be positive.")
    if center not in {"geometric", "arithmetic"}:
        raise ValueError("center must be 'geometric' or 'arithmetic'.")

    log_x = np.log10(x_values)
    half_window = window_log10 / 2.0
    smooth = np.full_like(y_values, np.nan, dtype=float)
    counts = np.zeros_like(y_values, dtype=int)
    finite_y = np.isfinite(y_values)

    for index, value in enumerate(log_x):
        mask = np.abs(log_x - value) <= half_window
        mask &= finite_y
        counts[index] = int(mask.sum())
        if counts[index] < min_points:
            continue
        window_y = y_values[mask]
        smooth[index] = (
            _robust_center(window_y)
            if robust
            else _weighted_mean(log_x[mask], window_y, value)
        )

    return pd.DataFrame(
        {
            "x": x_values,
            "y_raw": y_values,
            "y_smooth": smooth,
            "n_points_window": counts,
            "smoothing_method": "logspace_median" if robust else "logspace_weighted_mean",
            "smoothing_window_log10": window_log10,
        }
    )


def rolling_smooth_1d(
    x: ArrayLike,
    y: ArrayLike,
    window: int,
    method: Literal["mean", "median"] = "mean",
    center: bool = True,
) -> pd.DataFrame:
    """Smooth a time series with a rolling mean or median."""

    x_values, y_values = _as_xy(x, y)
    if window <= 0:
        raise ValueError("window must be positive.")
    series = pd.Series(y_values, dtype=float)
    rolling = series.rolling(window=window, min_periods=1, center=center)
    if method == "mean":
        smooth = rolling.mean()
    elif method == "median":
        smooth = rolling.median()
    else:
        raise ValueError("method must be 'mean' or 'median'.")
    counts = series.rolling(window=window, min_periods=1, center=center).count().astype(int)
    return pd.DataFrame(
        {
            "x": x_values,
            "y_raw": y_values,
            "y_smooth": smooth.to_numpy(dtype=float),
            "n_points_window": counts.to_numpy(dtype=int),
            "smoothing_method": f"rolling_{method}",
            "smoothing_window": window,
        }
    )


def savgol_smooth_1d(
    x: ArrayLike,
    y: ArrayLike,
    window_length: int,
    polyorder: int,
) -> pd.DataFrame:
    """Smooth a time series with a Savitzky-Golay filter."""

    x_values, y_values = _as_xy(x, y)
    if window_length % 2 == 0:
        raise ValueError("window_length must be odd for Savitzky-Golay smoothing.")
    if window_length <= polyorder:
        raise ValueError("window_length must be greater than polyorder.")
    if window_length > y_values.size:
        raise ValueError("window_length cannot exceed the number of samples.")
    finite = pd.Series(y_values, dtype=float).interpolate(limit_direction="both")
    smooth = signal.savgol_filter(finite.to_numpy(dtype=float), window_length, polyorder)
    return pd.DataFrame(
        {
            "x": x_values,
            "y_raw": y_values,
            "y_smooth": smooth,
            "n_points_window": window_length,
            "smoothing_method": "savgol",
            "smoothing_window": window_length,
            "smoothing_polyorder": polyorder,
        }
    )


def bootstrap_smooth_ci(
    x: ArrayLike,
    y: ArrayLike,
    group: ArrayLike | None = None,
    method: str = "logspace",
    n_boot: int = 200,
    ci: float = 0.95,
    seed: int = 123,
    **smoothing_kwargs: object,
) -> pd.DataFrame:
    """Bootstrap confidence intervals for a smoothed curve."""

    x_values, y_values = _as_xy(x, y)
    if n_boot <= 0:
        raise ValueError("n_boot must be positive.")
    if not 0 < ci < 1:
        raise ValueError("ci must be between 0 and 1.")
    rng = np.random.default_rng(seed)
    base = _smooth_dispatch(x_values, y_values, method, smoothing_kwargs)
    boot_values: list[np.ndarray] = []

    if group is None:
        indices = np.arange(x_values.size)
        for _ in range(n_boot):
            sample = rng.choice(indices, size=indices.size, replace=True)
            smoothed = _smooth_dispatch(
                x_values[sample],
                y_values[sample],
                method,
                smoothing_kwargs,
            )
            boot_values.append(smoothed["y_smooth"].to_numpy())
    else:
        group_values = np.asarray(group)
        if group_values.shape != x_values.shape:
            raise ValueError("group must have the same shape as x and y.")
        unique_groups = np.unique(group_values)
        for _ in range(n_boot):
            chosen = rng.choice(unique_groups, size=unique_groups.size, replace=True)
            mask = np.isin(group_values, chosen)
            aggregate = (
                pd.DataFrame({"x": x_values[mask], "y": y_values[mask]})
                .groupby("x", as_index=False)["y"]
                .mean()
                .sort_values("x")
            )
            boot_values.append(
                _smooth_dispatch(aggregate["x"], aggregate["y"], method, smoothing_kwargs)[
                    "y_smooth"
                ].to_numpy()
            )

    stacked = _pad_bootstrap(boot_values, base.shape[0])
    alpha = (1.0 - ci) / 2.0
    result = base[["x", "y_raw", "y_smooth"]].copy()
    result["ci_low"] = np.nanquantile(stacked, alpha, axis=0)
    result["ci_high"] = np.nanquantile(stacked, 1.0 - alpha, axis=0)
    result["n_boot"] = n_boot
    result["ci"] = ci
    return result


def smooth_contact_map(
    contact_map: ArrayLike,
    method: Literal["none", "gaussian", "median"] = "none",
    sigma: float = 1.0,
    size: int = 3,
    preserve_diagonal: bool = True,
) -> np.ndarray:
    """Smooth a contact map for visualization, preserving raw metrics by default."""

    matrix = np.asarray(contact_map, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("contact_map must be a square matrix.")
    if method == "none":
        return matrix.copy()
    diagonal = np.diag(matrix).copy()
    if method == "gaussian":
        if sigma <= 0:
            raise ValueError("sigma must be positive.")
        smoothed = ndimage.gaussian_filter(matrix, sigma=sigma)
    elif method == "median":
        if size <= 0:
            raise ValueError("size must be positive.")
        smoothed = ndimage.median_filter(matrix, size=size)
    else:
        raise ValueError("method must be 'none', 'gaussian', or 'median'.")
    smoothed = (smoothed + smoothed.T) / 2.0
    if preserve_diagonal:
        np.fill_diagonal(smoothed, diagonal)
    return smoothed


def coarse_bin_curve(
    x: ArrayLike,
    y: ArrayLike,
    bins: Literal["log", "linear"],
    n_bins: int,
    agg: Literal["mean", "median"] = "mean",
) -> pd.DataFrame:
    """Aggregate a noisy curve into coarse log or linear bins."""

    x_values, y_values = _as_xy(x, y)
    if n_bins <= 0:
        raise ValueError("n_bins must be positive.")
    if agg not in {"mean", "median"}:
        raise ValueError("agg must be 'mean' or 'median'.")
    if bins == "log":
        _validate_positive_x(x_values)
        edges = np.geomspace(float(np.nanmin(x_values)), float(np.nanmax(x_values)), n_bins + 1)
    elif bins == "linear":
        edges = np.linspace(float(np.nanmin(x_values)), float(np.nanmax(x_values)), n_bins + 1)
    else:
        raise ValueError("bins must be 'log' or 'linear'.")
    table = pd.DataFrame({"x": x_values, "y": y_values})
    table = table[np.isfinite(table["x"]) & np.isfinite(table["y"])].copy()
    table["bin"] = pd.cut(table["x"], bins=np.unique(edges), include_lowest=True)
    grouped = table.groupby("bin", observed=True)
    values = grouped["y"].mean() if agg == "mean" else grouped["y"].median()
    rows = []
    for interval, value in values.items():
        rows.append(
            {
                "bin_start": float(interval.left),
                "bin_end": float(interval.right),
                "x": float(np.sqrt(interval.left * interval.right))
                if bins == "log"
                else float((interval.left + interval.right) / 2.0),
                "y": float(value),
                "n_points": int(grouped.size().loc[interval]),
                "bins": bins,
                "agg": agg,
            }
        )
    return pd.DataFrame(rows)


def _as_xy(x: ArrayLike, y: ArrayLike) -> tuple[np.ndarray, np.ndarray]:
    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if x_values.shape != y_values.shape:
        raise ValueError("x and y must have the same shape.")
    if x_values.ndim != 1:
        raise ValueError("x and y must be one-dimensional.")
    return x_values, y_values


def _validate_positive_x(x_values: np.ndarray) -> None:
    if not np.all(np.isfinite(x_values)):
        raise ValueError("x must contain only finite values.")
    if np.any(x_values <= 0):
        raise ValueError("x must be strictly positive for log-space smoothing.")


def _robust_center(values: np.ndarray) -> float:
    if values.size >= 5:
        return float(stats.trim_mean(values, proportiontocut=0.1))
    return float(np.median(values))


def _weighted_mean(log_x: np.ndarray, y: np.ndarray, center: float) -> float:
    distances = np.abs(log_x - center)
    weights = 1.0 / np.maximum(distances, 1.0e-12)
    return float(np.average(y, weights=weights))


def _smooth_dispatch(
    x: ArrayLike,
    y: ArrayLike,
    method: str,
    kwargs: dict[str, object],
) -> pd.DataFrame:
    if method == "logspace":
        return logspace_smooth_1d(x, y, **kwargs)
    if method == "rolling":
        return rolling_smooth_1d(x, y, **kwargs)
    if method == "savgol":
        return savgol_smooth_1d(x, y, **kwargs)
    raise ValueError("method must be 'logspace', 'rolling', or 'savgol'.")


def _pad_bootstrap(values: list[np.ndarray], width: int) -> np.ndarray:
    padded = np.full((len(values), width), np.nan, dtype=float)
    for row_index, row in enumerate(values):
        count = min(width, row.size)
        padded[row_index, :count] = row[:count]
    return padded
