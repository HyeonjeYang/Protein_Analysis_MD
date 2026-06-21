"""Single-chain protein/IDR visualization panels."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idrptm.visualization.heatmaps import (
    observed_expected_contact_map,
    plot_heatmap,
    ptm_site_contact_profile,
)


def plot_rg_ree_hexbin(
    rg: pd.DataFrame,
    ree: pd.DataFrame,
    *,
    title: str = "Rg/Ree joint density",
) -> plt.Figure:
    """Plot a joint Rg-vs-Ree hexbin density."""

    table = _merge_frame_tables(rg, ree, "rg", "ree")
    fig, ax = plt.subplots(figsize=(5, 4.2))
    image = ax.hexbin(table["rg"], table["ree"], gridsize=30, mincnt=1, cmap="viridis")
    ax.set_xlabel("Rg (nm)")
    ax.set_ylabel("Ree (nm)")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Frame count")
    return fig


def plot_rg_ree_timeseries(
    rg: pd.DataFrame,
    ree: pd.DataFrame,
    *,
    show_smoothed: bool = True,
) -> plt.Figure:
    """Plot Rg and Ree time series with optional smoothed visual trends."""

    table = _merge_frame_tables(rg, ree, "rg", "ree")
    x_column = "time_ns" if "time_ns" in table else "frame"
    if "time_ps" in table and "time_ns" not in table:
        table["time_ns"] = table["time_ps"] / 1000.0
        x_column = "time_ns"
    fig, ax = plt.subplots(figsize=(6.2, 4))
    for column, label in (("rg", "Rg raw"), ("ree", "Ree raw")):
        ax.plot(table[x_column], table[column], alpha=0.45, label=label)
    if show_smoothed:
        for column, label in (
            ("rg_nm_smooth", "Rg smoothed visual trend"),
            ("ree_nm_smooth", "Ree smoothed visual trend"),
        ):
            if column in table:
                ax.plot(table[x_column], table[column], linewidth=2, label=label)
    ax.set_xlabel("Time (ns)" if x_column == "time_ns" else "Frame (index)")
    ax.set_ylabel("Distance (nm)")
    ax.set_title("Rg/Ree time series")
    ax.legend(frameon=False)
    return fig


def local_scaling_exponent(
    scaling: pd.DataFrame,
    *,
    distance_column: str = "mean_distance_nm",
) -> pd.DataFrame:
    """Compute local ``d log R(s) / d log s`` for visualization."""

    if distance_column not in scaling and "distance" in scaling:
        distance_column = "distance"
    table = scaling[["s", distance_column]].dropna().copy()
    table = table[(table["s"] > 0) & (table[distance_column] > 0)]
    log_s = np.log(table["s"].to_numpy(dtype=float))
    log_r = np.log(table[distance_column].to_numpy(dtype=float))
    table["local_nu"] = np.gradient(log_r, log_s) if len(table) > 1 else np.nan
    return table[["s", "local_nu"]]


def plot_local_scaling_exponent(local: pd.DataFrame) -> plt.Figure:
    """Plot the local scaling exponent as a visual diagnostic."""

    fig, ax = plt.subplots(figsize=(5.5, 3.8))
    ax.plot(local["s"], local["local_nu"], marker="o")
    ax.set_xscale("log")
    ax.set_xlabel("Sequence separation s (residues)")
    ax.set_ylabel("Local scaling exponent nu(s)")
    ax.set_title("Local scaling exponent (visual diagnostic)")
    return fig


def plot_contact_degree(contact_map: np.ndarray) -> plt.Figure:
    """Plot per-residue contact degree."""

    degree = np.asarray(contact_map, dtype=float).sum(axis=1)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(np.arange(1, degree.size + 1), degree)
    ax.set_xlabel("Residue index (residue)")
    ax.set_ylabel("Contact degree (dimensionless)")
    ax.set_title("Per-residue contact degree")
    return fig


def plot_observed_expected_contact_map(contact_map: np.ndarray) -> plt.Figure:
    """Plot an observed/expected-like contact map."""

    oe = observed_expected_contact_map(contact_map)
    return plot_heatmap(
        oe,
        title="Observed/expected-like contact map",
        colorbar_label="log contact enrichment (a.u.)",
        raw=True,
    )


def plot_ptm_centered_contact_profile(
    contact_map: np.ndarray,
    ptm_sites_1based: list[int],
    *,
    title: str = "PTM-centered contact profile",
) -> plt.Figure:
    """Plot contacts involving configured PTM sites."""

    profile = ptm_site_contact_profile(contact_map, ptm_sites_1based)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(profile["residue_index"], profile["contact_probability"])
    for site in ptm_sites_1based:
        ax.axvline(site, color="tab:red", linewidth=0.8)
    ax.set_xlabel("Residue index (residue)")
    ax.set_ylabel("Contact probability (dimensionless)")
    ax.set_title(title)
    return fig


def rg_ree_plotting_data(rg: pd.DataFrame, ree: pd.DataFrame) -> pd.DataFrame:
    """Return merged Rg/Ree plotting data for artifact saving."""

    return _merge_frame_tables(rg, ree, "rg", "ree")


def _merge_frame_tables(
    left: pd.DataFrame,
    right: pd.DataFrame,
    left_value: str,
    right_value: str,
) -> pd.DataFrame:
    keys = ["frame"]
    if "time_ps" in left and "time_ps" in right:
        keys.append("time_ps")
    table = left[keys + [left_value]].merge(right[keys + [right_value]], on=keys, how="inner")
    if "time_ps" in table:
        table["time_ns"] = table["time_ps"] / 1000.0
    for column in ("rg_nm_smooth", "ree_nm_smooth"):
        source = left if column.startswith("rg") else right
        if column in source:
            table = table.merge(source[["frame", column]], on="frame", how="left")
    return table
