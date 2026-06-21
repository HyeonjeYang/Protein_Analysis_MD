"""Multi-protein and phase-separation visualization helpers."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idrptm.plotting.plots import plot_cluster_size_distribution
from idrptm.visualization.heatmaps import plot_heatmap


def plot_density_profile(
    table: pd.DataFrame,
    *,
    z_column: str = "z_nm",
    density_column: str = "density",
    component_column: str | None = "component",
    smooth_column: str | None = None,
    dense_region: tuple[float, float] | None = None,
) -> plt.Figure:
    """Plot slab or droplet-axis density profile with optional visual smoothing."""

    fig, ax = plt.subplots(figsize=(6, 4))
    groupby = component_column if component_column in table else None
    groups = table.groupby(groupby, sort=True) if groupby else [("total", table)]
    for label, group in groups:
        ax.plot(
            group[z_column],
            group[density_column],
            alpha=0.45,
            marker="o",
            label=f"{label} raw",
        )
        if smooth_column and smooth_column in group:
            ax.plot(group[z_column], group[smooth_column], linewidth=2, label=f"{label} smoothed")
    if dense_region is not None:
        ax.axvspan(dense_region[0], dense_region[1], color="0.85", alpha=0.4, label="dense region")
    ax.set_xlabel("z (nm)")
    ax.set_ylabel("Concentration / density (a.u.)")
    ax.set_title("Density profile (phase behavior exploratory)")
    ax.legend(frameon=False)
    return fig


def plot_dense_dilute_timeseries(table: pd.DataFrame) -> plt.Figure:
    """Plot dense/dilute concentration time series."""

    fig, ax = plt.subplots(figsize=(6, 4))
    x_column = "time_ns" if "time_ns" in table else "frame"
    for column in ("dense_concentration", "dilute_concentration"):
        if column in table:
            ax.plot(table[x_column], table[column], label=column.replace("_", " "))
    ax.set_xlabel("Time (ns)" if x_column == "time_ns" else "Frame (index)")
    ax.set_ylabel("Concentration (a.u.)")
    ax.set_title("Dense/dilute concentration time series")
    ax.legend(frameon=False)
    return fig


def plot_density_projection(
    positions_xy: np.ndarray,
    *,
    bins: int = 50,
    plane: str = "xy",
    display_smoothed: bool = False,
) -> plt.Figure:
    """Plot a 2D projected density heatmap from bead/COM positions."""

    values = np.asarray(positions_xy, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("positions_xy must have shape (n_points, 2).")
    heatmap, x_edges, y_edges = np.histogram2d(values[:, 0], values[:, 1], bins=bins)
    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    mesh = ax.pcolormesh(x_edges, y_edges, heatmap.T, shading="auto", cmap="magma")
    ax.set_xlabel(f"{plane[0]} (nm)")
    ax.set_ylabel(f"{plane[1]} (nm)")
    suffix = "display-smoothed" if display_smoothed else "raw binned"
    ax.set_title(f"{plane.upper()} density projection ({suffix})")
    fig.colorbar(mesh, ax=ax, label="Count")
    return fig


def plot_partition_coefficients(table: pd.DataFrame) -> plt.Figure:
    """Plot component partition coefficients K = dense/dilute."""

    plot_table = table.copy()
    if "partition_coefficient" not in plot_table:
        plot_table["partition_coefficient"] = (
            plot_table["dense_concentration"] / plot_table["dilute_concentration"]
        )
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.bar(plot_table["component"].astype(str), plot_table["partition_coefficient"])
    ax.set_xlabel("Component")
    ax.set_ylabel("K = dense / dilute")
    ax.set_title("Component partitioning")
    ax.tick_params(axis="x", rotation=25)
    return fig


def plot_inter_chain_contact_heatmap(matrix: np.ndarray) -> plt.Figure:
    """Plot inter-chain contact matrix."""

    return plot_heatmap(
        matrix,
        title="Inter-chain contact matrix",
        colorbar_label="Contact probability (dimensionless)",
        raw=True,
    )


def plot_homotypic_heterotypic_heatmap(matrix: pd.DataFrame) -> plt.Figure:
    """Plot homotypic/heterotypic contact summary matrix."""

    return plot_heatmap(
        matrix,
        title="Homotypic/heterotypic contacts",
        colorbar_label="Contact probability (dimensionless)",
        x_label="Component",
        y_label="Component",
        raw=True,
    )


def phase_reliability_warning(n_frames: int, n_chains: int) -> str | None:
    """Return a cautionary warning for short/small phase-separation analyses."""

    if n_frames < 100 or n_chains < 10:
        return (
            "Phase-separation visualization is exploratory: trajectory is short or system is small "
            "for reliable dense/dilute inference."
        )
    return None


def cluster_distribution_figure(table: pd.DataFrame) -> plt.Figure:
    """Wrapper for cluster-size distribution plot."""

    return plot_cluster_size_distribution(table)
