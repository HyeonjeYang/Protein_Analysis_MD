"""Visualization helpers for exploratory decomposition outputs."""

from __future__ import annotations

import pandas as pd

from idrptm.plotting.plots import (
    plot_contact_eigenvectors,
    plot_contact_loading_heatmap,
    plot_delta_ev,
    plot_ev1_correlation,
    plot_explained_variance,
    plot_pca_centroid_shift,
    plot_pca_score_scatter,
    plot_pca_timeseries,
)


def decomposition_caption() -> str:
    """Return cautious terminology for reports."""

    return (
        "Contact-environment eigenvectors and PCA panels are exploratory, "
        "compartment-like contact decompositions; they are not chromosome compartments."
    )


def pc_score_distribution_shift(
    scores: pd.DataFrame,
    condition_column: str = "condition",
) -> pd.DataFrame:
    """Summarize PC-score centroid and spread by condition."""

    pc_columns = [column for column in scores.columns if column.startswith("PC")]
    rows = []
    for condition, group in scores.groupby(condition_column, sort=True):
        row: dict[str, float | str | int] = {
            "condition": str(condition),
            "n_frames": int(len(group)),
        }
        for column in pc_columns:
            row[f"{column}_mean"] = float(group[column].mean())
            row[f"{column}_std"] = float(group[column].std(ddof=1)) if len(group) > 1 else 0.0
        rows.append(row)
    return pd.DataFrame(rows)


__all__ = [
    "decomposition_caption",
    "pc_score_distribution_shift",
    "plot_contact_eigenvectors",
    "plot_contact_loading_heatmap",
    "plot_delta_ev",
    "plot_ev1_correlation",
    "plot_explained_variance",
    "plot_pca_centroid_shift",
    "plot_pca_score_scatter",
    "plot_pca_timeseries",
]
