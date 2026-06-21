"""PTM comparison visualization helpers."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idrptm.visualization.heatmaps import (
    RESIDUE_CLASSES,
    plot_heatmap,
    ptm_site_contact_profile,
    residue_class_contact_matrix,
)


def plot_ptm_delta_contact_map(
    wt_contact_map: np.ndarray,
    ptm_contact_map: np.ndarray,
    *,
    title: str = "PTM - WT delta contact map",
) -> plt.Figure:
    """Plot raw PTM-minus-WT contact-map differences."""

    delta = np.asarray(ptm_contact_map, dtype=float) - np.asarray(wt_contact_map, dtype=float)
    return plot_heatmap(
        delta,
        title=title,
        colorbar_label="Delta contact probability (dimensionless)",
        cmap="coolwarm",
        raw=True,
    )


def plot_ptm_site_profile(
    contact_map: np.ndarray,
    ptm_sites_1based: list[int],
    *,
    condition: str = "PTM",
) -> plt.Figure:
    """Plot contacts between PTM sites and all residues."""

    profile = ptm_site_contact_profile(contact_map, ptm_sites_1based)
    fig, ax = plt.subplots(figsize=(6, 3.5))
    ax.plot(profile["residue_index"], profile["contact_probability"], label=condition)
    for site in ptm_sites_1based:
        ax.axvline(site, color="tab:red", linewidth=0.8)
    ax.set_xlabel("Residue index (residue)")
    ax.set_ylabel("PTM-site contact probability (dimensionless)")
    ax.set_title("PTM-site contact profile")
    ax.legend(frameon=False)
    return fig


def delta_ptm_site_profile(
    wt_contact_map: np.ndarray,
    ptm_contact_map: np.ndarray,
    ptm_sites_1based: list[int],
) -> pd.DataFrame:
    """Return raw PTM-site contact-profile delta."""

    wt = ptm_site_contact_profile(wt_contact_map, ptm_sites_1based)
    ptm = ptm_site_contact_profile(ptm_contact_map, ptm_sites_1based)
    merged = wt.merge(ptm, on="residue_index", suffixes=("_wt", "_ptm"))
    merged["delta_contact_probability"] = (
        merged["contact_probability_ptm"] - merged["contact_probability_wt"]
    )
    return merged


def residue_class_contact_changes(
    wt_contact_map: np.ndarray,
    ptm_contact_map: np.ndarray,
    sequence: str,
) -> pd.DataFrame:
    """Return residue-class contact changes from raw maps."""

    wt = residue_class_contact_matrix(wt_contact_map, sequence)
    ptm = residue_class_contact_matrix(ptm_contact_map, sequence)
    delta = ptm - wt
    rows = []
    for class_i in RESIDUE_CLASSES:
        for class_j in RESIDUE_CLASSES:
            rows.append(
                {
                    "class_i": class_i,
                    "class_j": class_j,
                    "delta_contact_probability": float(delta.loc[class_i, class_j]),
                }
            )
    return pd.DataFrame(rows)


def plot_residue_class_contact_changes(delta_table: pd.DataFrame) -> plt.Figure:
    """Plot residue-class contact changes as a heatmap."""

    matrix = delta_table.pivot(
        index="class_i",
        columns="class_j",
        values="delta_contact_probability",
    )
    return plot_heatmap(
        matrix,
        title="Residue-class contact changes",
        colorbar_label="Delta contact probability (dimensionless)",
        cmap="coolwarm",
        x_label="Residue class",
        y_label="Residue class",
        raw=True,
    )
