"""Figure-generation helpers for reports."""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save_figure(fig: plt.Figure, output_base: str | Path) -> tuple[Path, Path]:
    """Save a Matplotlib figure as PNG and PDF."""

    base = Path(output_base)
    base.parent.mkdir(parents=True, exist_ok=True)
    png = base.with_suffix(".png")
    pdf = base.with_suffix(".pdf")
    fig.savefig(png, dpi=200, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)
    return png, pdf


def plot_distribution(
    table: pd.DataFrame,
    value_column: str,
    ylabel: str,
    title: str,
) -> plt.Figure:
    """Plot per-condition distributions as boxplots with replicate points."""

    conditions = sorted(table["condition"].unique())
    values = [
        table.loc[table["condition"] == condition, value_column].to_numpy()
        for condition in conditions
    ]
    fig, ax = plt.subplots(figsize=(max(5, 1.3 * len(conditions)), 4))
    ax.boxplot(values, tick_labels=conditions, showfliers=False)
    for index, condition_values in enumerate(values, start=1):
        jitter = np.linspace(-0.08, 0.08, len(condition_values)) if len(condition_values) else []
        ax.scatter(
            np.full(len(condition_values), index) + jitter,
            condition_values,
            s=16,
            alpha=0.65,
        )
    ax.set_ylabel(ylabel)
    ax.set_xlabel("Condition")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=30)
    return fig


def plot_matrix(matrix: np.ndarray, title: str, label: str, cmap: str = "viridis") -> plt.Figure:
    """Plot a square matrix with a colorbar."""

    fig, ax = plt.subplots(figsize=(5, 4.5))
    image = ax.imshow(matrix, origin="lower", cmap=cmap)
    ax.set_xlabel("Residue index (residue)")
    ax.set_ylabel("Residue index (residue)")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label=label)
    return fig


def plot_delta_matrix(matrix: np.ndarray, title: str) -> plt.Figure:
    """Plot a signed delta contact map."""

    vmax = float(np.nanmax(np.abs(matrix))) if matrix.size else 1.0
    vmax = vmax or 1.0
    fig, ax = plt.subplots(figsize=(5, 4.5))
    image = ax.imshow(matrix, origin="lower", cmap="coolwarm", vmin=-vmax, vmax=vmax)
    ax.set_xlabel("Residue index (residue)")
    ax.set_ylabel("Residue index (residue)")
    ax.set_title(title)
    fig.colorbar(image, ax=ax, label="Delta contact probability (dimensionless)")
    return fig


def plot_lines(
    table: pd.DataFrame,
    x: str,
    y: str,
    ylabel: str,
    title: str,
) -> plt.Figure:
    """Plot one line per condition."""

    fig, ax = plt.subplots(figsize=(5.5, 4))
    for condition, group in table.groupby("condition", sort=True):
        ax.plot(group[x], group[y], marker="o", label=condition)
    ax.set_xlabel(_axis_label(x))
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(frameon=False)
    return fig


def plot_ptm_sites(manifest: pd.DataFrame) -> plt.Figure:
    """Plot PTM site positions by condition."""

    fig, ax = plt.subplots(figsize=(7, max(2.5, 0.45 * len(manifest))))
    ytick_labels = []
    for row_index, row in manifest.reset_index(drop=True).iterrows():
        sequence = str(row.get("original_sequence", ""))
        length = len(sequence)
        condition = str(row.get("ptm_state") or row.get("condition") or row.get("variant_id"))
        ytick_labels.append(condition)
        ax.hlines(row_index, 1, max(length, 1), color="0.85", linewidth=4)
        for site in _parse_ptm_sites(str(row.get("ptm_sites_1based", ""))):
            ax.scatter(site, row_index, color="tab:red", s=45, zorder=3)
    ax.set_xlabel("Biological residue position (residue)")
    ax.set_yticks(range(len(ytick_labels)), ytick_labels)
    ax.set_title("PTM site annotation")
    ax.set_ylim(-0.8, len(ytick_labels) - 0.2)
    return fig


def plot_inter_protein_contact_map(matrix: np.ndarray, title: str) -> plt.Figure:
    """Plot an inter-protein residue contact probability map."""

    return plot_matrix(matrix, title, "Inter-protein contact probability (dimensionless)")


def plot_com_distance_distribution(table: pd.DataFrame) -> plt.Figure:
    """Plot COM distance distributions for chain pairs."""

    pairs = sorted({f"{row.chain_i}-{row.chain_j}" for row in table.itertuples()})
    values = []
    for pair in pairs:
        chain_i, chain_j = pair.split("-", maxsplit=1)
        values.append(
            table.loc[
                (table["chain_i"] == chain_i) & (table["chain_j"] == chain_j),
                "distance",
            ].to_numpy()
        )
    fig, ax = plt.subplots(figsize=(max(5, 1.4 * len(pairs)), 4))
    ax.boxplot(values, tick_labels=pairs, showfliers=False)
    ax.set_xlabel("Chain pair")
    ax.set_ylabel("COM distance (nm)")
    ax.set_title("COM distance distribution")
    ax.tick_params(axis="x", rotation=30)
    return fig


def plot_cluster_size_distribution(table: pd.DataFrame) -> plt.Figure:
    """Plot the largest cluster size distribution."""

    fig, ax = plt.subplots(figsize=(5, 4))
    bins = np.arange(
        int(table["largest_cluster_size"].min()),
        int(table["largest_cluster_size"].max()) + 2,
    )
    ax.hist(table["largest_cluster_size"], bins=bins, align="left", rwidth=0.8)
    ax.set_xlabel("Largest cluster size (chains)")
    ax.set_ylabel("Frame count")
    ax.set_title("Cluster size distribution")
    return fig


def plot_chain_resolved_rg(table: pd.DataFrame) -> plt.Figure:
    """Plot chain-resolved Rg distributions."""

    plot_table = table.rename(columns={"chain_id": "condition"})
    return plot_distribution(plot_table, "rg", "Rg (nm)", "Chain-resolved Rg")


def plot_cleavage_map(
    sequence_length: int,
    cleavage_sites: pd.DataFrame,
    title: str = "Cleavage map",
) -> plt.Figure:
    """Plot cleavage sites along the original sequence."""

    fig, ax = plt.subplots(figsize=(8, 1.8))
    ax.hlines(0, 1, sequence_length, color="0.75", linewidth=6)
    if not cleavage_sites.empty:
        ax.vlines(
            cleavage_sites["cut_after"],
            -0.35,
            0.35,
            color="tab:red",
            linewidth=2,
        )
    ax.set_xlim(1, max(sequence_length, 1))
    ax.set_ylim(-0.8, 0.8)
    ax.set_yticks([])
    ax.set_xlabel("Original residue position (residue)")
    ax.set_title(title)
    return fig


def plot_fragment_architecture(
    fragments: pd.DataFrame,
    title: str = "Fragment architecture",
) -> plt.Figure:
    """Plot fragment original-coordinate ranges."""

    fig, ax = plt.subplots(figsize=(8, max(2.0, 0.35 * len(fragments) + 1.0)))
    for row_index, row in fragments.reset_index(drop=True).iterrows():
        ax.hlines(
            row_index,
            int(row["original_start"]),
            int(row["original_end"]),
            linewidth=5,
            color="tab:blue",
        )
    ax.set_yticks(range(len(fragments)), fragments["fragment_id"].astype(str).tolist())
    ax.set_xlabel("Original residue position (residue)")
    ax.set_title(title)
    return fig


def plot_intact_vs_cleaved_contact_map(
    intact: np.ndarray,
    cleaved: np.ndarray,
    title: str = "Intact vs cleaved contact maps",
) -> plt.Figure:
    """Plot intact and cleaved contact maps side by side."""

    fig, axes = plt.subplots(1, 2, figsize=(9, 4))
    vmax = max(float(np.nanmax(intact)), float(np.nanmax(cleaved)), 1.0)
    for ax, matrix, label in zip(axes, [intact, cleaved], ["Intact", "Cleaved"], strict=True):
        image = ax.imshow(matrix, origin="lower", cmap="viridis", vmin=0, vmax=vmax)
        ax.set_xlabel("Residue index (residue)")
        ax.set_ylabel("Residue index (residue)")
        ax.set_title(label)
    fig.colorbar(image, ax=axes.ravel().tolist(), label="Contact probability (dimensionless)")
    fig.suptitle(title)
    return fig


def plot_cut_number_summary(
    table: pd.DataFrame,
    y: str,
    ylabel: str,
    title: str,
) -> plt.Figure:
    """Plot a summary observable against cleavage cut number."""

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.plot(table["cut_number"], table[y], marker="o")
    ax.set_xlabel("Cut number (count)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    return fig


def plot_energy_timeseries(table: pd.DataFrame) -> plt.Figure:
    """Plot potential energy and optional temperature over time."""

    fig, ax = plt.subplots(figsize=(6, 4))
    if "potential_energy_kj_mol" in table:
        ax.plot(table["time_ns"], table["potential_energy_kj_mol"], label="Potential energy")
    if "total_energy_kj_mol" in table:
        ax.plot(table["time_ns"], table["total_energy_kj_mol"], label="Total energy")
    ax.set_xlabel("Time (ns)")
    ax.set_ylabel("Energy (kJ/mol)")
    ax.set_title("Energy timeseries")
    ax.legend(frameon=False)
    return fig


def plot_event_schedule(events: pd.DataFrame) -> plt.Figure:
    """Plot cleavage event time versus cut number."""

    fig, ax = plt.subplots(figsize=(5.5, 4))
    ax.step(events["event_time_ns"], events["cut_number"], where="post")
    ax.set_xlabel("Event time (ns)")
    ax.set_ylabel("Cut number (count)")
    ax.set_title("Cleavage event schedule")
    return fig


def plot_summary_table(summary: pd.DataFrame) -> plt.Figure:
    """Render the scalar comparison summary as a compact table."""

    columns = [
        "condition",
        "n_replicates",
        "delta_mean_Rg",
        "delta_mean_Ree",
        "delta_flory_exponent",
    ]
    display = summary[columns].copy()
    for column in columns[2:]:
        display[column] = display[column].map(lambda value: f"{value:.4g}")
    fig, ax = plt.subplots(figsize=(9, max(2.0, 0.45 * len(display) + 1.0)))
    ax.axis("off")
    table = ax.table(
        cellText=display.values,
        colLabels=display.columns,
        loc="center",
        cellLoc="center",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.25)
    ax.set_title("WT-vs-PTM summary")
    return fig


def _parse_ptm_sites(site_text: str) -> list[int]:
    positions: list[int] = []
    for match in site_text.split(";"):
        if not match:
            continue
        digits = "".join(character for character in match if character.isdigit())
        if digits:
            positions.append(int(digits))
    return positions


def _axis_label(column: str) -> str:
    labels = {
        "s": "Sequence separation s (residues)",
        "lag": "Lag (frames)",
        "frame": "Frame (index)",
    }
    return labels.get(column, column)
