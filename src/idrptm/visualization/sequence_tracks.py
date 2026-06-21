"""Sequence-track visualization and per-sequence summary panels."""

from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def sequence_feature_table(sequence: str) -> pd.DataFrame:
    """Return per-residue sequence features used by visualization tracks."""

    rows = []
    for index, residue in enumerate(sequence, start=1):
        rows.append(
            {
                "residue_index": index,
                "residue": residue,
                "class": residue_class(residue),
                "charge": residue_charge(residue),
                "hydropathy": hydropathy(residue),
                "aromatic": float(residue.upper() in {"F", "Y", "W"}),
            }
        )
    return pd.DataFrame(rows)


def sequence_summary(sequence: str) -> dict[str, float]:
    """Return simple IDR sequence summary features."""

    if not sequence:
        return {
            "length": 0.0,
            "net_charge": 0.0,
            "NCPR": 0.0,
            "FCR": 0.0,
            "aromatic_fraction": 0.0,
        }
    charges = np.array([residue_charge(residue) for residue in sequence], dtype=float)
    charged = np.array([abs(value) > 0 for value in charges], dtype=float)
    aromatic = np.array([residue.upper() in {"F", "Y", "W"} for residue in sequence], dtype=float)
    return {
        "length": float(len(sequence)),
        "net_charge": float(charges.sum()),
        "NCPR": float(charges.sum() / len(sequence)),
        "FCR": float(charged.mean()),
        "aromatic_fraction": float(aromatic.mean()),
    }


def plot_sequence_tracks(
    sequence: str,
    *,
    ptm_sites: list[int] | None = None,
    cleavage_sites: list[int] | None = None,
    fragment_ranges: list[tuple[int, int, str]] | None = None,
    title: str = "Sequence tracks",
) -> plt.Figure:
    """Plot residue class, charge, hydropathy, and sequence annotations."""

    features = sequence_feature_table(sequence)
    fig, axes = plt.subplots(4, 1, figsize=(8, 5.6), sharex=True)
    class_codes = [CLASS_ORDER.get(value, 0) for value in features["class"]]
    axes[0].imshow(np.array([class_codes]), aspect="auto", cmap="tab20")
    axes[0].set_yticks([])
    axes[0].set_ylabel("Class")
    axes[1].plot(features["residue_index"], features["charge"], color="tab:blue")
    axes[1].axhline(0, color="0.7", linewidth=0.8)
    axes[1].set_ylabel("Charge (e)")
    axes[2].plot(features["residue_index"], features["hydropathy"], color="tab:green")
    axes[2].set_ylabel("Hydropathy")
    axes[3].hlines(0, 1, max(len(sequence), 1), color="0.75", linewidth=5)
    for start, end, label in fragment_ranges or []:
        axes[3].hlines(0.15, start, end, linewidth=5)
        axes[3].text((start + end) / 2, 0.32, label, ha="center", va="bottom", fontsize=7)
    for ax in axes:
        for site in ptm_sites or []:
            ax.axvline(site, color="tab:red", linewidth=0.9)
        for site in cleavage_sites or []:
            ax.axvline(site, color="tab:orange", linewidth=0.9, linestyle="--")
    axes[3].set_yticks([])
    axes[3].set_ylabel("Fragments")
    axes[3].set_xlabel("Residue index (residue)")
    axes[0].set_title(title)
    return fig


def plot_sequence_summary(sequence: str) -> plt.Figure:
    """Render a compact pI-like charge/composition summary panel."""

    summary = sequence_summary(sequence)
    labels = ["net_charge", "NCPR", "FCR", "aromatic_fraction"]
    values = [summary[label] for label in labels]
    fig, ax = plt.subplots(figsize=(5.5, 3.2))
    ax.bar(labels, values)
    ax.set_ylabel("Value")
    ax.set_title("Sequence summary")
    ax.tick_params(axis="x", rotation=25)
    return fig


def residue_class(residue: str) -> str:
    """Return a coarse residue class label."""

    code = residue.upper()
    if code in {"K", "R", "H"}:
        return "positive"
    if code in {"D", "E", "B", "O"}:
        return "negative"
    if code in {"F", "Y", "W"}:
        return "aromatic"
    if code in {"S", "T", "N", "Q", "C"}:
        return "polar"
    if code in {"A", "V", "I", "L", "M"}:
        return "hydrophobic"
    if code in {"P", "G"}:
        return "proline_glycine"
    return "other"


def residue_charge(residue: str) -> float:
    """Simple charge proxy for sequence-track visualization."""

    return {"D": -1.0, "E": -1.0, "K": 1.0, "R": 1.0, "H": 0.1, "B": -2.0, "O": -2.0}.get(
        residue.upper(),
        0.0,
    )


def hydropathy(residue: str) -> float:
    """Kyte-Doolittle hydropathy proxy."""

    values = {
        "I": 4.5,
        "V": 4.2,
        "L": 3.8,
        "F": 2.8,
        "C": 2.5,
        "M": 1.9,
        "A": 1.8,
        "G": -0.4,
        "T": -0.7,
        "S": -0.8,
        "W": -0.9,
        "Y": -1.3,
        "P": -1.6,
        "H": -3.2,
        "E": -3.5,
        "Q": -3.5,
        "D": -3.5,
        "N": -3.5,
        "K": -3.9,
        "R": -4.5,
    }
    return values.get(residue.upper(), 0.0)


CLASS_ORDER = {
    "positive": 1,
    "negative": 2,
    "aromatic": 3,
    "polar": 4,
    "hydrophobic": 5,
    "proline_glycine": 6,
    "other": 7,
}
