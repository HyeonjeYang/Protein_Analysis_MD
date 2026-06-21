"""Unified heatmap plotting helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from numpy.typing import ArrayLike

RESIDUE_CLASSES: dict[str, set[str]] = {
    "positive": {"K", "R", "H"},
    "negative": {"D", "E", "B", "O"},
    "aromatic": {"F", "Y", "W"},
    "polar": {"S", "T", "N", "Q", "C"},
    "hydrophobic": {"A", "V", "I", "L", "M"},
    "proline_glycine": {"P", "G"},
}


def expected_by_sequence_separation(contact_map: ArrayLike) -> np.ndarray:
    """Return expected contact probability for each sequence separation."""

    matrix = _as_square_matrix(contact_map)
    expected = np.zeros(matrix.shape[0], dtype=float)
    for separation in range(matrix.shape[0]):
        values = np.diag(matrix, k=separation)
        expected[separation] = float(np.nanmean(values)) if values.size else 0.0
    return expected


def observed_expected_contact_map(
    contact_map: ArrayLike,
    *,
    eps: float = 1.0e-6,
    min_sequence_separation: int = 1,
    method: str = "log_ratio",
) -> np.ndarray:
    """Compute an observed/expected-like contact map by sequence separation."""

    matrix = _as_square_matrix(contact_map)
    expected = expected_by_sequence_separation(matrix)
    residue_i, residue_j = np.indices(matrix.shape)
    expected_matrix = expected[np.abs(residue_i - residue_j)]
    if method == "log_ratio":
        oe = np.log((matrix + eps) / (expected_matrix + eps))
    elif method == "difference":
        oe = matrix - expected_matrix
    else:
        raise ValueError("method must be 'log_ratio' or 'difference'.")
    oe[np.abs(residue_i - residue_j) < min_sequence_separation] = np.nan
    return oe


def residue_class_contact_matrix(
    contact_map: ArrayLike,
    sequence: str,
    *,
    classes: dict[str, set[str]] | None = None,
) -> pd.DataFrame:
    """Aggregate contact probabilities by residue-class pair."""

    matrix = _as_square_matrix(contact_map)
    if len(sequence) != matrix.shape[0]:
        raise ValueError("sequence length must match contact_map dimensions.")
    class_map = classes or RESIDUE_CLASSES
    labels = list(class_map)
    output = pd.DataFrame(0.0, index=labels, columns=labels)
    counts = pd.DataFrame(0, index=labels, columns=labels)
    residue_labels = [_class_for_residue(residue, class_map) for residue in sequence]
    for i, class_i in enumerate(residue_labels):
        for j, class_j in enumerate(residue_labels):
            if class_i is None or class_j is None or i == j:
                continue
            output.loc[class_i, class_j] += matrix[i, j]
            counts.loc[class_i, class_j] += 1
    return output.divide(counts.replace(0, np.nan)).fillna(0.0)


def ptm_site_contact_profile(contact_map: ArrayLike, ptm_sites_1based: list[int]) -> pd.DataFrame:
    """Return contact probability from PTM sites to each residue."""

    matrix = _as_square_matrix(contact_map)
    if not ptm_sites_1based:
        return pd.DataFrame(columns=["residue_index", "contact_probability"])
    indices = np.array([site - 1 for site in ptm_sites_1based], dtype=int)
    if np.any(indices < 0) or np.any(indices >= matrix.shape[0]):
        raise ValueError("PTM site is outside the contact map.")
    profile = matrix[indices, :].mean(axis=0)
    return pd.DataFrame(
        {
            "residue_index": np.arange(1, matrix.shape[0] + 1, dtype=int),
            "contact_probability": profile,
        }
    )


def plot_heatmap(
    matrix: ArrayLike | pd.DataFrame,
    *,
    title: str,
    colorbar_label: str,
    x_label: str = "Residue index (residue)",
    y_label: str = "Residue index (residue)",
    cmap: str = "viridis",
    ptm_sites: list[int] | None = None,
    cleavage_sites: list[int] | None = None,
    fragment_boundaries: list[int] | None = None,
    raw: bool = True,
) -> plt.Figure:
    """Plot a heatmap with optional PTM, cleavage, and fragment annotations."""

    table_labels = None
    if isinstance(matrix, pd.DataFrame):
        values = matrix.to_numpy(dtype=float)
        table_labels = (matrix.index.astype(str).tolist(), matrix.columns.astype(str).tolist())
    else:
        values = np.asarray(matrix, dtype=float)
    if values.ndim != 2:
        raise ValueError("matrix must be two-dimensional.")
    fig, ax = plt.subplots(figsize=(5.2, 4.6))
    image = ax.imshow(values, origin="lower", cmap=cmap)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    suffix = "raw" if raw else "display-smoothed"
    ax.set_title(f"{title} ({suffix})")
    if table_labels is not None:
        rows, columns = table_labels
        ax.set_yticks(range(len(rows)), rows)
        ax.set_xticks(range(len(columns)), columns, rotation=30, ha="right")
    _annotate_sites(ax, ptm_sites or [], color="tab:red", linestyle="-", label="PTM")
    _annotate_sites(ax, cleavage_sites or [], color="tab:orange", linestyle="--", label="cleavage")
    for boundary in fragment_boundaries or []:
        ax.axvline(boundary - 0.5, color="white", linewidth=0.8, linestyle=":")
        ax.axhline(boundary - 0.5, color="white", linewidth=0.8, linestyle=":")
    fig.colorbar(image, ax=ax, label=colorbar_label)
    return fig


def save_matrix_data(matrix: ArrayLike | pd.DataFrame, output: str | Path) -> Path:
    """Save the raw matrix behind a heatmap."""

    path = Path(output)
    path.parent.mkdir(parents=True, exist_ok=True)
    if isinstance(matrix, pd.DataFrame):
        matrix.to_csv(path.with_suffix(".csv"))
        return path.with_suffix(".csv")
    np.save(path.with_suffix(".npy"), np.asarray(matrix, dtype=float))
    return path.with_suffix(".npy")


def _annotate_sites(
    ax: plt.Axes,
    sites_1based: list[int],
    *,
    color: str,
    linestyle: str,
    label: str,
) -> None:
    for index, site in enumerate(sites_1based):
        label_text = label if index == 0 else None
        ax.axvline(site - 1, color=color, linewidth=0.8, linestyle=linestyle, label=label_text)
        ax.axhline(site - 1, color=color, linewidth=0.8, linestyle=linestyle)
    if sites_1based:
        ax.legend(frameon=False, loc="upper right")


def _as_square_matrix(values: ArrayLike) -> np.ndarray:
    matrix = np.asarray(values, dtype=float)
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("contact map must be a square matrix.")
    return matrix


def _class_for_residue(residue: str, classes: dict[str, set[str]]) -> str | None:
    code = residue.upper()
    for label, members in classes.items():
        if code in members:
            return label
    return None
