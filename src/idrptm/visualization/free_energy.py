"""Free-energy landscape plotting helpers."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from idrptm.analysis.free_energy import FreeEnergySurface


def plot_free_energy_surface(
    surface: FreeEnergySurface,
    *,
    title: str = "Free-energy landscape",
) -> plt.Figure:
    """Plot a 2D exploratory free-energy surface from a computed grid."""

    return plot_free_energy_grid(
        surface.free_energy_kj_mol,
        surface.x_edges,
        surface.y_edges,
        x_label=str(surface.metadata.get("x_label", "x")),
        y_label=str(surface.metadata.get("y_label", "y")),
        title=title,
    )


def plot_free_energy_grid(
    free_energy_kj_mol: np.ndarray,
    x_edges: np.ndarray,
    y_edges: np.ndarray,
    *,
    x_label: str,
    y_label: str,
    title: str = "Free-energy landscape",
) -> plt.Figure:
    """Plot a free-energy grid loaded from analysis outputs."""

    fig, ax = plt.subplots(figsize=(5.2, 4.4))
    mesh = ax.pcolormesh(x_edges, y_edges, free_energy_kj_mol.T, shading="auto", cmap="magma")
    ax.set_xlabel(_axis_label(x_label))
    ax.set_ylabel(_axis_label(y_label))
    ax.set_title(f"{title} (exploratory)")
    fig.colorbar(mesh, ax=ax, label="F (kJ/mol)")
    return fig


def load_free_energy_grid(prefix: str | Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Load free-energy grid, x edges, and y edges from a prefix."""

    base = Path(prefix)
    return (
        np.load(base.with_name(f"{base.name}_free_energy.npy")),
        np.load(base.with_name(f"{base.name}_x_edges.npy")),
        np.load(base.with_name(f"{base.name}_y_edges.npy")),
    )


def _axis_label(name: str) -> str:
    labels = {
        "Rg": "Rg (nm)",
        "Ree": "Ree (nm)",
        "total_contacts": "Total contacts (count)",
        "PC1": "PC1 (a.u.)",
        "PC2": "PC2 (a.u.)",
    }
    return labels.get(name, name)
