"""Free-energy landscape helpers for exploratory visualization."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from numpy.typing import ArrayLike
from scipy import ndimage

KB_KJ_PER_MOL_K = 0.00831446261815324


@dataclass(frozen=True)
class FreeEnergySurface:
    """Binned two-dimensional free-energy surface."""

    x_edges: np.ndarray
    y_edges: np.ndarray
    counts: np.ndarray
    probability: np.ndarray
    free_energy_kj_mol: np.ndarray
    metadata: dict[str, object]


def free_energy_surface_2d(
    x: ArrayLike,
    y: ArrayLike,
    *,
    bins: int | tuple[int, int] = 50,
    temperature_k: float = 293.0,
    min_count: int = 1,
    display_smoothing_sigma: float | None = None,
    x_label: str = "x",
    y_label: str = "y",
) -> FreeEnergySurface:
    """Compute ``F(x, y) = -kBT ln P(x, y)`` from raw 2D histogram counts."""

    x_values = np.asarray(x, dtype=float)
    y_values = np.asarray(y, dtype=float)
    if x_values.shape != y_values.shape or x_values.ndim != 1:
        raise ValueError("x and y must be one-dimensional arrays with the same shape.")
    if temperature_k <= 0:
        raise ValueError("temperature_k must be positive.")
    if min_count <= 0:
        raise ValueError("min_count must be positive.")

    finite = np.isfinite(x_values) & np.isfinite(y_values)
    if not finite.any():
        raise ValueError("free_energy_surface_2d requires at least one finite sample.")
    counts, x_edges, y_edges = np.histogram2d(x_values[finite], y_values[finite], bins=bins)
    display_counts = counts.copy()
    if display_smoothing_sigma is not None:
        if display_smoothing_sigma <= 0:
            raise ValueError("display_smoothing_sigma must be positive.")
        display_counts = ndimage.gaussian_filter(display_counts, sigma=display_smoothing_sigma)
    probability = display_counts / display_counts.sum() if display_counts.sum() else display_counts
    with np.errstate(divide="ignore", invalid="ignore"):
        free_energy = -KB_KJ_PER_MOL_K * temperature_k * np.log(probability)
    free_energy[counts < min_count] = np.nan
    finite_energy = free_energy[np.isfinite(free_energy)]
    if finite_energy.size:
        free_energy = free_energy - float(finite_energy.min())
    metadata = {
        "x_label": x_label,
        "y_label": y_label,
        "temperature_K": temperature_k,
        "min_count": min_count,
        "raw_counts_preserved": True,
        "display_smoothing_sigma": display_smoothing_sigma,
        "exploratory": True,
        "warning": "Sampling-dependent surface; do not infer kinetic barriers from short runs.",
    }
    return FreeEnergySurface(
        x_edges=x_edges,
        y_edges=y_edges,
        counts=counts,
        probability=probability,
        free_energy_kj_mol=free_energy,
        metadata=metadata,
    )


def write_free_energy_surface(
    surface: FreeEnergySurface,
    output_dir: str | Path,
    prefix: str,
) -> dict[str, Path]:
    """Write raw counts, free energy grid, bin edges, and metadata."""

    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    paths = {
        "counts": root / f"{prefix}_counts.npy",
        "free_energy": root / f"{prefix}_free_energy.npy",
        "x_edges": root / f"{prefix}_x_edges.npy",
        "y_edges": root / f"{prefix}_y_edges.npy",
        "metadata": root / f"{prefix}_metadata.json",
        "grid_csv": root / f"{prefix}_grid.csv",
    }
    np.save(paths["counts"], surface.counts)
    np.save(paths["free_energy"], surface.free_energy_kj_mol)
    np.save(paths["x_edges"], surface.x_edges)
    np.save(paths["y_edges"], surface.y_edges)
    paths["metadata"].write_text(json.dumps(surface.metadata, indent=2) + "\n", encoding="utf-8")
    _surface_table(surface).to_csv(paths["grid_csv"], index=False)
    return paths


def run_free_energy_analysis(
    *,
    output_dir: str | Path,
    frame_table: pd.DataFrame,
    analysis_config: object,
    rg: ArrayLike | None = None,
    ree: ArrayLike | None = None,
) -> dict[str, Path]:
    """Run configured free-energy landscapes from available per-frame metrics."""

    config = getattr(analysis_config, "free_energy", {}) or {}
    if not isinstance(config, dict) or not config.get("enabled", False):
        return {}
    features = frame_table.copy()
    if rg is not None:
        features["Rg"] = np.asarray(rg, dtype=float)
    if ree is not None:
        features["Ree"] = np.asarray(ree, dtype=float)

    variables = config.get("variables", [["Rg", "Ree"]])
    if not isinstance(variables, list):
        raise ValueError("analysis.free_energy.variables must be a list of two-column lists.")
    outputs: dict[str, Path] = {}
    for variable_pair in variables:
        if not isinstance(variable_pair, list | tuple) or len(variable_pair) != 2:
            raise ValueError("Each free-energy variable pair must contain exactly two names.")
        x_name, y_name = str(variable_pair[0]), str(variable_pair[1])
        if x_name not in features or y_name not in features:
            continue
        prefix = f"free_energy_{_slug(x_name)}_{_slug(y_name)}"
        surface = free_energy_surface_2d(
            features[x_name].to_numpy(dtype=float),
            features[y_name].to_numpy(dtype=float),
            bins=int(config.get("bins", 50)),
            temperature_k=float(config.get("temperature_K", 293.0)),
            min_count=int(config.get("min_count", 1)),
            display_smoothing_sigma=(
                float(config["display_smoothing_sigma"])
                if config.get("display_smoothing_sigma") is not None
                else None
            ),
            x_label=x_name,
            y_label=y_name,
        )
        for name, path in write_free_energy_surface(surface, output_dir, prefix).items():
            outputs[f"{prefix}_{name}"] = path
    return outputs


def _surface_table(surface: FreeEnergySurface) -> pd.DataFrame:
    x_centers = (surface.x_edges[:-1] + surface.x_edges[1:]) / 2.0
    y_centers = (surface.y_edges[:-1] + surface.y_edges[1:]) / 2.0
    rows: list[dict[str, float]] = []
    for i, x_value in enumerate(x_centers):
        for j, y_value in enumerate(y_centers):
            rows.append(
                {
                    "x": float(x_value),
                    "y": float(y_value),
                    "count": float(surface.counts[i, j]),
                    "probability": float(surface.probability[i, j]),
                    "free_energy_kj_mol": float(surface.free_energy_kj_mol[i, j]),
                }
            )
    return pd.DataFrame(rows)


def _slug(value: str) -> str:
    slug = "".join(character.lower() if character.isalnum() else "_" for character in value)
    return slug.strip("_")
