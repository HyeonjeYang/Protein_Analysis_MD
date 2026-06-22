"""Publication-draft visualization helpers for protein_analysis_md."""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from matplotlib.figure import Figure

from idrptm.plotting.plots import save_figure
from idrptm.visualization.pymol import (
    PyMOLExportResult,
    PyMOLRunExport,
    export_pymol_project,
)


@dataclass(frozen=True)
class VisualizationArtifact:
    """Paths written for one visualization."""

    png: Path
    pdf: Path | None
    data: Path
    metadata: Path | None = None


def save_visualization(
    fig: Figure,
    data: pd.DataFrame | np.ndarray | dict[str, object],
    output_base: str | Path,
    *,
    metadata: dict[str, object] | None = None,
) -> VisualizationArtifact:
    """Save figure plus underlying plotting data without overwriting raw analysis data."""

    figure_paths = save_figure(fig, output_base)
    png = next(path for path in figure_paths if path.suffix == ".png")
    pdf = next((path for path in figure_paths if path.suffix == ".pdf"), None)
    base = Path(output_base)
    if isinstance(data, pd.DataFrame):
        data_path = base.with_suffix(".csv")
        data.to_csv(data_path, index=False)
    elif isinstance(data, np.ndarray):
        data_path = base.with_suffix(".npy")
        np.save(data_path, data)
    else:
        data_path = base.with_suffix(".json")
        data_path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    metadata_path = None
    if metadata is not None:
        metadata_path = base.with_suffix(".metadata.json")
        metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return VisualizationArtifact(png=png, pdf=pdf, data=data_path, metadata=metadata_path)


__all__ = [
    "PyMOLExportResult",
    "PyMOLRunExport",
    "VisualizationArtifact",
    "export_pymol_project",
    "save_visualization",
]
