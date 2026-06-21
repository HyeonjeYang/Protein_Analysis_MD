from __future__ import annotations

import numpy as np

from idrptm.analysis.free_energy import free_energy_surface_2d, write_free_energy_surface
from idrptm.visualization.free_energy import plot_free_energy_surface


def test_free_energy_rg_ree_grid_and_raw_counts(tmp_path) -> None:
    surface = free_energy_surface_2d(
        np.array([1.0, 1.1, 1.2, 1.3]),
        np.array([2.0, 2.1, 2.1, 2.3]),
        bins=3,
        temperature_k=293.0,
    )
    paths = write_free_energy_surface(surface, tmp_path, "free_energy_Rg_Ree")
    fig = plot_free_energy_surface(surface)

    assert surface.counts.shape == (3, 3)
    assert paths["counts"].is_file()
    assert paths["free_energy"].is_file()
    assert fig.axes
