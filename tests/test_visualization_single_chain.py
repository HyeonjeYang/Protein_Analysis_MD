from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idrptm.visualization import save_visualization
from idrptm.visualization.single_chain import (
    local_scaling_exponent,
    plot_contact_degree,
    plot_local_scaling_exponent,
    plot_observed_expected_contact_map,
    plot_ptm_centered_contact_profile,
    plot_rg_ree_hexbin,
    plot_rg_ree_timeseries,
    rg_ree_plotting_data,
)


def _rg() -> pd.DataFrame:
    return pd.DataFrame(
        {"frame": [0, 1, 2], "time_ps": [0.0, 10.0, 20.0], "rg": [1.0, 1.2, 1.1]}
    )


def _ree() -> pd.DataFrame:
    return pd.DataFrame(
        {"frame": [0, 1, 2], "time_ps": [0.0, 10.0, 20.0], "ree": [2.4, 2.6, 2.5]}
    )


def _contact_map() -> np.ndarray:
    return np.array([[0.0, 0.2, 0.1], [0.2, 0.0, 0.4], [0.1, 0.4, 0.0]])


def test_rg_ree_plots_and_artifact_save(tmp_path) -> None:
    table = rg_ree_plotting_data(_rg(), _ree())
    fig = plot_rg_ree_hexbin(_rg(), _ree())
    artifact = save_visualization(fig, table, tmp_path / "rg_ree")

    assert artifact.png.is_file()
    assert artifact.pdf.is_file()
    assert artifact.data.is_file()

    fig = plot_rg_ree_timeseries(_rg(), _ree())
    assert fig.axes
    plt.close(fig)


def test_ps_rs_contact_and_ptm_profile_plots() -> None:
    scaling = pd.DataFrame({"s": [1, 2, 3], "mean_distance_nm": [0.5, 0.8, 1.1]})
    local = local_scaling_exponent(scaling)
    figs = [
        plot_local_scaling_exponent(local),
        plot_contact_degree(_contact_map()),
        plot_observed_expected_contact_map(_contact_map()),
        plot_ptm_centered_contact_profile(_contact_map(), [2]),
    ]

    for fig in figs:
        assert fig.axes
        plt.close(fig)
