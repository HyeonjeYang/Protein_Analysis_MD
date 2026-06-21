from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from idrptm.visualization.phase import (
    cluster_distribution_figure,
    phase_reliability_warning,
    plot_density_profile,
    plot_density_projection,
    plot_inter_chain_contact_heatmap,
    plot_partition_coefficients,
)


def test_density_profile_and_projection_plots() -> None:
    density = pd.DataFrame(
        {
            "z_nm": [0.0, 1.0, 2.0],
            "density": [0.1, 0.5, 0.2],
            "component": ["A", "A", "A"],
        }
    )
    positions = np.array([[0.0, 0.0], [1.0, 1.0], [1.5, 0.5]])
    figs = [plot_density_profile(density), plot_density_projection(positions, bins=3)]

    for fig in figs:
        assert fig.axes
        plt.close(fig)


def test_phase_cluster_partition_and_contact_plots() -> None:
    cluster = pd.DataFrame({"largest_cluster_size": [1, 2, 2, 3]})
    partition = pd.DataFrame(
        {
            "component": ["A", "B"],
            "dense_concentration": [2.0, 1.0],
            "dilute_concentration": [1.0, 2.0],
        }
    )
    figs = [
        cluster_distribution_figure(cluster),
        plot_partition_coefficients(partition),
        plot_inter_chain_contact_heatmap(np.array([[0.0, 0.4], [0.4, 0.0]])),
    ]

    assert phase_reliability_warning(10, 2) is not None
    for fig in figs:
        assert fig.axes
        plt.close(fig)
