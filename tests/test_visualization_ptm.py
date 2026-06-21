from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pytest

from idrptm.visualization.ptm import (
    delta_ptm_site_profile,
    plot_ptm_delta_contact_map,
    plot_ptm_site_profile,
    plot_residue_class_contact_changes,
    residue_class_contact_changes,
)


def _maps() -> tuple[np.ndarray, np.ndarray]:
    wt = np.array([[0.0, 0.2, 0.1], [0.2, 0.0, 0.3], [0.1, 0.3, 0.0]])
    ptm = wt + np.array([[0.0, 0.1, -0.05], [0.1, 0.0, 0.2], [-0.05, 0.2, 0.0]])
    return wt, ptm


def test_ptm_delta_contact_map_uses_raw_maps() -> None:
    wt, ptm = _maps()
    fig = plot_ptm_delta_contact_map(wt, ptm)
    profile = delta_ptm_site_profile(wt, ptm, [2])

    assert profile.loc[0, "delta_contact_probability"] == pytest.approx(0.1)
    assert fig.axes
    plt.close(fig)


def test_ptm_profile_and_residue_class_changes_plot() -> None:
    wt, ptm = _maps()
    changes = residue_class_contact_changes(wt, ptm, "KDF")
    figs = [
        plot_ptm_site_profile(ptm, [2]),
        plot_residue_class_contact_changes(changes),
    ]

    for fig in figs:
        assert fig.axes
        plt.close(fig)
