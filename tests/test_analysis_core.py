from __future__ import annotations

import numpy as np
import pytest

from idrptm.analysis import (
    com_msd,
    contact_lifetime,
    contact_map_from_positions,
    fit_flory_exponent,
    internal_distance_scaling,
    p_of_s,
    ree_timeseries,
    rg_timeseries,
)


def test_rg_timeseries_uniform_and_weighted() -> None:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[10.0, 0.0, 0.0], [12.0, 0.0, 0.0]],
        ]
    )

    assert np.allclose(rg_timeseries(positions), [1.0, 1.0])

    weighted = np.array([[[0.0, 0.0, 0.0], [4.0, 0.0, 0.0]]])
    assert rg_timeseries(weighted, masses=[1.0, 3.0])[0] == pytest.approx(np.sqrt(3.0))


def test_ree_timeseries_uses_terminal_distance() -> None:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 4.0, 0.0]],
            [[1.0, 1.0, 0.0], [2.0, 1.0, 0.0], [4.0, 5.0, 0.0]],
        ]
    )

    assert np.allclose(ree_timeseries(positions), [5.0, 5.0])


def test_contact_map_from_positions_returns_contact_frequencies() -> None:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [0.5, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0], [0.75, 0.0, 0.0]],
        ]
    )

    contact_map = contact_map_from_positions(positions, cutoff=1.0)

    assert contact_map[0, 0] == 0.0
    assert contact_map[0, 1] == pytest.approx(0.5)
    assert contact_map[0, 2] == pytest.approx(0.5)
    assert contact_map[1, 2] == 0.0
    assert np.allclose(contact_map, contact_map.T)


def test_p_of_s_averages_contact_map_diagonals() -> None:
    contact_map = np.array(
        [
            [0.0, 1.0, 0.2],
            [1.0, 0.0, 0.6],
            [0.2, 0.6, 0.0],
        ]
    )

    ps = p_of_s(contact_map)

    assert ps["s"].tolist() == [1, 2]
    assert ps["n_pairs"].tolist() == [2, 1]
    assert ps.loc[ps["s"] == 1, "p"].item() == pytest.approx(0.8)
    assert ps.loc[ps["s"] == 2, "p"].item() == pytest.approx(0.2)


def test_internal_distance_scaling_on_linear_chain() -> None:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
        ]
    )

    scaling = internal_distance_scaling(positions)

    assert scaling["s"].tolist() == [1, 2]
    assert scaling["distance"].tolist() == pytest.approx([1.5, 3.0])
    assert scaling["mean_square_distance_nm2"].tolist() == pytest.approx([2.5, 9.0])
    assert scaling["rms_distance_nm"].tolist() == pytest.approx([np.sqrt(2.5), 3.0])
    assert scaling["n_pairs"].tolist() == [2, 1]


def test_fit_flory_exponent_recovers_power_law() -> None:
    s = np.array([1.0, 2.0, 4.0, 8.0])
    distances = 2.0 * np.sqrt(s)

    fit = fit_flory_exponent(s=s, distances=distances)

    assert fit.nu == pytest.approx(0.5)
    assert fit.prefactor == pytest.approx(2.0)
    assert fit.n_points == 4


def test_contact_lifetime_returns_normalized_intermittent_correlation() -> None:
    contacts = np.array([1, 1, 0, 1])

    lifetime = contact_lifetime(contacts, max_lag=2)

    assert lifetime["lag"].tolist() == [0, 1, 2]
    assert lifetime.loc[lifetime["lag"] == 0, "correlation"].item() == pytest.approx(1.0)
    assert lifetime.loc[lifetime["lag"] == 1, "correlation"].item() == pytest.approx(4.0 / 9.0)
    assert lifetime.loc[lifetime["lag"] == 2, "correlation"].item() == pytest.approx(2.0 / 3.0)


def test_com_msd_tracks_center_of_mass_displacement() -> None:
    positions = np.array(
        [
            [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
            [[1.0, 0.0, 0.0], [3.0, 0.0, 0.0]],
            [[3.0, 0.0, 0.0], [5.0, 0.0, 0.0]],
        ]
    )

    msd = com_msd(positions, max_lag=2)

    assert msd["lag"].tolist() == [0, 1, 2]
    assert msd["msd"].tolist() == pytest.approx([0.0, 2.5, 9.0])
    assert msd["n_origins"].tolist() == [3, 2, 1]
