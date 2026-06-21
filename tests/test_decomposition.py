from __future__ import annotations

from types import SimpleNamespace

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

from idrptm.analysis.compare import compare_decomposition_outputs
from idrptm.analysis.decomposition import (
    contact_eigendecomposition,
    contact_map_trajectory,
    contact_pca,
    coordinate_pca,
    distance_map_trajectory,
    distance_pca,
    ev1_correlation,
    feature_pca,
    run_decomposition_analysis,
)
from idrptm.analysis.io import TrajectoryData
from idrptm.analysis.pipeline import analyze_trajectory_data
from idrptm.plotting.plots import (
    plot_contact_eigenvectors,
    plot_contact_loading_heatmap,
    plot_delta_ev,
    plot_ev1_correlation,
    plot_explained_variance,
    plot_pca_centroid_shift,
    plot_pca_score_scatter,
    plot_pca_timeseries,
)
from idrptm.schema import AnalysisConfig


def _positions() -> np.ndarray:
    return np.array(
        [
            [[0.0, 0.0, 0.0], [0.4, 0.0, 0.0], [2.0, 0.0, 0.0], [2.4, 0.0, 0.0]],
            [[0.0, 0.0, 0.0], [0.5, 0.1, 0.0], [2.2, 0.0, 0.0], [2.5, 0.1, 0.0]],
            [[0.1, 0.0, 0.0], [0.6, 0.0, 0.0], [2.1, 0.2, 0.0], [2.6, 0.2, 0.0]],
            [[0.0, 0.1, 0.0], [0.5, 0.2, 0.0], [2.3, 0.1, 0.0], [2.7, 0.2, 0.0]],
        ],
        dtype=float,
    )


def _block_contact_map() -> np.ndarray:
    matrix = np.full((6, 6), 0.05, dtype=float)
    matrix[:3, :3] = 0.9
    matrix[3:, 3:] = 0.8
    np.fill_diagonal(matrix, 0.0)
    return matrix


def test_coordinate_pca_returns_expected_shapes() -> None:
    result, representatives = coordinate_pca(
        _positions(),
        time_ps=np.array([0.0, 10.0, 20.0, 30.0]),
        n_components=2,
    )

    assert result.scores[["PC1", "PC2"]].shape == (4, 2)
    assert result.loadings.shape == (2, 12)
    assert result.explained_variance["component"].tolist() == ["PC1", "PC2"]
    assert set(representatives["representative"]) == {"minimum", "maximum", "near_zero"}


def test_contact_pca_returns_scores_loadings_and_explained_variance() -> None:
    contacts = contact_map_trajectory(_positions(), cutoff=0.8, min_sequence_separation=1)

    result, top_contacts = contact_pca(contacts, min_sequence_separation=1, n_components=2)

    assert result.scores[["PC1", "PC2"]].shape[0] == contacts.shape[0]
    assert result.loadings.shape[0] == 2
    assert not top_contacts.empty


def test_distance_and_feature_pca_handle_available_features() -> None:
    distances = distance_map_trajectory(_positions())
    distance_result, _ = distance_pca(distances, min_sequence_separation=1, n_components=2)
    feature_result = feature_pca(
        pd.DataFrame(
            {
                "frame": [0, 1, 2, 3],
                "rg": [1.0, 1.1, 1.2, 1.3],
                "ree": [2.0, 2.1, 2.2, 2.3],
            }
        ),
        n_components=2,
    )

    assert distance_result.scores[["PC1", "PC2"]].shape == (4, 2)
    assert feature_result.loadings.shape == (2, 2)


def test_contact_eigendecomposition_separates_synthetic_blocks() -> None:
    result = contact_eigendecomposition(
        _block_contact_map(),
        sequence="KKKDDD",
        min_sequence_separation=1,
        n_eigs=2,
    )

    ev1 = result.contact_eigs["EV1"].to_numpy()
    assert result.observed_expected.shape == (6, 6)
    assert result.correlation.shape == (6, 6)
    assert ev1[:3].mean() > ev1[3:].mean()


def test_contact_eigendecomposition_handles_zeros_and_masks_short_separations() -> None:
    result = contact_eigendecomposition(
        np.zeros((4, 4)),
        eps=1.0e-6,
        min_sequence_separation=2,
    )

    assert np.isnan(result.observed_expected[0, 1])
    assert np.isfinite(result.correlation).all()


def test_contact_eigendecomposition_orientation_is_deterministic() -> None:
    first = contact_eigendecomposition(_block_contact_map(), sequence="KKKDDD")
    second = contact_eigendecomposition(_block_contact_map(), sequence="KKKDDD")

    np.testing.assert_allclose(first.contact_eigs["EV1"], second.contact_eigs["EV1"])


def test_ev1_correlation_and_comparison_outputs(tmp_path) -> None:
    wt = contact_eigendecomposition(_block_contact_map(), sequence="KKKDDD").contact_eigs
    variant = wt.copy()
    variant["EV1"] = variant["EV1"] * 0.8
    assert ev1_correlation(wt, variant) == pytest.approx(1.0)

    wt_dir = tmp_path / "wt"
    ptm_dir = tmp_path / "ptm"
    wt_dir.mkdir()
    ptm_dir.mkdir()
    wt.to_csv(wt_dir / "contact_eigs.csv", index=False)
    variant.to_csv(ptm_dir / "contact_eigs.csv", index=False)
    pd.DataFrame({"PC1": [0.0, 1.0], "PC2": [0.0, 0.0]}).to_parquet(
        wt_dir / "feature_pca_scores.parquet"
    )
    pd.DataFrame({"PC1": [2.0, 3.0], "PC2": [1.0, 1.0]}).to_parquet(
        ptm_dir / "feature_pca_scores.parquet"
    )
    runs = (
        SimpleNamespace(condition="WT", analysis_dir=wt_dir),
        SimpleNamespace(
            condition="pSer",
            analysis_dir=ptm_dir,
            cleavage_state="cut_1",
            cut_number=1,
            event_time_ns=5.0,
        ),
    )

    outputs = compare_decomposition_outputs(runs, "WT", tmp_path)

    assert (tmp_path / "decomposition_comparison.csv").is_file()
    assert (tmp_path / "delta_ev.csv").is_file()
    assert (tmp_path / "pca_centroid_shift.csv").is_file()
    assert "decomposition_comparison_csv" in outputs
    comparison = pd.read_csv(tmp_path / "decomposition_comparison.csv")
    assert comparison.loc[0, "cut_number"] == 1

    figs = [
        plot_delta_ev(pd.read_csv(tmp_path / "delta_ev.csv")),
        plot_ev1_correlation(comparison),
        plot_pca_centroid_shift(pd.read_csv(tmp_path / "pca_centroid_shift.csv")),
    ]
    for fig in figs:
        assert fig.axes
        plt.close(fig)


def test_decomposition_pipeline_writes_outputs(tmp_path) -> None:
    trajectory = TrajectoryData(
        positions=_positions(),
        topology_path="top.pdb",
        trajectory_path="trajectory.dcd",
        time_ps=np.array([0.0, 10.0, 20.0, 30.0]),
    )
    config = AnalysisConfig(
        observables=["rg", "ree", "contacts", "ps", "scaling"],
        contact_cutoff_nm=0.8,
        decomposition={
            "enabled": True,
            "feature_pca": {"enabled": True, "standardize_features": True, "n_components": 2},
            "contact_eigs": {"enabled": True, "min_sequence_separation": 1, "n_eigs": 2},
        },
    )

    result = analyze_trajectory_data(
        trajectory,
        output_dir=tmp_path / "analysis",
        analysis_config=config,
    )

    assert (result.output_dir / "feature_pca_scores.parquet").is_file()
    assert (result.output_dir / "contact_eigs.csv").is_file()
    assert "feature_pca_scores" in result.outputs


def test_run_decomposition_analysis_and_plots(tmp_path) -> None:
    contact_map = _block_contact_map()
    outputs = run_decomposition_analysis(
        output_dir=tmp_path,
        positions=_positions(),
        frame_table=pd.DataFrame({"frame": [0, 1, 2, 3]}),
        contact_map=contact_map[:4, :4],
        analysis_config=SimpleNamespace(
            contact_cutoff_nm=0.8,
            decomposition={
                "enabled": True,
                "coordinate_pca": {"enabled": True, "n_components": 2},
                "contact_pca": {"enabled": True, "min_sequence_separation": 1, "n_components": 2},
                "contact_eigs": {"enabled": True, "min_sequence_separation": 1, "n_eigs": 2},
            },
        ),
    )

    assert (tmp_path / "pca_scores.parquet").is_file()
    assert (tmp_path / "top_loading_contacts.csv").is_file()
    assert "contact_eigs" in outputs

    scores = pd.read_parquet(tmp_path / "pca_scores.parquet")
    variance = pd.read_csv(tmp_path / "pca_explained_variance.csv")
    eigs = pd.read_csv(tmp_path / "contact_eigs.csv")
    figs = [
        plot_pca_score_scatter(scores),
        plot_pca_timeseries(scores),
        plot_explained_variance(variance),
        plot_contact_eigenvectors(eigs),
        plot_contact_loading_heatmap(
            np.load(tmp_path / "contact_pca_loadings.npy"),
            pd.read_csv(tmp_path / "top_loading_contacts.csv"),
        ),
    ]
    for fig in figs:
        assert fig.axes
        plt.close(fig)
