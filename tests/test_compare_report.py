from __future__ import annotations

import csv
import json

import numpy as np
import pandas as pd
import pytest

from idrptm.analysis.compare import compare_project, detect_wt_condition
from idrptm.plotting.report import generate_report


def _write_analysis(
    analysis_dir,
    *,
    rg_values: list[float],
    ree_values: list[float],
    flory: float,
    contact_map: np.ndarray,
    ps_values: dict[int, float],
    scaling_values: dict[int, float],
) -> None:
    analysis_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame({"frame": range(len(rg_values)), "rg": rg_values}).to_parquet(
        analysis_dir / "timeseries_rg.parquet",
        index=False,
    )
    pd.DataFrame({"frame": range(len(ree_values)), "ree": ree_values}).to_parquet(
        analysis_dir / "timeseries_ree.parquet",
        index=False,
    )
    np.save(analysis_dir / "contact_map.npy", contact_map)
    pd.DataFrame(
        [{"s": separation, "p": value, "n_pairs": 2} for separation, value in ps_values.items()]
    ).to_parquet(analysis_dir / "ps.parquet", index=False)
    pd.DataFrame(
        [
            {"s": separation, "distance": value, "n_pairs": 2}
            for separation, value in scaling_values.items()
        ]
    ).to_parquet(analysis_dir / "scaling.parquet", index=False)
    (analysis_dir / "summary.json").write_text(
        json.dumps({"flory_fit": {"nu": flory}}) + "\n",
        encoding="utf-8",
    )


def _write_project(project_dir) -> None:
    rows = [
        {
            "variant_id": "toy__WT",
            "condition": "WT",
            "replicate": "WT_1",
            "ptm_state": "WT",
            "ptm_sites_1based": "",
            "original_sequence": "AST",
            "metadata_path": "runs/toy__WT/metadata.yaml",
        },
        {
            "variant_id": "toy__pSer2_rep1",
            "condition": "pSer2",
            "replicate": "rep1",
            "ptm_state": "pSer2",
            "ptm_sites_1based": "pSer:S2->B",
            "original_sequence": "AST",
            "metadata_path": "runs/toy__pSer2_rep1/metadata.yaml",
        },
        {
            "variant_id": "toy__pSer2_rep2",
            "condition": "pSer2",
            "replicate": "rep2",
            "ptm_state": "pSer2",
            "ptm_sites_1based": "pSer:S2->B",
            "original_sequence": "AST",
            "metadata_path": "runs/toy__pSer2_rep2/metadata.yaml",
        },
    ]
    manifest = project_dir / "manifest.csv"
    with manifest.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    wt_map = np.array(
        [
            [0.0, 0.2, 0.1],
            [0.2, 0.0, 0.3],
            [0.1, 0.3, 0.0],
        ]
    )
    _write_analysis(
        project_dir / "runs" / "toy__WT" / "analysis",
        rg_values=[1.0, 1.2],
        ree_values=[3.0, 3.2],
        flory=0.5,
        contact_map=wt_map,
        ps_values={1: 0.2, 2: 0.1},
        scaling_values={1: 1.0, 2: 1.4},
    )
    _write_analysis(
        project_dir / "runs" / "toy__pSer2_rep1" / "analysis",
        rg_values=[1.4, 1.6],
        ree_values=[3.4, 3.6],
        flory=0.6,
        contact_map=wt_map + 0.1,
        ps_values={1: 0.4, 2: 0.2},
        scaling_values={1: 1.2, 2: 1.7},
    )
    _write_analysis(
        project_dir / "runs" / "toy__pSer2_rep2" / "analysis",
        rg_values=[1.2, 1.4],
        ree_values=[3.6, 3.8],
        flory=0.8,
        contact_map=wt_map + 0.2,
        ps_values={1: 0.6, 2: 0.4},
        scaling_values={1: 1.3, 2: 1.9},
    )


def test_detect_wt_condition_by_name() -> None:
    assert detect_wt_condition(["pSer2", "toy__WT"]) == "toy__WT"


def test_compare_project_aggregates_replicates_and_deltas(tmp_path) -> None:
    _write_project(tmp_path)

    comparison = compare_project(tmp_path)

    assert comparison.wt_condition == "WT"
    pser = comparison.summary_table.set_index("condition").loc["pSer2"]
    assert pser["n_replicates"] == 2
    assert pser["mean_Rg"] == pytest.approx(1.4)
    assert pser["std_Rg"] == pytest.approx(np.sqrt(0.02))
    assert pser["sem_Rg"] == pytest.approx(0.1)
    assert pser["delta_mean_Rg"] == pytest.approx(0.3)
    assert pser["delta_mean_Ree"] == pytest.approx(0.5)
    assert pser["delta_flory_exponent"] == pytest.approx(0.2)

    delta_ps = comparison.ps_delta.set_index("s")
    assert delta_ps.loc[1, "delta_p"] == pytest.approx(0.3)
    assert delta_ps.loc[2, "delta_p"] == pytest.approx(0.2)
    assert comparison.delta_maps["pSer2"][0, 1] == pytest.approx(0.15)
    assert (tmp_path / "comparison" / "comparison_summary.csv").exists()
    assert (tmp_path / "comparison" / "delta_contact_map_pSer2_minus_WT.npy").exists()


def test_generate_report_writes_markdown_and_png_figures_by_default(tmp_path) -> None:
    _write_project(tmp_path)

    report = generate_report(tmp_path)

    assert report.report_path.exists()
    text = report.report_path.read_text(encoding="utf-8")
    assert "WT condition" in text
    assert "Delta Contact Map Pser2" in text
    figure_names = {path.name for path in report.figure_paths}
    assert "rg_distribution.png" in figure_names
    assert "delta_contact_map_pSer2.png" in figure_names
    assert "summary_table.png" in figure_names
    assert not any(name.endswith(".pdf") for name in figure_names)


def test_generate_report_can_write_pdf_when_requested(tmp_path, monkeypatch) -> None:
    monkeypatch.setenv("PAMD_FIGURE_FORMATS", "png,pdf")
    _write_project(tmp_path)

    report = generate_report(tmp_path)

    figure_names = {path.name for path in report.figure_paths}
    assert "rg_distribution.png" in figure_names
    assert "rg_distribution.pdf" in figure_names
    assert "delta_contact_map_pSer2.pdf" in figure_names
    assert "summary_table.pdf" in figure_names
