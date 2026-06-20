from __future__ import annotations

import csv
import json

import numpy as np
import pytest
import yaml

from idrptm.analysis.multichain import (
    chain_com_msd,
    cluster_size_timeseries,
    com_distance_timeseries,
    inter_chain_contact_map,
    inter_protein_contact_lifetime,
    intra_chain_contact_map,
    per_chain_ree,
    per_chain_rg,
)
from idrptm.calvados_adapter import prepare_from_config
from idrptm.design import write_design_outputs
from idrptm.schema import (
    CalvadosConfig,
    ProteinConfig,
    PTMConfig,
    PTMSite,
    SystemComponentConfig,
    SystemSetConfig,
    WorkflowConfig,
)


def _write_base_residues(path) -> None:
    path.write_text(
        "\n".join(
            [
                "one,three,MW,lambdas,sigmas,q,bondlength",
                "A,ALA,71.07,0.274329796904,0.504,0.0,0.38",
                "S,SER,87.08,0.462541681161,0.518,0.0,0.38",
                "T,THR,101.11,0.371316297627,0.562,0.0,0.38",
                "G,GLY,57.05,0.705884373326,0.45,0.0,0.38",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _multi_config(residue_csv=None, work_dir=None) -> WorkflowConfig:
    return WorkflowConfig(
        project="mixed_system",
        proteins=[
            ProteinConfig(
                name="alpha",
                sequence="ASG",
                ptm=PTMConfig(
                    mode="single_site_scan",
                    include_wt=True,
                    sites=[PTMSite(position=2, residue="S", ptm="pSer")],
                ),
            ),
            ProteinConfig(
                name="beta",
                sequence="ATG",
                ptm=PTMConfig(
                    mode="explicit",
                    include_wt=True,
                    sites=[PTMSite(position=2, residue="T", ptm="pThr")],
                ),
            ),
        ],
        system_sets=[
            SystemSetConfig(
                name="alpha_beta_mixed",
                placement="random",
                components=[
                    SystemComponentConfig(
                        protein="alpha",
                        ptm_state="WT",
                        component_name="alpha_wt",
                    ),
                    SystemComponentConfig(
                        protein="alpha",
                        ptm_state="pSer2",
                        component_name="alpha_pser",
                    ),
                    SystemComponentConfig(
                        protein="beta",
                        ptm_state="pThr2",
                        copies=2,
                        component_name="beta_pthr",
                    ),
                ],
            )
        ],
        calvados=(
            CalvadosConfig(residue_parameters=residue_csv)
            if residue_csv
            else CalvadosConfig()
        ),
        runner={"work_dir": work_dir or "runs"},
    )


def _two_chain_positions() -> tuple[np.ndarray, np.ndarray]:
    positions = np.array(
        [
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 0.5, 0.0],
                [1.0, 0.5, 0.0],
            ],
            [
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 5.0, 0.0],
                [1.0, 5.0, 0.0],
            ],
        ],
        dtype=float,
    )
    chain_ids = np.array(["A", "A", "B", "B"])
    return positions, chain_ids


def test_multi_protein_design_writes_mixed_system_manifest(tmp_path) -> None:
    result = write_design_outputs(_multi_config(work_dir=tmp_path), output_dir=tmp_path)

    rows = list(csv.DictReader(result.manifest_path.open(encoding="utf-8")))
    assert len(rows) == 1
    row = rows[0]
    components = json.loads(row["components_json"])

    assert row["system_name"] == "alpha_beta_mixed"
    assert row["is_multi_component"] == "1"
    assert row["placement"] == "random"
    assert row["total_molecule_copies"] == "4"
    assert [component["component_name"] for component in components] == [
        "alpha_wt",
        "alpha_pser",
        "beta_pthr",
    ]
    assert [component["simulation_sequence"] for component in components] == [
        "ASG",
        "ABG",
        "AOG",
    ]
    fasta_text = result.fasta_paths[0].read_text(encoding="utf-8")
    assert ">alpha_wt" in fasta_text
    assert ">alpha_pser" in fasta_text
    assert ">beta_pthr" in fasta_text


def test_single_protein_config_normalizes_to_proteins(tmp_path) -> None:
    config = WorkflowConfig(
        project="single_protein",
        protein=ProteinConfig(name="solo", sequence="AST"),
    )

    result = write_design_outputs(config, output_dir=tmp_path)
    rows = list(csv.DictReader(result.manifest_path.open(encoding="utf-8")))

    assert len(config.proteins) == 1
    assert rows[0]["variant_id"] == "solo__WT"
    assert rows[0]["component_count"] == "1"


def test_prepare_writes_multi_component_calvados_inputs(tmp_path) -> None:
    residue_csv = tmp_path / "base_residues.csv"
    _write_base_residues(residue_csv)

    result = prepare_from_config(
        _multi_config(residue_csv=residue_csv, work_dir=tmp_path / "prepared")
    )

    assert len(result.run_directories) == 1
    run_dir = result.run_directories[0].path
    config_yaml = yaml.safe_load((run_dir / "config.yaml").read_text(encoding="utf-8"))
    components_yaml = yaml.safe_load((run_dir / "components.yaml").read_text(encoding="utf-8"))
    metadata = json.loads((run_dir / "metadata.json").read_text(encoding="utf-8"))

    assert config_yaml["topol"] == "random"
    assert set(components_yaml["system"]) == {"alpha_wt", "alpha_pser", "beta_pthr"}
    assert components_yaml["system"]["beta_pthr"]["nmol"] == 2
    assert components_yaml["system"]["alpha_pser"]["molecule_type"] == "protein"
    assert metadata["is_multi_component"] is True
    assert len(metadata["components"]) == 3


def test_synthetic_two_chain_analysis_functions() -> None:
    positions, chain_ids = _two_chain_positions()

    rg = per_chain_rg(positions, chain_ids)
    ree = per_chain_ree(positions, chain_ids)
    intra = intra_chain_contact_map(positions, chain_ids, cutoff=1.1)
    inter = inter_chain_contact_map(positions, chain_ids, cutoff=0.6)
    com_distance = com_distance_timeseries(positions, chain_ids)
    msd = chain_com_msd(positions, chain_ids, max_lag=1)
    clusters = cluster_size_timeseries(positions, chain_ids, cutoff=0.6)
    lifetime = inter_protein_contact_lifetime(positions, chain_ids, cutoff=0.6, max_lag=1)

    assert rg.groupby("chain_id")["rg"].mean().to_dict() == pytest.approx({"A": 0.5, "B": 0.5})
    assert ree.groupby("chain_id")["ree"].mean().to_dict() == pytest.approx({"A": 1.0, "B": 1.0})
    assert set(intra) == {"A", "B"}
    assert intra["A"][0, 1] == pytest.approx(1.0)
    assert inter[0, 2] == pytest.approx(0.5)
    assert com_distance["distance"].tolist() == pytest.approx([0.5, 5.0])
    assert msd[(msd["chain_id"] == "A") & (msd["lag"] == 1)]["msd"].item() == pytest.approx(0.0)
    assert msd[(msd["chain_id"] == "B") & (msd["lag"] == 1)]["msd"].item() == pytest.approx(20.25)
    assert clusters["largest_cluster_size"].tolist() == [2, 1]
    assert lifetime["lag"].tolist() == [0, 1]
    assert lifetime.loc[0, "correlation"] == pytest.approx(1.0)
    assert lifetime.loc[1, "raw_probability"] == pytest.approx(0.0)
