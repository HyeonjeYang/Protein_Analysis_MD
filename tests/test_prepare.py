from __future__ import annotations

import json

import yaml

from idrptm.calvados_adapter import prepare_from_config
from idrptm.schema import (
    CalvadosConfig,
    PTMConfig,
    PTMSite,
    SequenceConfig,
    SimulationConfig,
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
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def _workflow_config(residue_csv, work_dir) -> WorkflowConfig:
    return WorkflowConfig(
        project="prepare_test",
        sequence=SequenceConfig(name="prep_seq", sequence="AST"),
        ptm=PTMConfig(
            mode="explicit",
            include_wt=True,
            sites=[PTMSite(position=2, residue="S", ptm="pSer")],
        ),
        calvados=CalvadosConfig(
            model="pCALVADOS2",
            residue_parameters=residue_csv,
            ph=7.4,
            temperature_k=298.0,
            simulation=SimulationConfig(save_every_steps=100, n_frames=5),
        ),
        runner={"work_dir": work_dir},
    )


def test_prepare_creates_one_run_directory_per_manifest_row(tmp_path) -> None:
    residue_csv = tmp_path / "base_residues.csv"
    output_dir = tmp_path / "prepared"
    _write_base_residues(residue_csv)

    result = prepare_from_config(_workflow_config(residue_csv, output_dir))

    assert len(result.run_directories) == 2
    assert result.manifest_path.exists()
    for run_dir in result.run_directories:
        assert (run_dir.path / "input.fasta").exists()
        assert (run_dir.path / "residues.csv").exists()
        assert (run_dir.path / "config.yaml").exists()
        assert (run_dir.path / "components.yaml").exists()
        assert (run_dir.path / "run.py").exists()
        assert (run_dir.path / "metadata.json").exists()

    ptm_run = output_dir / "runs" / "prep_seq__pSer2"
    metadata = json.loads((ptm_run / "metadata.json").read_text(encoding="utf-8"))
    assert metadata["original_sequence"] == "AST"
    assert metadata["simulation_sequence"] == "ABT"
    assert metadata["residue_parameters"]["ph"] == 7.4
    assert metadata["residue_parameters"]["source_file"] == str(residue_csv.resolve())
    assert metadata["residue_parameters"]["ptm_charges"]["pSer"]["pka"] == 6.01
    assert metadata["simulation"]["dt_ps"] == 0.01
    assert metadata["simulation"]["save_every_steps"] == 100
    assert metadata["simulation"]["frame_interval_ns"] == 0.001
    assert metadata["simulation"]["n_frames"] == 5
    assert metadata["simulation"]["total_steps"] == 500
    assert metadata["simulation"]["total_time_ns"] == 0.005
    assert metadata["upstream_policy"]["mutates_calvados_source"] is False

    config_yaml = yaml.safe_load((ptm_run / "config.yaml").read_text(encoding="utf-8"))
    components_yaml = yaml.safe_load((ptm_run / "components.yaml").read_text(encoding="utf-8"))
    assert config_yaml["sysname"] == "prep_seq__pSer2"
    assert config_yaml["pH"] == 7.4
    assert config_yaml["wfreq"] == 100
    assert config_yaml["steps"] == 500
    assert config_yaml["runtime"] == 0
    assert config_yaml["platform"] == "CPU"
    assert config_yaml["restart"] == "checkpoint"
    assert config_yaml["frestart"] == "restart.chk"
    assert components_yaml["defaults"]["fresidues"] == "residues.csv"
    assert components_yaml["defaults"]["ffasta"] == "input.fasta"
    assert "from calvados import sim" in (ptm_run / "run.py").read_text(encoding="utf-8")


def test_prepare_dry_run_does_not_write_files(tmp_path) -> None:
    output_dir = tmp_path / "dry_run"

    config = WorkflowConfig(
        project="dry_run_test",
        sequence=SequenceConfig(name="dry_seq", sequence="AST"),
        ptm=PTMConfig(mode="wt", include_wt=True),
        runner={"work_dir": output_dir},
    )

    result = prepare_from_config(config, dry_run=True)

    assert result.dry_run is True
    assert len(result.run_directories) == 1
    assert result.run_directories[0].path == output_dir / "runs" / "dry_seq__WT"
    assert not output_dir.exists()
