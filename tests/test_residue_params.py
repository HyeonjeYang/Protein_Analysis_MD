from __future__ import annotations

import csv

import pytest

from idrptm.residue_params import (
    phosphorylated_charge,
    resolve_residue_source,
    write_residue_parameters,
)


def _write_base_residues(path) -> None:
    path.write_text(
        "\n".join(
            [
                "one,three,MW,lambdas,sigmas,q,bondlength",
                "A,ALA,71.07,0.274329796904,0.504,0.0,0.38",
                "S,SER,87.08,0.462541681161,0.518,0.0,0.38",
                "T,THR,101.11,0.371316297627,0.562,0.0,0.38",
                "B,SEP,165.04,0.0,0.0,0.0,0.38",
            ]
        )
        + "\n",
        encoding="utf-8",
    )


def test_write_residue_parameters_adds_and_overrides_phosphorylated_rows(tmp_path) -> None:
    source = tmp_path / "base_residues.csv"
    output = tmp_path / "run" / "residues.csv"
    _write_base_residues(source)

    result = write_residue_parameters(source_csv=source, output_csv=output, ph=7.4)

    rows = {row["one"]: row for row in csv.DictReader(output.open(encoding="utf-8"))}
    assert rows["B"]["three"] == "SEP"
    assert rows["O"]["three"] == "TPO"
    assert float(rows["B"]["q"]) == pytest.approx(phosphorylated_charge("pSer", 7.4))
    assert float(rows["O"]["q"]) == pytest.approx(phosphorylated_charge("pThr", 7.4))
    assert result.metadata()["source_file"] == str(source.resolve())
    assert result.metadata()["ptm_charges"]["pSer"]["pka"] == 6.01
    assert result.metadata()["ptm_charges"]["pThr"]["pka"] == 6.30


def test_residue_source_can_come_from_environment(tmp_path, monkeypatch) -> None:
    source = tmp_path / "base_residues.csv"
    _write_base_residues(source)
    monkeypatch.setenv("IDRPTM_CALVADOS_RESIDUES", str(source))

    assert resolve_residue_source() == source.resolve()
