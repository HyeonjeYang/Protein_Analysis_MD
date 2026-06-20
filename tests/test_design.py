from __future__ import annotations

import csv

from idrptm.design import write_design_outputs
from idrptm.schema import (
    CalvadosConfig,
    PTMConfig,
    PTMSite,
    RunnerConfig,
    SequenceConfig,
    WorkflowConfig,
)


def _workflow_config() -> WorkflowConfig:
    return WorkflowConfig(
        project="deterministic_scan",
        sequence=SequenceConfig(name="toy_idr", sequence="MSTG"),
        ptm=PTMConfig(
            mode="single_site_scan",
            include_wt=True,
            sites=[
                PTMSite(position=3, residue="T", ptm="pThr"),
                PTMSite(position=2, residue="S", ptm="pSer"),
            ],
        ),
        calvados=CalvadosConfig(model="pCALVADOS2"),
        runner=RunnerConfig(work_dir="runs"),
    )


def test_manifest_is_deterministic(tmp_path) -> None:
    first = write_design_outputs(_workflow_config(), output_dir=tmp_path / "first")
    second = write_design_outputs(_workflow_config(), output_dir=tmp_path / "second")

    first_manifest = first.manifest_path.read_text(encoding="utf-8")
    second_manifest = second.manifest_path.read_text(encoding="utf-8")

    assert first_manifest == second_manifest

    rows = list(csv.DictReader(first.manifest_path.open(encoding="utf-8")))
    assert [row["variant_id"] for row in rows] == [
        "toy_idr__WT",
        "toy_idr__pSer2",
        "toy_idr__pThr3",
    ]
    assert [row["simulation_sequence"] for row in rows] == ["MSTG", "MBTG", "MSOG"]
    assert rows[1]["ptm_sites_0based"] == "pSer:1:S->B"
    assert rows[2]["ptm_sites_0based"] == "pThr:2:T->O"
