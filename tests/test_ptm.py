from __future__ import annotations

import pytest

from idrptm.ptm import PTMRequest, apply_ptms


def test_pser_on_serine_succeeds() -> None:
    simulation_sequence, applied = apply_ptms(
        "ASG",
        (PTMRequest(biological_position=2, ptm="pSer", expected_residue="S"),),
    )

    assert simulation_sequence == "ABG"
    assert applied[0].zero_based_index == 1
    assert applied[0].simulation_code == "B"


def test_pthr_on_threonine_succeeds() -> None:
    simulation_sequence, applied = apply_ptms(
        "ATG",
        (PTMRequest(biological_position=2, ptm="pThr", expected_residue="T"),),
    )

    assert simulation_sequence == "AOG"
    assert applied[0].zero_based_index == 1
    assert applied[0].simulation_code == "O"


def test_pser_on_non_serine_fails() -> None:
    with pytest.raises(ValueError, match="pSer requires S"):
        apply_ptms(
            "ATG",
            (PTMRequest(biological_position=2, ptm="pSer"),),
        )


def test_pthr_on_non_threonine_fails() -> None:
    with pytest.raises(ValueError, match="pThr requires T"):
        apply_ptms(
            "ASG",
            (PTMRequest(biological_position=2, ptm="pThr"),),
        )
