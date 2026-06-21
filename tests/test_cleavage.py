from __future__ import annotations

import csv

from idrptm.cleavage import fragments_from_cuts, generate_cleavage_states, protease_candidate_cuts
from idrptm.design import write_design_outputs
from idrptm.ptm import AppliedPTM
from idrptm.schema import CleavageSet, ProteinConfig, PTMConfig, PTMSite, WorkflowConfig


def test_trypsin_rule_cuts_after_k_or_r_not_before_p() -> None:
    cuts = protease_candidate_cuts("AKRPQRK", "trypsin_simple")

    assert cuts == (2, 6)


def test_random_cleavage_is_reproducible_with_seed() -> None:
    cleavage = CleavageSet(
        name="random_two",
        mode="random",
        n_cuts=2,
        seed=17,
        min_fragment_length=2,
    )

    first = generate_cleavage_states("AAAAAAAAAA", cleavage)
    second = generate_cleavage_states("AAAAAAAAAA", cleavage)

    assert first[0].cuts == second[0].cuts
    assert len(first[0].cuts) == 2


def test_sequential_cleavage_produces_ordered_states() -> None:
    cleavage = CleavageSet(
        name="trypsin_series",
        mode="sequential",
        protease="trypsin_simple",
        order="n_to_c",
        n_cuts=2,
    )

    states = generate_cleavage_states("AKRAK", cleavage)

    assert [state.cut_number for state in states] == [1, 2]
    assert [state.cuts for state in states] == [(2,), (2, 3)]


def test_ptm_site_mapping_survives_cleavage() -> None:
    ptm = AppliedPTM(
        biological_position=2,
        zero_based_index=1,
        ptm="pSer",
        source_residue="S",
        simulation_code="B",
    )

    products = fragments_from_cuts("ABTK", [2], ptm_sites=(ptm,))

    assert products[0].ptm_sites_1based == [2]
    assert products[1].ptm_sites_1based == []


def test_fragment_ranges_are_original_residue_numbering() -> None:
    products = fragments_from_cuts("AKRTA", [2, 4])

    assert [(product.original_start, product.original_end) for product in products] == [
        (1, 2),
        (3, 4),
        (5, 5),
    ]
    assert [product.original_indices for product in products] == [[1, 2], [3, 4], [5]]


def test_design_writes_cleavage_outputs_and_fragment_metadata(tmp_path) -> None:
    config = WorkflowConfig(
        project="cleavage_design",
        protein=ProteinConfig(
            name="fraggy",
            sequence="ASAKR",
            ptm=PTMConfig(
                mode="explicit",
                include_wt=False,
                sites=[PTMSite(position=2, residue="S", ptm="pSer")],
            ),
            cleavage_sets=[
                CleavageSet(
                    name="manual_cut",
                    mode="manual",
                    manual_cuts=[2],
                    individual_fragments=True,
                    fragment_mixture=True,
                )
            ],
        ),
    )

    result = write_design_outputs(config, output_dir=tmp_path)

    assert result.cleavage_sites_path is not None
    assert result.fragments_fasta_path is not None
    assert result.cleavage_manifest_path is not None
    assert result.cleavage_sites_path.exists()
    assert result.fragments_fasta_path.exists()
    assert result.cleavage_manifest_path.exists()

    fragments = list(csv.DictReader(result.cleavage_manifest_path.open(encoding="utf-8")))
    assert {row["original_start"] for row in fragments} >= {"1", "3"}
    assert any(row["ptm_sites_1based"] == "2" for row in fragments)
    fasta_text = result.fragments_fasta_path.read_text(encoding="utf-8")
    assert "AB" in fasta_text
