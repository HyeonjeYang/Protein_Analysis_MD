from __future__ import annotations

from protein_analysis_md.recipe import Experiment, Protein


def test_recipe_compiles_direct_sequence(tmp_path) -> None:
    exp = Experiment(name="recipe", outdir=tmp_path / "runs" / "recipe")
    exp.add_protein(Protein.from_sequence("seq", "AST"))
    exp.use_preset(simulation="smoke_single_chain", analysis="minimal", report="minimal")
    exp.add_ptm_state("WT")
    exp.add_cleavage_state("intact")

    locked = exp.compile()

    assert locked.lock_yaml.exists()
    assert locked.workflow.project == "recipe"
    assert locked.workflow.proteins[0].name == "seq"


def test_recipe_writes_yaml(tmp_path) -> None:
    exp = Experiment(name="yaml_recipe", outdir=tmp_path / "runs" / "yaml_recipe")
    exp.add_protein(Protein.from_fasta("seq", "seq.fasta"))

    output = exp.write_yaml(tmp_path / "recipe.yaml")

    assert output.exists()
    assert "yaml_recipe" in output.read_text(encoding="utf-8")
