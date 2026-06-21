from __future__ import annotations

from protein_analysis_md.recipe import Experiment, Protein

experiment = Experiment(name="flk_smoke", outdir="runs/flk_smoke")
experiment.add_protein(
    Protein.from_uniprot(
        query="FLK",
        reviewed_only=True,
        organism="Homo sapiens",
        interactive_select=True,
        region={"mode": "prompt_if_long", "max_length": 300},
    )
)
experiment.use_preset(
    simulation="smoke_single_chain",
    analysis="standard_idr",
    report="standard",
)
experiment.add_ptm_state("WT")
experiment.add_cleavage_state("intact")


if __name__ == "__main__":
    experiment.write_yaml("configs/flk_smoke.recipe.yaml")
    experiment.compile(force=True)
