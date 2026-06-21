# Python Recipes

Recipes use `protein_analysis_md.recipe.Experiment` and `Protein` to build the
same concise config structure available in YAML.

```python
from protein_analysis_md.recipe import Experiment, Protein

exp = Experiment(name="example", outdir="runs/example")
exp.add_protein(Protein.from_sequence("seq", "MSSSSPST"))
exp.use_preset(simulation="smoke_single_chain", analysis="standard_idr")
exp.add_ptm_state("WT")
exp.add_cleavage_state("intact")
exp.compile()
```

Run recipes directly with Python or through `pamd run-recipe recipes/file.py`.
