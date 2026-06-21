# Repository Check

This report is a practical readiness scan, not legal advice. It warns about
possible public-release issues but does not claim legal certainty.

## Required Files

- LICENSE: present
- NOTICE.md: present
- THIRD_PARTY_LICENSES.md: present
- CITATION.cff: present
- README.md lines: 88
- README <= 120 lines: True

## Gitignore

- Missing required patterns: none

## Large Files

### Files > 5 MB

- none

### Files > 50 MB

- none


## Generated Or MD Output Files

### Generated files found in working tree

- none

### Tracked generated files

- none


## Possible Upstream/License References

### Reference hits

- `LICENSE: GPL`
- `CITATION.cff: CALVADOS`
- `pyproject.toml: CALVADOS`
- `REPO_CHECK.md: CALVADOS`
- `README.md: CALVADOS`
- `AUDIT.md: CALVADOS`
- `.gitignore: CALVADOS`
- `THIRD_PARTY_LICENSES.md: CALVADOS`
- `NOTICE.md: CALVADOS`
- `tests/test_smoothing.py: Config(`
- `tests/test_smoothing_policy.py: Config(`
- `tests/test_recipe.py: Simulation`
- `tests/test_imports.py: CALVADOS`
- `tests/test_residue_params.py: CALVADOS`
- `tests/test_cleavage.py: Config(`
- `tests/test_ptm.py: Simulation`
- `tests/test_storage.py: CALVADOS`
- `tests/test_environment_check.py: CALVADOS`
- `tests/test_simulation_config.py: Config(`
- `tests/test_runner.py: Components`
- `tests/test_analysis_io_pipeline.py: CALVADOS`
- `tests/test_decomposition.py: Config(`
- `tests/test_config_compile.py: CALVADOS`
- `tests/test_prepare.py: CALVADOS`
- `tests/test_design.py: CALVADOS`
- `tests/test_multi_protein.py: CALVADOS`
- `docs/cleavage_models.md: Simulation`
- `docs/architecture.md: CALVADOS`
- `docs/parameter_reference.md: Simulation`
- `docs/limitations.md: CALVADOS`
- `docs/phase_separation.md: CALVADOS`
- `docs/license_notes.md: CALVADOS`
- `docs/phase_separation_plots.md: Simulation`
- `docs/config_architecture.md: CALVADOS`
- `docs/config_quickstart.md: Simulation`
- `docs/scientific_scope.md: Simulation`
- `docs/staged_temporal_cleavage.md: Simulation`
- `docs/presets.md: Simulation`
- `docs/python_recipes.md: Simulation`
- `recipes/flk_smoke.py: Simulation`
- `configs/example_ptm_scan.yaml: CALVADOS`
- `configs/example_cleavage_manual.yaml: CALVADOS`
- `configs/example_single_idr.yaml: CALVADOS`
- `configs/example_cleavage_enzyme.yaml: CALVADOS`
- `configs/example_cleavage_poisson.yaml: CALVADOS`
- `configs/example_phase_separation.yaml: CALVADOS`
- `configs/example_multi_protein.yaml: CALVADOS`
- `configs/flk_smoke.yaml: CALVADOS`
- `configs/example_cleavage.yaml: CALVADOS`
- `data/residues/README.md: CALVADOS`
- `data/sequences/Q99895_CTRC_HUMAN.fasta: GPL`
- `data/sequences/Q99895_CTRC_HUMAN.metadata.json: GPL`
- `src/protein_analysis_md.egg-info/PKG-INFO: CALVADOS`
- `src/protein_analysis_md.egg-info/SOURCES.txt: CALVADOS`
- `src/idrptm/configuration.py: CALVADOS`
- `src/idrptm/runner.py: CALVADOS`
- `src/idrptm/cleavage.py: Simulation`
- `src/idrptm/ptm.py: CALVADOS`
- `src/idrptm/presets.py: CALVADOS`
- `src/idrptm/calvados_adapter.py: CALVADOS`
- `src/idrptm/design.py: CALVADOS`
- `src/idrptm/__init__.py: CALVADOS`
- `src/idrptm/hpc.py: Config(`
- `src/idrptm/residue_params.py: CALVADOS`
- `src/idrptm/cli.py: CALVADOS`
- `src/idrptm/storage.py: CALVADOS`
- `src/idrptm/environment.py: Simulation`
- `src/idrptm/recipe.py: Config(`
- `src/idrptm/schema.py: CALVADOS`
- `src/idrptm/analysis/cleavage.py: Simulation`
- `src/idrptm/analysis/io.py: CALVADOS`
- `src/idrptm/analysis/__init__.py: CALVADOS`
- `src/idrptm/analysis/energy.py: CALVADOS`
- `src/idrptm/analysis/pipeline.py: CALVADOS`
- `src/idrptm/analysis/decomposition.py: Config(`
- `src/idrptm/analysis/compare.py: Components`
- `src/idrptm/plotting/plots.py: Components`
- `src/protein_analysis_md/environment_check.py: CALVADOS`
- `src/protein_analysis_md/calvados_adapter.py: CALVADOS`
- `src/protein_analysis_md/__init__.py: CALVADOS`
- `src/protein_analysis_md/repo_check.py: CALVADOS`


## Secret Scan

### Potential secret-pattern hits

- none


## Import Check

- protein_analysis_md: ok
- idrptm: ok

## Suggested Next Command

```bash
pytest
ruff check .
```
