# idr-ptm-md

`idr-ptm-md` is a Python research workflow framework for sequence- and
post-translational-modification-dependent coarse-grained molecular dynamics of
intrinsically disordered regions and proteins.

The project is intentionally a workflow layer. It prepares sequence/PTM designs,
CALVADOS-compatible run directories, local or HPC execution scaffolds, trajectory
analysis, WT-vs-PTM comparisons, and report/figure outputs. It does not rewrite a
molecular dynamics engine.

## CALVADOS Backend

Simulations must use [KULL-Centre/CALVADOS](https://github.com/KULL-Centre/CALVADOS)
as the backend. This repository does not vendor, fork, patch, or modify upstream
CALVADOS source code. Users should install CALVADOS separately and cite CALVADOS
and the relevant CALVADOS model/software papers for simulation results.

## MVP Scientific Scope

Stage 1 is a skeleton only. The intended MVP scientific scope is:

- single-chain IDR/protein sequence simulations
- WT vs phosphorylated sequence variants
- supported PTMs: pSer and pThr only
- planned observables: Rg, Ree, contact map, P(s), R(s), scaling exponent,
  contact lifetime, and MSD

## Scientific Limitations

- No LLPS or multi-chain production pipeline in the MVP.
- No pTyr, acetylation, methylation, ubiquitination, or other PTMs until
  parameters are explicitly added and validated.
- Stage 1 commands are placeholders that parse arguments and print clear intent.
- The framework does not validate new residue parameters or scientific
  conclusions by itself.

## Repository Layout

```text
src/idrptm/                 Python package
src/idrptm/analysis/        Analysis API placeholders
src/idrptm/plotting/        Plot/report API placeholders
configs/                    Example workflow YAML files
tests/                      Import and CLI smoke tests
```

## CLI

The package exposes both `idrptm` and `idr-ptm-md` console scripts:

```bash
idrptm --help
idrptm init --output work/example
idrptm design configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan
idrptm prepare --config configs/example_ptm_scan.yaml --output-dir runs
idrptm run --config configs/example_ptm_scan.yaml --dry-run
idrptm analyze --config configs/example_ptm_scan.yaml
idrptm compare --reference wt.csv --variant ptm.csv
idrptm report --config configs/example_ptm_scan.yaml
```

`idrptm design` currently writes:

- `manifest.csv`
- one simulation-sequence FASTA file per variant under `fasta/`
- one per-run metadata stub under `runs/<variant_id>/metadata.yaml`

Design metadata preserves both `original_sequence` and `simulation_sequence`.
Configured residue positions use one-based biological numbering; internal
metadata also records zero-based indices.

## Development

```bash
python -m pip install -e ".[dev]"
python -m pytest
ruff check .
```

## Roadmap

1. Stage 1: repository skeleton, typed package layout, config examples, CLI
   placeholders, citation/license files, and smoke tests.
2. Stage 2: sequence ingestion, PTM state generation for pSer/pThr, and
   CALVADOS run-directory generation.
3. Stage 3: local and SLURM execution scaffolds with reproducible manifests.
4. Stage 4: trajectory analysis for Rg, Ree, contacts, P(s), R(s), scaling,
   lifetime, and MSD.
5. Stage 5: WT-vs-PTM comparison tables, plots, and reports.
6. Stage 6: parameter validation workflow before expanding beyond pSer/pThr.
