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
idrptm prepare configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan
idrptm prepare configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan --dry-run
idrptm run --config configs/example_ptm_scan.yaml --dry-run
idrptm analyze runs/example_ptm_scan/runs/phosphorylation_scan_fragment__WT
idrptm compare runs/example_ptm_scan
idrptm report runs/example_ptm_scan
```

`idrptm design` currently writes:

- `manifest.csv`
- one simulation-sequence FASTA file per variant under `fasta/`
- one per-run metadata stub under `runs/<variant_id>/metadata.yaml`

Design metadata preserves both `original_sequence` and `simulation_sequence`.
Configured residue positions use one-based biological numbering; internal
metadata also records zero-based indices.

`idrptm prepare` creates one CALVADOS run directory per manifest row. Each run
directory contains `input.fasta`, `residues.csv`, `config.yaml`,
`components.yaml`, `run.py`, and `metadata.json`.

The base CALVADOS residue CSV must be supplied either as
`calvados.residue_parameters` in the workflow config or through
`IDRPTM_CALVADOS_RESIDUES` / `CALVADOS_RESIDUES_CSV`. The source file is only
read; per-run residue tables are written into the generated run directories.

Simulation-time controls live under `calvados.simulation`. The MVP uses the
CALVADOS default integrator with `dt_ps: 0.01` only. You can specify either
`save_every_steps` plus `n_frames`, or `total_time_ns` plus
`frame_interval_ns`; the framework derives the CALVADOS `wfreq`, `steps`,
frame interval, and total simulated time and records those values in each
run's `metadata.json`.

The pure analysis core works on synthetic or trajectory-derived coordinate
arrays without importing CALVADOS. Implemented observables include Rg, Ree,
contact maps, P(s), internal-distance scaling, Flory exponent fitting, contact
lifetime, and center-of-mass MSD.

`idrptm analyze` expects a prepared CALVADOS run directory containing `top.pdb`
and `trajectory.dcd`, unless explicit paths are supplied with `--topology` and
`--trajectory`. It writes `timeseries_rg.parquet`, `timeseries_ree.parquet`,
`contact_map.npy`, `ps.parquet`, `scaling.parquet`, optional `msd.parquet`,
optional `contact_lifetime.parquet`, and `summary.json`.

`idrptm compare PROJECT_DIR` detects the WT condition by name and compares each
PTM condition against WT. It writes scalar summaries, aggregate P(s), delta P(s),
condition-average contact maps, and delta contact maps under `comparison/`.

`idrptm report PROJECT_DIR` writes `report/report.md` plus PNG and PDF figures
for Rg, Ree, contact maps, delta contact maps, P(s), R(s), PTM site annotation,
and the scalar summary table.

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
4. Stage 4: pure trajectory-analysis core for Rg, Ree, contacts, P(s), R(s),
   scaling, lifetime, and MSD.
5. Stage 5: CALVADOS trajectory loading and `idrptm analyze`.
6. Stage 6: WT-vs-PTM comparison tables, plots, and reports.
7. Stage 7: parameter validation workflow before expanding beyond pSer/pThr.
