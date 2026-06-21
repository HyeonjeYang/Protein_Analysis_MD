# protein_analysis_md

`protein_analysis_md` is a Python research workflow framework for sequence- and
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
src/protein_analysis_md/    Public Python package namespace
src/idrptm/                 Backward-compatible historical package
src/idrptm/analysis/        Analysis API
src/idrptm/plotting/        Plot/report API
configs/                    Example workflow YAML files
tests/                      Import and CLI smoke tests
```

## CLI

The package exposes `pamd` as the primary console script. Legacy `idrptm` and
`idr-ptm-md` entry points are retained as wrappers.

Scientific parameters should live in YAML configs or Python recipes. CLI flags
are reserved for execution controls such as `--force`, `--dry-run`, `--phase`,
and `--all-runs`.

### Minimal YAML

```bash
pamd compile configs/flk_smoke.yaml --force
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --all-runs --execute
pamd analyze runs/flk_smoke/runs/<RUN_ID>
pamd report runs/flk_smoke
```

`pamd compile` writes `project.lock.yaml`, `config_resolved.json`,
`storage_estimate.json`, and a manifest preview when sequence resolution is
available without interactive selection. Downstream commands can target the
compiled project directory and will use the lock file.

### Python Recipe

```bash
python recipes/flk_smoke.py
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --all-runs
```

Recipe files use `protein_analysis_md.recipe.Experiment` and `Protein` to define
the same scientific parameters as YAML.

### Interactive Wizard

```bash
pamd wizard
pamd compile configs/my_project.yaml
```

The wizard writes a config file only; it does not run MD.

```bash
pamd --help
pamd compile configs/flk_smoke.yaml --force
pamd search-uniprot FLK --reviewed --organism "Homo sapiens"
pamd fetch-sequence FLK --reviewed --organism "Homo sapiens" --interactive
pamd estimate-size runs/flk_smoke
pamd design configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan
pamd design configs/example_multi_protein.yaml --output-dir runs/example_multi_protein
pamd design configs/example_cleavage.yaml --output-dir runs/example_cleavage
pamd prepare configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan
pamd prepare configs/example_ptm_scan.yaml --output-dir runs/example_ptm_scan --dry-run
pamd run runs/example_ptm_scan --all-runs --dry-run
pamd status runs/example_ptm_scan
pamd resume runs/example_ptm_scan --dry-run
pamd analyze runs/example_ptm_scan/runs/phosphorylation_scan_fragment__WT
pamd compare runs/example_ptm_scan
pamd report runs/example_ptm_scan
```

`pamd design` currently writes:

- `manifest.csv`
- one simulation-sequence FASTA file per variant under `fasta/`
- one per-run metadata stub under `runs/<variant_id>/metadata.yaml`

Design metadata preserves both `original_sequence` and `simulation_sequence`.
Configured residue positions use one-based biological numbering; internal
metadata also records zero-based indices.

New workflows should use `proteins:` to define one or more protein/IDR sequence
inputs, each with its own PTM design settings. Legacy single-sequence configs
using `sequence:` plus top-level `ptm:` remain supported, and a single
`protein:` entry is also accepted. `system_sets:` defines explicit combinations
of protein PTM states, copy numbers, molecule types, and placement modes. This
supports single-protein runs, multi-protein runs, and mixed WT/PTM systems such
as WT and phosphorylated copies of the same protein in one CALVADOS system.

Proteolysis support is implemented as pre-simulation sequence-state design, not
dynamic bond breaking during MD. `cleavage_sets:` on a protein can generate an
intact state, individual fragment simulations, fragment-mixture simulations, and
sequential cleavage series. Built-in MVP rules include simple trypsin, Lys-C,
Arg-C, high-specificity chymotrypsin, CNBr, and TEV-style cleavage. Design
outputs include `cleavage_sites.csv`, `fragments.fasta`, and
`cleavage_manifest.csv` with original residue ranges and preserved PTM mapping.

`pamd prepare` creates one CALVADOS run directory per manifest row. Each run
directory contains `input.fasta`, `residues.csv`, `config.yaml`,
`components.yaml`, `run.py`, and `metadata.json`.

For multi-component systems, `components.yaml` contains one CALVADOS component
per configured system component. If `calvados.topol` is left as `center`, the
adapter uses the `system_sets[].placement` value, such as `grid`, `random`, or
`slab`, rather than centering all components on top of each other.

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

When per-residue chain metadata are available, analysis also supports
chain-resolved Rg/Ree, intra-chain and inter-chain contact maps, chain COM
distances, per-chain COM MSD, cluster-size time series, and inter-protein
contact lifetime.

For cleaved systems, helper analysis APIs cover fragment-resolved Rg/Ree,
inter-fragment contacts, fragment cluster size, intact-vs-cleaved delta contact
maps, and contact maps projected back onto original sequence coordinates.

`pamd analyze` expects a prepared CALVADOS run directory containing `top.pdb`
and `trajectory.dcd`, unless explicit paths are supplied with `--topology` and
`--trajectory`. It writes `timeseries_rg.parquet`, `timeseries_ree.parquet`,
`contact_map.npy`, `ps.parquet`, `scaling.parquet`, optional `msd.parquet`,
optional `contact_lifetime.parquet`, and `summary.json`.

Internal analysis units are explicit: coordinates and distance observables are
stored in nm, timestep metadata in ps, report times in ns, energies in kJ/mol,
temperature in K, ionic strength in M, charge in elementary charge, mass in amu,
MSD in nm^2, and contact probabilities/P(s) as dimensionless values. Each
analysis artifact is accompanied by a `*.units.json` sidecar, and every
`summary.json` contains a `units` block. The MDAnalysis trajectory reader
converts Angstrom coordinates to nm before analysis and records both the input
and canonical coordinate units.

`pamd compare PROJECT_DIR` detects the WT condition by name and compares each
PTM condition against WT. It writes scalar summaries, aggregate P(s), delta P(s),
condition-average contact maps, and delta contact maps under `comparison/`.

`pamd report PROJECT_DIR` writes `report/report.md` plus PNG and PDF figures
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
5. Stage 5: CALVADOS trajectory loading and `pamd analyze`.
6. Stage 6: WT-vs-PTM comparison tables, plots, and reports.
7. Stage 7: parameter validation workflow before expanding beyond pSer/pThr.
