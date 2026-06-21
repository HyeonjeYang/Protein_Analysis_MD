# protein_analysis_md

`protein_analysis_md` is a CALVADOS-based workflow framework for
coarse-grained protein/IDR MD setup, perturbation design, trajectory analysis,
and reporting. It is not an MD engine: CALVADOS and its simulation dependencies
must be installed separately.

## What This Repo Does

- Ingest protein/IDR sequences from direct input, FASTA, or optional UniProt use.
- Design WT, pSer, pThr, multi-protein, and pre-cleaved sequence states.
- Prepare CALVADOS-compatible run directories without modifying CALVADOS.
- Provide local/HPC execution scaffolding.
- Analyze trajectories for Rg, Ree, contacts, P(s), R(s), MSD, lifetimes, and
  related comparison metrics.
- Optionally run exploratory PCA and contact-environment eigen analyses.
- Generate WT-vs-perturbation reports and figures.
- Reports can include optional smoothed visual trends for noisy curves such as
  P(s), R(s), and energy, while preserving raw data and using raw values for
  default quantitative metrics.
- Provide environment and repository readiness diagnostics.

## What This Repo Does Not Do

- It does not implement or rewrite a molecular dynamics engine.
- It does not vendor, fork, patch, or redistribute CALVADOS or OpenMM source.
- It does not make generated trajectories scientifically sufficient by itself.
- It does not validate unsupported PTMs or residue parameters.
- It does not guarantee legal or licensing conclusions.

## Installation

```bash
python -m pip install -e ".[dev]"
```

Install CALVADOS separately according to the CALVADOS project instructions.
Generated trajectories and run outputs should not be committed to git.

## Minimal Usage

```bash
pamd compile configs/flk_smoke.yaml
pamd estimate-size configs/flk_smoke.yaml
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --phase all
pamd analyze runs/flk_smoke
pamd report runs/flk_smoke
```

Useful diagnostics:

```bash
pamd env-check
pamd repo-check
```

## Repository Structure

```text
configs/                  Small example workflow configs
docs/                     Design notes, scientific scope, and references
recipes/                  Small Python recipe examples
src/protein_analysis_md/  Public package namespace
src/idrptm/               Historical compatibility namespace
tests/                    Unit and CLI tests
```

Generated data should live in ignored directories such as `runs/`, `outputs/`,
`trajectories/`, `analysis_outputs/`, `reports/`, `logs/`, or local scratch
directories.

## License And Attribution

This repository is licensed GPL-3.0-only. CALVADOS is treated as an external
backend and is not vendored here. OpenMM is used through CALVADOS or optional
local environments. See `NOTICE.md`, `THIRD_PARTY_LICENSES.md`, `CITATION.cff`,
and `AUDIT.md` for attribution and audit notes.

Parts of this repository may have been developed with AI-assisted coding tools
and reviewed by the maintainer.

## Scientific Caution

Results are coarse-grained, model-dependent predictions. Short smoke runs are
pipeline tests, not scientific sampling. Production claims require appropriate
model choice, parameter provenance, convergence checks, controls, and domain
expert review.
