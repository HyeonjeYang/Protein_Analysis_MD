# protein_analysis_md

Status: **v0.1.1 pre-alpha research framework**.

`protein_analysis_md` is a Python workflow layer for sequence-, PTM-, and
cleavage-aware coarse-grained protein/IDR molecular dynamics. It prepares,
runs, analyzes, compares, and reports CALVADOS-compatible simulations.

Developed by Present0206 with implementation help from Codex.

## Scope

- Uses KULL-Centre/CALVADOS as an external simulation backend.
- Does not implement, vendor, fork, patch, or redistribute CALVADOS/OpenMM.
- Supports MVP phosphorylation states: WT, pSer, and pThr.
- Supports pre-cleaved sequence states and fragment-mixture setup.
- Includes analysis and figures for Rg, Ree, contacts, P(s), R(s), `<R^2(s)>`, MSD,
  lifetimes, PCA/decomposition, and comparison workflows.
- Provides local/HPC scaffolding, tmux-friendly local launch, PyMOL export, and
  a local static HTML dashboard.

This is a research prototype. Short runs are pipeline checks, not scientific
sampling. Production conclusions require appropriate model choice, parameter
provenance, convergence checks, controls, and domain expert review.

## Install

```bash
git clone https://github.com/HyeonjeYang/protein_analysis_md.git
cd protein_analysis_md
conda env create -f environment.yml
conda activate protein_analysis_md
python -m pip install -e ".[dev]"
pamd env-check
```

For analysis-only development, `python -m pip install -e ".[dev]"` is enough.
For simulation, install CALVADOS/OpenMM in the active environment and keep their
license/citation requirements.

## Quick Start

```bash
pamd compile configs/flk_smoke.yaml
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --all-runs --phase all --execute
pamd analyze runs/flk_smoke/runs/<RUN_ID> --config runs/flk_smoke/project.lock.yaml
pamd compare runs/flk_smoke
pamd report runs/flk_smoke
pamd dashboard runs/flk_smoke --open
```

Useful commands:

```bash
pamd watch runs/<PROJECT_DIR> --follow
pamd launch-local runs/<PROJECT_DIR> --backend auto --terminal tmux
pamd finalize runs/<PROJECT_DIR>
pamd pymol runs/<PROJECT_DIR>
pamd pack runs/<PROJECT_DIR>
pamd repo-check
```

Reports write PNG figures by default. Set `PAMD_FIGURE_FORMATS=png,pdf` before
`pamd report` or `pamd finalize` if PDF copies are needed.

## Documentation

- `docs/config_quickstart.md` - concise config workflow
- `docs/parameter_reference.md` - key parameters
- `docs/scientific_scope.md` and `docs/limitations.md` - scientific boundaries
- `docs/visualization_gallery.md` - figure and report outputs
- `CHANGELOG.md` - release history

Generated trajectories, figures, reports, logs, dashboards, and local
environments are ignored by git and should not be committed.

## Repository Layout

```text
configs/                  Example workflow configs
docs/                     User and developer notes
recipes/                  Small Python recipe examples
src/protein_analysis_md/  Public package namespace
src/idrptm/               Historical compatibility namespace
tests/                    Unit and CLI tests
```

## License And Citation

This repository is GPL-3.0-only. CALVADOS and OpenMM remain external
dependencies. See `LICENSE`, `NOTICE.md`, `THIRD_PARTY_LICENSES.md`,
`CITATION.cff`, and `AUDIT.md`.
