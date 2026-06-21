# v0.1.0 Release Notes

## Tag

- Git tag: `v0.1.0`
- Release title: `protein_analysis_md v0.1.0`
- Status: pre-alpha research workflow framework

## Short Description

`protein_analysis_md` v0.1.0 is a GPL-3.0-only Python workflow framework for
CALVADOS-backed coarse-grained protein/IDR MD setup, execution scaffolding,
trajectory analysis, WT-vs-perturbation comparison, and report generation.

CALVADOS and OpenMM are external dependencies. This repository does not vendor,
fork, patch, or redistribute upstream CALVADOS/OpenMM source code.

## Suggested GitHub Topics

- `calvados`
- `molecular-dynamics`
- `coarse-grained`
- `protein-analysis`
- `intrinsically-disordered-proteins`
- `phosphorylation`
- `openmm`
- `python`

## Release Body

```markdown
## protein_analysis_md v0.1.0

Initial public pre-alpha release.

Highlights:
- CALVADOS-compatible run-directory preparation without modifying CALVADOS.
- WT, pSer, pThr, multi-protein, and pre-cleavage workflow scaffolds.
- Local/HPC execution helpers, tmux launch, progress watch, and finalize.
- Analysis, comparison, report figures, PyMOL export, and local HTML dashboard.
- GPL-3.0-only license, citation metadata, audit notes, and CI.

Scientific caution:
This is a research prototype. Short runs are pipeline checks, not scientific
sampling. Production conclusions require appropriate model choice, parameter
provenance, convergence checks, controls, and domain expert review.
```

## Release Assets

Upload or rely on:

- GitHub auto-generated `Source code (zip)`.
- GitHub auto-generated `Source code (tar.gz)`.
- Optional Python artifacts if built locally:
  - `dist/protein_analysis_md-0.1.0.tar.gz`
  - `dist/protein_analysis_md-0.1.0-py3-none-any.whl`

Do not upload:

- `runs/`
- `work/`
- `.venv/`
- DCD/XTC/TRR trajectories
- PDB/checkpoint/log/report/dashboard outputs
- generated figures/tables
- user-provided CALVADOS residue tables unless redistribution rights are
  explicitly reviewed

## Optional Build Commands

```bash
python -m pip install build
python -m build
```

Run checks before tagging:

```bash
python -m ruff check src tests
python -m pytest -q
pamd repo-check
```
