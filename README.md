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

### Download The Repository

Clone with Git:

```bash
git clone https://github.com/HyeonjeYang/Protein_Analysis_MD.git
cd Protein_Analysis_MD
```

Or download the GitHub ZIP archive, unzip it, and open a terminal in the
repository root.

### Recommended Conda Install

The same `environment.yml` is intended for Linux, macOS, Windows PowerShell, and
Windows WSL2. It installs OpenMM and CALVADOS as external dependencies; their
source is not copied into this repository.

```bash
conda env create -f environment.yml
conda activate protein_analysis_md
python -m pip install -e ".[dev]"
pamd env-check
```

If CALVADOS is installed in the active environment, `pamd prepare` will try to
use `calvados/data/residues.csv` automatically. To force a specific residue
table, set `IDRPTM_CALVADOS_RESIDUES`:

```bash
export IDRPTM_CALVADOS_RESIDUES=/path/to/calvados/data/residues.csv
```

PowerShell:

```powershell
$env:IDRPTM_CALVADOS_RESIDUES = "C:\path\to\calvados\data\residues.csv"
```

### Linux

Use the recommended conda install above. On shared Linux/HPC systems, create the
environment in your home or project scratch space, run `pamd env-check`, and use
`pamd hpc-script` only after `pamd prepare` has produced run directories.

### macOS

Use the recommended conda install above on Intel or Apple Silicon Macs. Local
OpenMM platforms may include `CPU` and sometimes `OpenCL`; CUDA is not expected
on ordinary Macs. For long local runs, keep the Mac awake and plugged in.

### Windows

WSL2 with Ubuntu is the preferred Windows route for long CALVADOS runs:

```powershell
wsl --install -d Ubuntu
```

Then open Ubuntu and use the Linux commands above.

Native Windows PowerShell can also run the workflow if conda, OpenMM, and
CALVADOS install successfully:

```powershell
conda env create -f environment.yml
conda activate protein_analysis_md
python -m pip install -e ".[dev]"
pamd env-check
```

For production MD, prefer WSL2 or a Linux server/HPC if native Windows OpenMM or
CALVADOS dependencies are unstable.

### Pip-Only Development Install

For analysis-only development without CALVADOS/OpenMM execution:

```bash
python -m pip install -e ".[dev]"
```

Generated trajectories and run outputs should not be committed to git.

## Minimal Usage

```bash
pamd compile configs/flk_smoke.yaml
pamd estimate-size runs/flk_smoke
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --all-runs --phase all --execute
pamd analyze runs/flk_smoke/runs/<RUN_ID> --config runs/flk_smoke/project.lock.yaml
pamd compare runs/flk_smoke
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
