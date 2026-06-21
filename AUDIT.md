# Repository Audit

This audit is a practical public-release cleanup review for
`protein_analysis_md`. It is not legal advice. When licensing or provenance is
uncertain, the uncertainty is recorded here rather than treated as resolved.

## Current Structure Summary

- Package/distribution name: `protein-analysis-md` in `pyproject.toml`.
- Public import namespace: `protein_analysis_md`.
- Historical/backward-compatible namespace: `idrptm`.
- License metadata: `GPL-3.0-only` in `pyproject.toml`.
- License file: `LICENSE` exists and contains the GNU GPL v3 text.
- Source directories:
  - `src/idrptm/`: current implementation body.
  - `src/protein_analysis_md/`: public compatibility wrappers and entry points.
- Tests directory: `tests/`.
- Config examples: `configs/`.
- Documentation: `docs/`.
- Recipes: `recipes/`.
- Generated/local output directories observed:
  - `runs/` with chymotrypsin smoke/long outputs.
  - `work/` with local CALVADOS residue CSV, local YAMLs, and a virtual
    environment-like `work/calvados_venv/`.
  - `.pytest_cache/`, `.ruff_cache/`, `.venv/`, and generated egg-info.

## Compact Public Tree

The public tree should focus on:

```text
.
├── CITATION.cff
├── LICENSE
├── NOTICE.md
├── README.md
├── THIRD_PARTY_LICENSES.md
├── configs/
├── docs/
├── environment.yml
├── pyproject.toml
├── recipes/
├── src/
│   ├── idrptm/
│   └── protein_analysis_md/
└── tests/
```

`runs/`, `work/`, caches, trajectories, checkpoints, analysis outputs, and
large binary files should remain untracked and ignored.

## Possible Upstream-Derived Files

No tracked file was identified as a direct vendored CALVADOS or OpenMM source
file during this audit. The tracked implementation appears to generate
CALVADOS-compatible inputs rather than copying CALVADOS source code.

Files and areas that intentionally reference CALVADOS/OpenMM concepts and
should remain reviewed for provenance:

- `src/idrptm/calvados_adapter.py`: generates `run.py`, `config.yaml`,
  `components.yaml`, and `residues.csv` for prepared runs. It should continue
  to generate minimal adapter scaffolding rather than vendoring upstream files.
- `src/idrptm/residue_params.py`: writes per-run residue parameters and refers
  to CALVADOS/pCALVADOS-style residue tables. User-provided CALVADOS residue
  CSVs should not be committed.
- `src/idrptm/schema.py`: maps framework timing fields to CALVADOS-like
  `wfreq`, `steps`, `runtime`, platform, and restart settings.
- `tests/test_prepare.py` and related tests: validate generated file shapes
  but do not appear to include upstream CALVADOS code.

Generated run directories under `runs/` contain CALVADOS-like file names
(`run.py`, `config.yaml`, `components.yaml`, `residues.csv`) and trajectory
outputs. They are ignored and should not be committed.

## Large And Generated Files Found

Tracked files larger than 5 MB: none found.

Local untracked/ignored large files larger than 5 MB were found under
`work/calvados_venv/`, including dynamic libraries and OpenMM data files. This
directory is local environment state and should stay ignored.

Generated MD artifacts were found under `runs/`, including:

- DCD trajectories.
- OpenMM/CALVADOS-style restart checkpoints.
- Log files.
- Parquet analysis outputs.
- NumPy contact-map outputs.
- SQLite project databases.
- Report/comparison outputs.

These files are user data or generated outputs. They were not deleted.

## License Status

- Current repository license is GPL-3.0-only.
- This is a conservative choice because the framework is built around
  CALVADOS workflows and may generate CALVADOS-compatible scaffolding, while
  CALVADOS itself remains an external dependency.
- CALVADOS and OpenMM are not vendored by this repository.
- Users should install CALVADOS separately and follow CALVADOS/OpenMM citation
  and license requirements for their own environment.
- The exact license obligations for generated run scripts and user-provided
  CALVADOS residue parameter tables are not determined here. Keep generated
  runs and upstream parameter tables out of git unless their provenance and
  redistribution terms are explicitly reviewed.

## Recommended Cleanup Steps

1. Keep `LICENSE`, `NOTICE.md`, `THIRD_PARTY_LICENSES.md`, and `CITATION.cff`.
2. Keep GPL-3.0-only unless a maintainer performs a legal review supporting a
   different license.
3. Keep CALVADOS/OpenMM as external dependencies.
4. Do not commit upstream CALVADOS source, OpenMM source/data, generated run
   directories, trajectories, checkpoints, logs, analysis outputs, or local
   virtual environments.
5. Update `.gitignore` to cover common MD outputs, generated arrays/tables,
   local environments, logs, and caches.
6. Shorten `README.md`; move long design notes to `docs/`.
7. Add `pamd env-check` so users can see whether commands are executing on a
   local laptop or a Remote-SSH host.
8. Add `pamd repo-check` for repeatable public-release checks.
9. Consider migrating implementation modules from `src/idrptm/` to
   `src/protein_analysis_md/` in a compatibility-preserving way. Until then,
   document `idrptm` as a historical namespace.

## Uncertainties

- This audit did not compare repository text against the full current CALVADOS
  repository or OpenMM source tree. It is a local repository audit only.
- The provenance of local ignored files under `runs/` and `work/` was not
  independently verified.
- The source of `work/calvados_residues.csv` is not reviewed here. It should
  remain ignored unless redistribution rights are confirmed.
- Some generated `run.py` templates may resemble common CALVADOS usage
  patterns because they call the CALVADOS API. That does not by itself prove
  copying, but it should remain reviewed if template text grows closer to
  upstream examples.
- No legal conclusion is made about derivative-work status. GPL-3.0-only is
  retained as a conservative project policy.
