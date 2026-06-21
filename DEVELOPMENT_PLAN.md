# Development Plan

## Repository Audit

- Current repository root: `protein_analysis_md` workspace at `Protein_Analysis_MD`.
- Historical package: `idrptm`.
- New public package wrapper: `protein_analysis_md`.
- New CLI: `pamd`; legacy `idrptm` and `idr-ptm-md` entry points remain.
- Existing tests cover import, CLI, design, prepare, PTM, cleavage, analysis, comparison, and reporting.
- CALVADOS is not vendored and is expected to be installed separately.

## Implementation Direction

1. Preserve existing `idrptm` behavior while exposing `protein_analysis_md`.
2. Keep CALVADOS execution isolated in `calvados_adapter.py` generated run scripts.
3. Implement sequence acquisition, PTM, cleavage, storage estimation, preparation, analysis, and reports as wrapper workflow layers.
4. Keep unit tests independent of internet access and CALVADOS.
5. Treat UniProt and CALVADOS smoke runs as optional integration steps.

## Current Known External Requirements

- Real MD execution requires a compatible CALVADOS installation, OpenMM, and a base CALVADOS residue parameter CSV.
- UniProt fetching requires network access and can be cached under `data/sequences/`.
- Short smoke simulations are pipeline checks, not scientifically meaningful sampling.
