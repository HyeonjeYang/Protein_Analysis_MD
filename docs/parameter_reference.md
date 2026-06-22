# Parameter Reference

Use concise config blocks:

- `project`: name and output directory.
- `input`: protein, PTM, mutation, and cleavage definitions.
- `protocol`: simulation preset plus environment/simulation/storage overrides.
- `analysis`: analysis preset plus overrides.
- `report`: report preset and output format controls. Standard reports write
  PNG figures by default; set `PAMD_FIGURE_FORMATS=png,pdf` for PDF copies.
- `sweep`: optional dimensions for planned parameter sweeps.

Downstream lock files are explicit `WorkflowConfig` YAML files for reproducible
execution.
