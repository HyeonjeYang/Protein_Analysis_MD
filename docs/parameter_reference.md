# Parameter Reference

Use concise config blocks:

- `project`: name and output directory.
- `input`: protein, PTM, mutation, and cleavage definitions.
- `protocol`: simulation preset plus environment/simulation/storage overrides.
- `analysis`: analysis preset plus overrides.
- `report`: report preset and output format controls.
- `sweep`: optional dimensions for planned parameter sweeps.

Downstream lock files are explicit `WorkflowConfig` YAML files for reproducible
execution.
