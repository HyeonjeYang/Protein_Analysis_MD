# Config Architecture

`protein_analysis_md` supports concise user configs and explicit locked configs.
Scientific parameters live in YAML or Python recipes. CLI flags are reserved for
execution controls such as `--force`, `--dry-run`, `--phase`, and `--all-runs`.

Run `pamd compile CONFIG.yaml` to normalize inputs, resolve presets, validate the
workflow, derive CALVADOS timing, and write `project.lock.yaml` plus
`config_resolved.json`. Downstream commands can then target the project directory.
