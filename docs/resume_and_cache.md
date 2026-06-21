# Resume And Cache

`pamd status PROJECT_DIR` summarizes prepared, completed, failed, analyzed, and
checkpointed runs. `pamd resume PROJECT_DIR` plans or executes incomplete runs.
`pamd clean PROJECT_DIR --yes` removes derived cache/temp files only.

Analysis cache metadata is stored in `analysis/cache_manifest.json` and is keyed
by trajectory file metadata, analysis config, software version, and units version.
