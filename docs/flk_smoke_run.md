# FLK Smoke Run

Recommended commands:

```bash
pamd search-uniprot FLK --reviewed --organism "Homo sapiens"
pamd fetch-sequence FLK --reviewed --organism "Homo sapiens" --interactive
pamd compile configs/flk_smoke.yaml
pamd estimate-size runs/flk_smoke
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke/runs/<RUN_ID> --phase all --execute
pamd analyze runs/flk_smoke/runs/<RUN_ID> --config runs/flk_smoke/project.lock.yaml
pamd compare runs/flk_smoke
pamd report runs/flk_smoke
```

FLK is ambiguous. Select a reviewed human candidate explicitly before using the
sequence for a smoke run.
