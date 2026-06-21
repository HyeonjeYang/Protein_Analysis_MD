# FLK Smoke Run

Recommended commands:

```bash
pamd search-uniprot FLK --reviewed --organism "Homo sapiens"
pamd fetch-sequence FLK --reviewed --organism "Homo sapiens" --interactive
pamd estimate-size configs/flk_smoke.yaml
pamd design configs/flk_smoke.yaml
pamd prepare configs/flk_smoke.yaml
pamd run runs/flk_smoke/<RUN_ID> --phase all
pamd analyze runs/flk_smoke/<RUN_ID>
pamd compare runs/flk_smoke
pamd report runs/flk_smoke
```

FLK is ambiguous. Select a reviewed human candidate explicitly before using the
sequence for a smoke run.
