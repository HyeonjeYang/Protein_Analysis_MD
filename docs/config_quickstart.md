# Config Quickstart

Minimal short-config workflow:

```bash
pamd compile configs/flk_smoke.yaml
pamd estimate-size runs/flk_smoke
pamd prepare runs/flk_smoke
pamd run runs/flk_smoke --all-runs --phase all --execute
pamd analyze runs/flk_smoke/runs/<RUN_ID> --config runs/flk_smoke/project.lock.yaml
pamd compare runs/flk_smoke
pamd report runs/flk_smoke
```

Use `execution.require_remote_for_md: true` when a project must only run on a
Remote-SSH/server host:

```yaml
execution:
  require_remote_for_md: true
  expected_hostname_contains: "server"
```

Use `pamd env-check` before launching long simulations to confirm where commands
are executing.
