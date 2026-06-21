from __future__ import annotations

import json
import sys

from idrptm.runner import execute_local_runs, plan_local_run, write_planned_status


def test_run_dry_status_and_execute(tmp_path) -> None:
    run_dir = tmp_path / "run"
    run_dir.mkdir()
    (run_dir / "config.yaml").write_text("{}\n", encoding="utf-8")
    (run_dir / "components.yaml").write_text("{}\n", encoding="utf-8")
    (run_dir / "run.py").write_text(
        "from pathlib import Path\n"
        "Path('completed.txt').write_text('ok\\n', encoding='utf-8')\n",
        encoding="utf-8",
    )

    plan = plan_local_run(run_dir, python_executable=sys.executable)
    status_path = write_planned_status(plan)
    planned = json.loads(status_path.read_text(encoding="utf-8"))

    assert planned["status"] == "planned"
    assert planned["phase"] == "all"

    (result,) = execute_local_runs((plan,))

    assert result.status == "completed"
    assert (run_dir / "completed.txt").read_text(encoding="utf-8") == "ok\n"
    completed = json.loads(result.status_json.read_text(encoding="utf-8"))
    assert completed["status"] == "completed"
    assert completed["elapsed_wall_s"] >= 0
    assert completed["execution_environment"]["hostname"]
    parameters_txt = (run_dir / "parameters.txt").read_text(encoding="utf-8")
    assert "paramdict" in parameters_txt
    assert "run_status.elapsed_wall_s" in parameters_txt
