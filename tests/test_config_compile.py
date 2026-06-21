from __future__ import annotations

from textwrap import dedent

import yaml

from idrptm.configuration import compile_config_file
from idrptm.presets import merge_overrides, simulation_preset
from idrptm.schema import load_config


def test_short_config_compiles_to_locked_config(tmp_path) -> None:
    config = tmp_path / "short.yaml"
    outdir = tmp_path / "runs" / "short"
    config.write_text(
        dedent(
            f"""
            project:
              name: short
              outdir: {outdir}
            input:
              protein:
                source: direct
                name: seq
                sequence: "AST"
              ptm:
                mode: none
              cleavage:
                mode: none
            protocol:
              preset: smoke_single_chain
              simulation:
                production:
                  total_time_ns: 1.0
                  frame_interval_ns: 0.2
            analysis:
              preset: minimal
            report:
              preset: minimal
            execution:
              require_remote_for_md: true
              expected_hostname_contains: server
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    locked = compile_config_file(config)
    workflow = load_config(locked.project_dir)

    assert locked.lock_yaml.exists()
    assert locked.resolved_json.exists()
    assert workflow.project == "short"
    assert workflow.calvados.simulation.n_frames == 5
    assert workflow.calvados.simulation.save_every_steps == 20000
    assert workflow.execution.require_remote_for_md is True
    assert workflow.execution.expected_hostname_contains == "server"
    assert workflow.storage_estimate["runs"][0]["run_id"] == "seq__WT"


def test_legacy_explicit_config_compiles(tmp_path) -> None:
    config = tmp_path / "legacy.yaml"
    outdir = tmp_path / "runs" / "legacy"
    config.write_text(
        dedent(
            f"""
            project: legacy
            sequence:
              name: legacy_seq
              sequence: "AST"
            ptm:
              mode: wt
              include_wt: true
            calvados:
              platform: CPU
              simulation:
                total_time_ns: 0.2
                frame_interval_ns: 0.1
            runner:
              work_dir: {outdir}
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    locked = compile_config_file(config)
    workflow = load_config(locked.lock_yaml)

    assert workflow.project == "legacy"
    assert workflow.calvados.simulation.platform == "CPU"
    assert workflow.calvados.simulation.n_frames == 2


def test_preset_override_rules() -> None:
    preset = simulation_preset("smoke_single_chain")
    merged = merge_overrides(
        preset,
        {"production": {"total_time_ns": 2.0}, "replicates": 2},
    )

    assert merged["production"]["total_time_ns"] == 2.0
    assert merged["production"]["frame_interval_ns"] == 0.1
    assert merged["replicates"] == 2


def test_locked_yaml_is_plain_workflow_yaml(tmp_path) -> None:
    config = tmp_path / "plain.yaml"
    outdir = tmp_path / "runs" / "plain"
    config.write_text(
        dedent(
            f"""
            project:
              name: plain
              outdir: {outdir}
            input:
              protein:
                source: direct
                name: plain_seq
                sequence: "AST"
            protocol:
              preset: smoke_single_chain
            """
        ).strip()
        + "\n",
        encoding="utf-8",
    )

    locked = compile_config_file(config)
    payload = yaml.safe_load(locked.lock_yaml.read_text(encoding="utf-8"))

    assert payload["project"] == "plain"
    assert payload["compiled"]["simulation_preset"] == "smoke_single_chain"
