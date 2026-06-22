#!/usr/bin/env python3
"""Unified simulation launcher.

Generate a config YAML from parameters, compile, prepare, and launch.

Examples
--------
  python run_simulation.py --reps 25 --total-time-ns 400
  python run_simulation.py --reps 50 --total-time-ns 100 --name chymo_test
  python run_simulation.py --reps 10 --total-time-ns 50 --dry-run
"""

from __future__ import annotations

import argparse
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import yaml

REPO_ROOT = Path(__file__).resolve().parent
CONFIGS_DIR = REPO_ROOT / "configs"
PAMD = REPO_ROOT / ".venv" / "bin" / "pamd"
DEFAULT_FASTA = REPO_ROOT / "data" / "sequences" / "Q99895_CTRC_HUMAN.fasta"
DEFAULT_PROTEIN_NAME = "CTRC_HUMAN_Q99895"


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Generate a config and launch a CALVADOS simulation pipeline.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--name", default=None, help="Project name (auto-generated if omitted)")
    p.add_argument("--fasta", type=Path, default=DEFAULT_FASTA, help="FASTA file path")
    p.add_argument("--protein-name", default=None, help="Internal protein identifier (auto from FASTA stem)")
    p.add_argument("--reps", type=int, default=25, help="Number of replicates")
    p.add_argument("--total-time-ns", type=float, default=100.0, help="Production simulation time (ns)")
    p.add_argument("--frame-interval-ns", type=float, default=None, help="Trajectory frame interval (ns, auto if omitted)")
    p.add_argument("--equil-time-ns", type=float, default=10.0, help="Equilibration time (ns)")
    p.add_argument("--temperature-k", type=float, default=298.0, help="Temperature (K)")
    p.add_argument("--ph", type=float, default=7.4, help="pH")
    p.add_argument("--ionic-m", type=float, default=0.15, help="Ionic strength (M)")
    p.add_argument("--seed", type=int, default=None, help="Random seed (YYYYMMDD today if omitted)")
    p.add_argument("--simulation-parallel", default="auto", help="Concurrent simulations")
    p.add_argument("--analysis-parallel", default="auto", help="Concurrent analyses")
    p.add_argument("--session-name", default=None, help="tmux/byobu session name")
    p.add_argument("--force", action="store_true", help="Overwrite existing config/lock files")
    p.add_argument("--dry-run", action="store_true", help="Print commands without executing")
    return p.parse_args()


def auto_frame_interval(total_time_ns: float) -> float:
    """Target ~1000 frames, clamped to [0.05, 0.5] ns."""
    raw = total_time_ns / 1000.0
    return round(max(0.05, min(0.5, raw)), 3)


def build_config(args: argparse.Namespace) -> tuple[str, dict]:
    fasta = args.fasta.resolve()
    fasta_rel = Path("..") / fasta.relative_to(REPO_ROOT)

    protein_name = args.protein_name or fasta.stem
    total_ns = args.total_time_ns
    name = args.name or f"{protein_name}_{args.reps}rep_{int(total_ns)}ns"
    frame_interval = args.frame_interval_ns or auto_frame_interval(total_ns)
    seed = args.seed or int(datetime.now(timezone.utc).strftime("%Y%m%d"))

    config = {
        "project": {
            "name": name,
            "outdir": f"runs/{name}",
        },
        "input": {
            "protein": {
                "source": "fasta",
                "name": protein_name,
                "fasta": str(fasta_rel),
                "charge_termini": "both",
            },
            "ptm": {"mode": "none"},
            "cleavage": {"mode": "none"},
        },
        "protocol": {
            "preset": "production_single_chain",
            "simulation": {
                "replicates": args.reps,
                "model": "CALVADOS2",
                "box_nm": [50.0, 50.0, 50.0],
                "platform": "CPU",
                "runtime_hours": 0,
                "random_seed": seed,
                "production": {
                    "total_time_ns": total_ns,
                    "frame_interval_ns": frame_interval,
                },
                "equilibration": {
                    "total_time_ns": args.equil_time_ns,
                    "frame_interval_ns": 1.0,
                    "save_trajectory": False,
                },
            },
            "environment": {
                "temperature_K": args.temperature_k,
                "pH": args.ph,
                "ionic_M": args.ionic_m,
            },
        },
        "analysis": {
            "preset": "standard_idr",
            "overrides": {
                "contact_cutoff_nm": 1.0,
                "min_sequence_separation": 2,
                "max_lag": 200,
                "free_energy": {
                    "enabled": True,
                    "variables": [["Rg", "Ree"]],
                    "bins": 50,
                    "temperature_K": args.temperature_k,
                    "min_count": 1,
                },
            },
        },
        "report": {"preset": "standard"},
        "execution": {"require_remote_for_md": False},
    }
    return name, config


def run_cmd(cmd: list[str], dry_run: bool) -> None:
    print("$", " ".join(str(c) for c in cmd))
    if not dry_run:
        result = subprocess.run(cmd)
        if result.returncode != 0:
            sys.exit(result.returncode)


def main() -> None:
    args = parse_args()
    name, config = build_config(args)

    config_path = CONFIGS_DIR / f"{name}.yaml"
    project_dir = REPO_ROOT / config["project"]["outdir"]
    sim = config["protocol"]["simulation"]

    print(f"Project  : {name}")
    print(f"Config   : {config_path}")
    print(f"Out      : {project_dir}")
    print(f"Reps     : {sim['replicates']}")
    print(f"Time     : {sim['production']['total_time_ns']} ns  (interval {sim['production']['frame_interval_ns']} ns)")
    print(f"Equil    : {sim['equilibration']['total_time_ns']} ns")
    print(f"Seed     : {sim['random_seed']}")
    print(f"Parallel : sim={args.simulation_parallel}  analysis={args.analysis_parallel}")
    print()

    if not args.dry_run:
        if config_path.exists() and not args.force:
            print(f"Config already exists: {config_path}  (use --force to overwrite)")
            sys.exit(1)
        CONFIGS_DIR.mkdir(parents=True, exist_ok=True)
        config_path.write_text(yaml.safe_dump(config, sort_keys=False), encoding="utf-8")
        print(f"Wrote {config_path}\n")

    force_flag = ["--force"] if args.force else []

    run_cmd([str(PAMD), "compile", str(config_path)] + force_flag, args.dry_run)
    run_cmd([str(PAMD), "prepare", str(config_path)], args.dry_run)

    launch_cmd = [
        str(PAMD), "launch-local", str(project_dir),
        "--simulation-parallel", str(args.simulation_parallel),
        "--analysis-parallel", str(args.analysis_parallel),
        "--yes",
    ]
    if args.session_name:
        launch_cmd += ["--session-name", args.session_name]
    run_cmd(launch_cmd, args.dry_run)


if __name__ == "__main__":
    main()
