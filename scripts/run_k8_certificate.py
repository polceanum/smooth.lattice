#!/usr/bin/env python3
"""Run a reproducible k=8 sums-only unrank plus independent rank certificate."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_PRIMES = "2,3,5,7,11,13,17,19"
DEFAULT_N = 1_000_000_000_000


def load_harness():
    path = ROOT / "scripts" / "run_xplusy_vs_layer5.py"
    spec = importlib.util.spec_from_file_location("run_xplusy_vs_layer5", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def exps_csv(exps: Any) -> str | None:
    if isinstance(exps, list) and all(isinstance(item, int) for item in exps):
        return ",".join(str(item) for item in exps)
    return None


def ok(run: dict[str, Any]) -> bool:
    return int(run.get("returncode", 1)) == 0


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    row = report["summary"]
    with (out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(row.keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerow(row)

    lines = [
        "# Certified k=8 Sums-Only Unrank",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- P: `{report['config']['primes']}`",
        f"- N: `{report['config']['N']}`",
        f"- Solver return code: `{row['solver_returncode']}`",
        f"- Auditor return code: `{row['audit_returncode']}`",
        f"- Rank certified: `{row['rank_certified']}`",
        f"- Certified count_le: `{row['certified_count_le']}`",
        f"- Exponents: `{row['exps']}`",
        "",
        "| solver wall s | audit wall s | solver reported s | cert reported s | groupA | groupB | ambiguous | digits |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|",
        (
            f"| {float(row['solver_wall_seconds']):.6f} | "
            f"{float(row['audit_wall_seconds']):.6f} | "
            f"{float(row['solver_reported_seconds']):.6f} | "
            f"{float(row['cert_reported_seconds']):.6f} | "
            f"{row['groupA']} | {row['groupB']} | "
            f"{row['ambiguous_possible']} | {row['digits']} |"
        ),
        "",
        "The solver is the exploratory high-k sums-only MITM path. The certificate is",
        "the independent interval-log auditor with exact big-integer resolution of",
        "boundary ambiguities. This artifact supports correctness of this fixed",
        "instance; it is not a broad best-known claim.",
    ]
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", default=DEFAULT_PRIMES)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    args = parser.parse_args()

    harness = load_harness()
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"k8_certificate_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.k8_certificate.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {"primes": args.primes, "N": args.N},
        "runs": [],
    }

    build = harness.run_measured("build_core", ["./build.sh"])
    report["runs"].append(build)
    if not ok(build):
        report["summary"] = {"status": "build_failed"}
        write_outputs(out_dir, report)
        print(f"build failed; wrote partial report to {out_dir}")
        return 1

    solver = harness.run_measured(
        "sums_only_mitm",
        ["bin/smooth_sums_only_scalable", "nth", args.primes, str(args.N)],
    )
    report["runs"].append(solver)

    exps = exps_csv(solver.get("metrics", {}).get("exps"))
    if ok(solver) and exps is not None:
        audit = harness.run_measured(
            "interval_audit_exps",
            ["bin/smooth_interval_audit_exps", args.primes, exps, str(args.N)],
        )
    else:
        audit = {
            "name": "interval_audit_exps",
            "command": [],
            "returncode": 1,
            "elapsed_seconds": "",
            "stdout": "",
            "stderr": "solver did not return parseable exponents",
            "metrics": {},
        }
    report["runs"].append(audit)

    sm = solver.get("metrics", {})
    am = audit.get("metrics", {})
    report["summary"] = {
        "primes": args.primes,
        "N": args.N,
        "solver_returncode": solver.get("returncode"),
        "audit_returncode": audit.get("returncode"),
        "solver_wall_seconds": solver.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "solver_reported_seconds": sm.get("seconds", ""),
        "cert_reported_seconds": am.get("cert_seconds", ""),
        "rank_certified": bool(am.get("rank_certified", False)),
        "certified_count_le": am.get("certified_count_le", ""),
        "groupA": am.get("groupA", ""),
        "groupB": am.get("groupB", ""),
        "ambiguous_possible": am.get("ambiguous_possible", ""),
        "ambiguous_resolved": am.get("ambiguous_resolved", ""),
        "digits": sm.get("digits", ""),
        "exps": sm.get("exps", ""),
    }
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if ok(solver) and ok(audit) and report["summary"]["rank_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
