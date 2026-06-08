#!/usr/bin/env python3
"""Run the certified five-prime X+Y versus layer-compressed benchmark suite."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_N = 1_000_000_000_000
DEFAULT_PRIME_SETS = [
    "2,3,5,7,11",
    "2,3,5,7,13",
    "2,3,5,11,13",
    "2,3,7,11,13",
    "2,5,7,11,13",
    "3,5,7,11,13",
]


def load_single_harness():
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


def case_slug(primes: str) -> str:
    return "P" + primes.replace(",", "_")


def result_ok(result: dict[str, Any]) -> bool:
    return int(result.get("returncode", 1)) == 0


def run_case(single, primes: str, n: int, target_gap: int, heap_limit: int, skip_audit: bool) -> dict[str, Any]:
    runs = []
    xplusy = single.run_measured(
        "xplusy_adaptive",
        ["bin/smooth_xplusy_baseline", "bench", primes, str(n), str(heap_limit)],
    )
    runs.append(xplusy)

    layer = single.run_measured(
        "layer_compressed",
        ["bin/smooth_layer_compressed_general", "nth", primes, str(n), str(target_gap)],
    )
    runs.append(layer)

    layer_exps = exps_csv(layer.get("metrics", {}).get("exps"))
    if layer_exps and not skip_audit:
        audit = single.run_measured(
            "interval_audit",
            ["bin/smooth_interval_audit_exps_k6", primes, layer_exps, str(n)],
        )
        runs.append(audit)

    comparison = single.compute_comparison(runs)
    return {
        "primes": primes,
        "N": n,
        "slug": case_slug(primes),
        "runs": runs,
        "comparison": comparison,
    }


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in case["runs"]}
    xplusy = by_name.get("xplusy_adaptive", {})
    layer = by_name.get("layer_compressed", {})
    audit = by_name.get("interval_audit", {})
    comparison = case.get("comparison", {})
    return {
        "primes": case["primes"],
        "N": case["N"],
        "xplusy_returncode": xplusy.get("returncode"),
        "layer_returncode": layer.get("returncode"),
        "audit_returncode": audit.get("returncode", ""),
        "xplusy_wall_seconds": xplusy.get("elapsed_seconds", ""),
        "layer_wall_seconds": layer.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "xplusy_max_rss_kb": xplusy.get("max_rss_kb", ""),
        "layer_max_rss_kb": layer.get("max_rss_kb", ""),
        "audit_max_rss_kb": audit.get("max_rss_kb", ""),
        "xplusy_reported_total_over_layer_reported": comparison.get("xplusy_reported_total_over_layer_reported", ""),
        "xplusy_wall_over_layer_wall": comparison.get("xplusy_wall_over_layer_wall", ""),
        "layer_exps": layer.get("metrics", {}).get("exps", ""),
        "rank_certified": comparison.get("rank_certified", False),
        "certified_count_le": comparison.get("certified_count_le", ""),
        "layer_wins_wall": bool(comparison.get("xplusy_wall_over_layer_wall", 0) > 1.0),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["xplusy_returncode"] == 0 and row["layer_returncode"] == 0]
    certified = [row for row in completed if row["rank_certified"]]
    wins = [row for row in certified if row["layer_wins_wall"]]
    ratios = [float(row["xplusy_wall_over_layer_wall"]) for row in certified if row["xplusy_wall_over_layer_wall"] != ""]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "certified_cases": len(certified),
        "layer_wall_wins": len(wins),
        "all_completed": len(completed) == len(rows),
        "all_certified": len(certified) == len(rows),
        "layer_wall_wins_all_certified": len(wins) == len(rows),
        "min_wall_ratio": min(ratios) if ratios else None,
        "max_wall_ratio": max(ratios) if ratios else None,
        "mean_wall_ratio": sum(ratios) / len(ratios) if ratios else None,
    }


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    rows = report["summary_rows"]
    with (out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    lines = [
        "# Certified Five-Prime Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- N: `{report['config']['N']}`",
        f"- Cases: `{report['aggregate']['cases']}`",
        f"- Certified cases: `{report['aggregate']['certified_cases']}`",
        f"- Layer wall-time wins: `{report['aggregate']['layer_wall_wins']}`",
        f"- Mean wall ratio X+Y/layer: `{report['aggregate']['mean_wall_ratio']}`",
        "",
        "The `xplusy_adaptive` rows measure adaptive Cartesian-sum value selection only.",
        "The `layer_compressed` rows measure full exponent-vector unranking.",
        "",
        "| P | layer wall s | X+Y wall s | ratio | layer RSS KB | X+Y RSS KB | certified | exps |",
        "|---|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['primes']}` | {float(row['layer_wall_seconds']):.6f} | "
            f"{float(row['xplusy_wall_seconds']):.6f} | "
            f"{float(row['xplusy_wall_over_layer_wall']):.6f} | "
            f"{row['layer_max_rss_kb']} | {row['xplusy_max_rss_kb']} | "
            f"{row['rank_certified']} | `{row['layer_exps']}` |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--prime-set", action="append", dest="prime_sets")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--target-gap", type=int, default=50_000)
    parser.add_argument("--heap-limit", type=int, default=0)
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    single = load_single_harness()
    prime_sets = args.prime_sets or DEFAULT_PRIME_SETS
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"five_prime_suite_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.five_prime_suite.v1",
        "metadata": single.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "N": args.N,
            "prime_sets": prime_sets,
            "target_gap": args.target_gap,
            "heap_limit": args.heap_limit,
            "skip_audit": args.skip_audit,
        },
        "build_runs": [],
        "cases": [],
    }

    build_core = single.run_measured("build_core", ["./build.sh"])
    report["build_runs"].append(build_core)
    if not result_ok(build_core):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"build failed; wrote partial report to {out_dir}")
        return 1

    build_xy = single.build_baseline(args.cxx, args.cxxflags)
    report["build_runs"].append(build_xy)
    if not result_ok(build_xy):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "xplusy_build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"X+Y build failed; wrote partial report to {out_dir}")
        return 1

    for primes in prime_sets:
        report["cases"].append(run_case(single, primes, args.N, args.target_gap, args.heap_limit, args.skip_audit))

    rows = [summarize_case(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["all_completed"] and report["aggregate"]["all_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
