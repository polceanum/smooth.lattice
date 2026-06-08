#!/usr/bin/env python3
"""Run the certified five-prime suite against full materialized X+Y unranking."""
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


def build_full_xplusy(single, cxx: str, cxxflags: str) -> dict[str, Any]:
    single.BIN.mkdir(exist_ok=True)
    return single.run_measured(
        "build_xplusy_full_unrank",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/smooth_xplusy_full_unrank.cpp",
            "-o",
            "bin/smooth_xplusy_full_unrank",
        ],
    )


def audit_exps(single, name: str, primes: str, exps: Any, n: int) -> dict[str, Any] | None:
    exps_text = exps_csv(exps)
    if exps_text is None:
        return None
    return single.run_measured(
        name,
        ["bin/smooth_interval_audit_exps_k6", primes, exps_text, str(n)],
    )


def run_case(single, primes: str, n: int, target_gap: int, skip_audit: bool) -> dict[str, Any]:
    runs = []
    xplusy = single.run_measured(
        "xplusy_full_unrank",
        ["bin/smooth_xplusy_full_unrank", "nth", primes, str(n), str(target_gap)],
    )
    runs.append(xplusy)

    layer = single.run_measured(
        "layer_compressed",
        ["bin/smooth_layer_compressed_general", "nth", primes, str(n), str(target_gap)],
    )
    runs.append(layer)

    if not skip_audit:
        layer_audit = audit_exps(single, "layer_interval_audit", primes, layer.get("metrics", {}).get("exps"), n)
        if layer_audit is not None:
            runs.append(layer_audit)
        xplusy_audit = audit_exps(single, "xplusy_interval_audit", primes, xplusy.get("metrics", {}).get("exps"), n)
        if xplusy_audit is not None:
            runs.append(xplusy_audit)

    return {
        "primes": primes,
        "N": n,
        "slug": case_slug(primes),
        "runs": runs,
        "comparison": compute_comparison(runs),
    }


def compute_comparison(runs: list[dict[str, Any]]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in runs}
    xplusy = by_name.get("xplusy_full_unrank", {})
    layer = by_name.get("layer_compressed", {})
    layer_audit = by_name.get("layer_interval_audit", {})
    xplusy_audit = by_name.get("xplusy_interval_audit", {})
    comparison: dict[str, Any] = {}
    if result_ok(xplusy) and result_ok(layer):
        xplusy_wall = float(xplusy["elapsed_seconds"])
        layer_wall = float(layer["elapsed_seconds"])
        if layer_wall > 0:
            comparison["xplusy_wall_over_layer_wall"] = xplusy_wall / layer_wall
        layer_sec = layer["metrics"].get("seconds")
        xplusy_sec = xplusy["metrics"].get("seconds")
        if layer_sec and xplusy_sec:
            comparison["xplusy_reported_over_layer_reported"] = float(xplusy_sec) / float(layer_sec)
        comparison["same_exps"] = xplusy["metrics"].get("exps") == layer["metrics"].get("exps")
        comparison["same_rank"] = xplusy["metrics"].get("N") == layer["metrics"].get("N")
        comparison["same_primes"] = xplusy["metrics"].get("P") == layer["metrics"].get("P")
    comparison["layer_rank_certified"] = bool(layer_audit.get("metrics", {}).get("rank_certified", False))
    comparison["xplusy_rank_certified"] = bool(xplusy_audit.get("metrics", {}).get("rank_certified", False))
    comparison["layer_certified_count_le"] = layer_audit.get("metrics", {}).get("certified_count_le", "")
    comparison["xplusy_certified_count_le"] = xplusy_audit.get("metrics", {}).get("certified_count_le", "")
    return comparison


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in case["runs"]}
    xplusy = by_name.get("xplusy_full_unrank", {})
    layer = by_name.get("layer_compressed", {})
    layer_audit = by_name.get("layer_interval_audit", {})
    xplusy_audit = by_name.get("xplusy_interval_audit", {})
    comparison = case.get("comparison", {})
    return {
        "primes": case["primes"],
        "N": case["N"],
        "xplusy_returncode": xplusy.get("returncode"),
        "layer_returncode": layer.get("returncode"),
        "layer_audit_returncode": layer_audit.get("returncode", ""),
        "xplusy_audit_returncode": xplusy_audit.get("returncode", ""),
        "xplusy_wall_seconds": xplusy.get("elapsed_seconds", ""),
        "layer_wall_seconds": layer.get("elapsed_seconds", ""),
        "layer_audit_wall_seconds": layer_audit.get("elapsed_seconds", ""),
        "xplusy_audit_wall_seconds": xplusy_audit.get("elapsed_seconds", ""),
        "xplusy_max_rss_kb": xplusy.get("max_rss_kb", ""),
        "layer_max_rss_kb": layer.get("max_rss_kb", ""),
        "xplusy_reported_seconds": xplusy.get("metrics", {}).get("seconds", ""),
        "layer_reported_seconds": layer.get("metrics", {}).get("seconds", ""),
        "xplusy_A": xplusy.get("metrics", {}).get("A", ""),
        "xplusy_B": xplusy.get("metrics", {}).get("B", ""),
        "xplusy_band": xplusy.get("metrics", {}).get("band", ""),
        "layer_A": layer.get("metrics", {}).get("A", ""),
        "layer_Base": layer.get("metrics", {}).get("Base", ""),
        "layer_band": layer.get("metrics", {}).get("band", ""),
        "xplusy_wall_over_layer_wall": comparison.get("xplusy_wall_over_layer_wall", ""),
        "xplusy_reported_over_layer_reported": comparison.get("xplusy_reported_over_layer_reported", ""),
        "same_exps": comparison.get("same_exps", False),
        "layer_exps": layer.get("metrics", {}).get("exps", ""),
        "xplusy_exps": xplusy.get("metrics", {}).get("exps", ""),
        "layer_rank_certified": comparison.get("layer_rank_certified", False),
        "xplusy_rank_certified": comparison.get("xplusy_rank_certified", False),
        "layer_wins_wall": bool(comparison.get("xplusy_wall_over_layer_wall", 0) > 1.0),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["xplusy_returncode"] == 0 and row["layer_returncode"] == 0]
    certified = [
        row for row in completed
        if row["layer_rank_certified"] and row["xplusy_rank_certified"] and row["same_exps"]
    ]
    wins = [row for row in certified if row["layer_wins_wall"]]
    ratios = [float(row["xplusy_wall_over_layer_wall"]) for row in certified if row["xplusy_wall_over_layer_wall"] != ""]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "same_exps_certified_cases": len(certified),
        "layer_wall_wins": len(wins),
        "all_completed": len(completed) == len(rows),
        "all_same_exps_certified": len(certified) == len(rows),
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
        "# Certified Full X+Y Five-Prime Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- N: `{report['config']['N']}`",
        f"- Cases: `{report['aggregate']['cases']}`",
        f"- Same-exponent certified cases: `{report['aggregate']['same_exps_certified_cases']}`",
        f"- Layer wall-time wins: `{report['aggregate']['layer_wall_wins']}`",
        f"- Mean wall ratio full X+Y/layer: `{report['aggregate']['mean_wall_ratio']}`",
        "",
        "Both compared rows return full exponent vectors. Each successful case audits both outputs",
        "with the independent interval-rank auditor and checks that the exponent vectors match.",
        "",
        "| P | layer wall s | full X+Y wall s | ratio | layer RSS KB | X+Y RSS KB | same exps | both certified | exps |",
        "|---|---:|---:|---:|---:|---:|---|---|---|",
    ]
    for row in rows:
        both_certified = row["layer_rank_certified"] and row["xplusy_rank_certified"]
        ratio = row["xplusy_wall_over_layer_wall"]
        lines.append(
            f"| `{row['primes']}` | {float(row['layer_wall_seconds']):.6f} | "
            f"{float(row['xplusy_wall_seconds']):.6f} | "
            f"{float(ratio):.6f} | "
            f"{row['layer_max_rss_kb']} | {row['xplusy_max_rss_kb']} | "
            f"{row['same_exps']} | {both_certified} | `{row['layer_exps']}` |"
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
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    single = load_single_harness()
    prime_sets = args.prime_sets or DEFAULT_PRIME_SETS
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"full_xplusy_suite_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.full_xplusy_suite.v1",
        "metadata": single.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "N": args.N,
            "prime_sets": prime_sets,
            "target_gap": args.target_gap,
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

    build_xy = build_full_xplusy(single, args.cxx, args.cxxflags)
    report["build_runs"].append(build_xy)
    if not result_ok(build_xy):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "xplusy_full_build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"full X+Y build failed; wrote partial report to {out_dir}")
        return 1

    for primes in prime_sets:
        report["cases"].append(run_case(single, primes, args.N, args.target_gap, args.skip_audit))

    rows = [summarize_case(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    ok = report["aggregate"]["all_completed"] and report["aggregate"]["all_same_exps_certified"]
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
