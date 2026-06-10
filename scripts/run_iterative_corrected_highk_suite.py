#!/usr/bin/env python3
"""Run certified high-k iterative residual-correction benchmarks.

The suite compares the adaptive sums-only MITM unranker with the iterative
residual-corrected analytic-band mode on higher-dimensional fixed-prime cases.
Each corrected result is independently interval-audited.
"""
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
DEFAULT_CASES = [
    {
        "label": "k8_first8",
        "primes": "2,3,5,7,11,13,17,19",
        "rank_radius": 25,
        "max_candidates": 200_000,
        "refine_steps": 4,
    },
    {
        "label": "k10_first10",
        "primes": "2,3,5,7,11,13,17,19,23,29",
        "rank_radius": 25,
        "max_candidates": 200_000,
        "refine_steps": 4,
    },
    {
        "label": "k12_first12",
        "primes": "2,3,5,7,11,13,17,19,23,29,31,37",
        "rank_radius": 25,
        "max_candidates": 200_000,
        "refine_steps": 4,
    },
]


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


def ok(run: dict[str, Any] | None) -> bool:
    return run is not None and int(run.get("returncode", 1)) == 0


def exps_csv(exps: Any) -> str | None:
    if isinstance(exps, list) and all(isinstance(item, int) for item in exps):
        return ",".join(str(item) for item in exps)
    if isinstance(exps, str) and exps:
        stripped = exps.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            stripped = stripped[1:-1]
        return stripped.replace(" ", "")
    return None


def build_sums(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
    harness.BIN.mkdir(exist_ok=True)
    return harness.run_measured(
        "build_sums_only_scalable",
        [
            cxx,
            *cxxflags.split(),
            "src/smooth_sums_only_scalable.cpp",
            "-o",
            "bin/smooth_sums_only_scalable",
        ],
    )


def build_auditor(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
    harness.BIN.mkdir(exist_ok=True)
    return harness.run_measured(
        "build_interval_auditor",
        [
            cxx,
            *cxxflags.split(),
            "src/smooth_interval_audit_exps_k6.cpp",
            "-o",
            "bin/smooth_interval_audit_exps",
        ],
    )


def run_case(harness, case: dict[str, Any], n: int, skip_audit: bool) -> dict[str, Any]:
    primes = case["primes"]
    adaptive = harness.run_measured(
        "sums_adaptive",
        ["bin/smooth_sums_only_scalable", "nth", primes, str(n)],
    )
    corrected = harness.run_measured(
        "sums_iterative_corrected",
        [
            "bin/smooth_sums_only_scalable",
            "analytic-band-corrected",
            primes,
            str(n),
            str(case["rank_radius"]),
            str(case["max_candidates"]),
            str(case["refine_steps"]),
        ],
    )
    runs = [adaptive, corrected]
    audit = None
    corrected_exps = corrected.get("metrics", {}).get("exps")
    if not skip_audit and ok(corrected):
        exps = exps_csv(corrected_exps)
        if exps is not None:
            audit = harness.run_measured(
                "interval_audit",
                ["bin/smooth_interval_audit_exps", primes, exps, str(n)],
            )
            runs.append(audit)
    return {
        **case,
        "N": n,
        "runs": runs,
        "comparison": compare(adaptive, corrected, audit, n),
    }


def compare(
    adaptive: dict[str, Any],
    corrected: dict[str, Any],
    audit: dict[str, Any] | None,
    n: int,
) -> dict[str, Any]:
    same_exps = (
        ok(adaptive)
        and ok(corrected)
        and adaptive.get("metrics", {}).get("exps") == corrected.get("metrics", {}).get("exps")
    )
    wall_ratio = ""
    reported_ratio = ""
    if ok(adaptive) and ok(corrected):
        corrected_wall = float(corrected["elapsed_seconds"])
        if corrected_wall > 0:
            wall_ratio = float(adaptive["elapsed_seconds"]) / corrected_wall
        corrected_sec = corrected.get("metrics", {}).get("seconds", "")
        adaptive_sec = adaptive.get("metrics", {}).get("seconds", "")
        if corrected_sec != "" and adaptive_sec != "" and float(corrected_sec) > 0:
            reported_ratio = float(adaptive_sec) / float(corrected_sec)
    certified = bool(
        audit
        and ok(audit)
        and audit.get("metrics", {}).get("rank_certified", False)
        and audit.get("metrics", {}).get("certified_count_le", "") == n
    )
    return {
        "same_exps": same_exps,
        "certified": certified,
        "adaptive_wall_over_corrected_wall": wall_ratio,
        "adaptive_reported_over_corrected_reported": reported_ratio,
    }


def summarize(case: dict[str, Any]) -> dict[str, Any]:
    runs = {run["name"]: run for run in case["runs"]}
    adaptive = runs.get("sums_adaptive", {})
    corrected = runs.get("sums_iterative_corrected", {})
    audit = runs.get("interval_audit", {})
    comparison = case["comparison"]
    return {
        "label": case["label"],
        "primes": case["primes"],
        "k": len([p for p in case["primes"].split(",") if p]),
        "N": case["N"],
        "rank_radius": case["rank_radius"],
        "refine_steps": case["refine_steps"],
        "adaptive_returncode": adaptive.get("returncode", ""),
        "corrected_returncode": corrected.get("returncode", ""),
        "audit_returncode": audit.get("returncode", ""),
        "same_exps": comparison["same_exps"],
        "certified": comparison["certified"],
        "adaptive_wall_seconds": adaptive.get("elapsed_seconds", ""),
        "corrected_wall_seconds": corrected.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "adaptive_reported_seconds": adaptive.get("metrics", {}).get("seconds", ""),
        "corrected_reported_seconds": corrected.get("metrics", {}).get("seconds", ""),
        "corrected_build_seconds": corrected.get("metrics", {}).get("build", ""),
        "corrected_count_seconds": corrected.get("metrics", {}).get("count_phase", ""),
        "corrected_band_seconds": corrected.get("metrics", {}).get("band_phase", ""),
        "adaptive_max_rss_kb": adaptive.get("max_rss_kb", ""),
        "corrected_max_rss_kb": corrected.get("max_rss_kb", ""),
        "audit_max_rss_kb": audit.get("max_rss_kb", ""),
        "adaptive_wall_over_corrected_wall": comparison["adaptive_wall_over_corrected_wall"],
        "adaptive_reported_over_corrected_reported": comparison["adaptive_reported_over_corrected_reported"],
        "center_rank_error": corrected.get("metrics", {}).get("center_rank_error", ""),
        "final_rank_error": corrected.get("metrics", {}).get("final_rank_error", ""),
        "band_count": corrected.get("metrics", {}).get("band_count", ""),
        "cands": corrected.get("metrics", {}).get("cands", ""),
        "target_inside": corrected.get("metrics", {}).get("target_inside", ""),
        "recovered": corrected.get("metrics", {}).get("recovered", ""),
        "cert_seconds": audit.get("metrics", {}).get("cert_seconds", ""),
        "ambiguous_possible": audit.get("metrics", {}).get("ambiguous_possible", ""),
        "ambiguous_resolved": audit.get("metrics", {}).get("ambiguous_resolved", ""),
        "exps": corrected.get("metrics", {}).get("exps", ""),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["adaptive_returncode"] == 0 and row["corrected_returncode"] == 0]
    certified = [row for row in completed if row["same_exps"] and row["certified"]]
    ratios = [
        float(row["adaptive_wall_over_corrected_wall"])
        for row in certified
        if row["adaptive_wall_over_corrected_wall"] != ""
    ]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "same_exps_certified_cases": len(certified),
        "corrected_wall_wins": sum(1 for row in certified if float(row["adaptive_wall_over_corrected_wall"]) > 1.0),
        "all_same_exps_certified": len(certified) == len(rows),
        "min_adaptive_wall_over_corrected_wall": min(ratios) if ratios else None,
        "max_adaptive_wall_over_corrected_wall": max(ratios) if ratios else None,
        "mean_adaptive_wall_over_corrected_wall": sum(ratios) / len(ratios) if ratios else None,
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
        "# Iterative Corrected High-k Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- N: `{report['config']['N']}`",
        f"- Cases: `{report['aggregate']['cases']}`",
        f"- Same-exponent certified cases: `{report['aggregate']['same_exps_certified_cases']}`",
        f"- Corrected wall-time wins: `{report['aggregate']['corrected_wall_wins']}`",
        f"- Mean adaptive/corrected wall ratio: `{report['aggregate']['mean_adaptive_wall_over_corrected_wall']}`",
        "",
        "The corrected row applies a shared iterative policy: four exact residual",
        "corrections and a rank-radius-25 final exact-sorted band.",
        "Every corrected output is independently interval-rank certified.",
        "",
        "| label | k | radius | refine | adaptive wall s | corrected wall s | ratio | band | cands | audit s | certified | exps |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['rank_radius']} | {row['refine_steps']} | "
            f"{float(row['adaptive_wall_seconds']):.6f} | {float(row['corrected_wall_seconds']):.6f} | "
            f"{float(row['adaptive_wall_over_corrected_wall']):.6f} | {row['band_count']} | "
            f"{row['cands']} | {float(row['audit_wall_seconds']):.6f} | {row['certified']} | "
            f"`{row['exps']}` |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def parse_case(raw: str) -> dict[str, Any]:
    parts = raw.split(":")
    if len(parts) != 5:
        raise argparse.ArgumentTypeError("case must be label:primes_csv:rank_radius:max_candidates:refine_steps")
    return {
        "label": parts[0],
        "primes": parts[1],
        "rank_radius": int(parts[2]),
        "max_candidates": int(parts[3]),
        "refine_steps": int(parts[4]),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--case", action="append", type=parse_case)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    harness = load_harness()
    cases = args.case or DEFAULT_CASES
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"iterative_corrected_highk_{utc_stamp()}")
    build_runs = [
        build_sums(harness, args.cxx, args.cxxflags),
        build_auditor(harness, args.cxx, args.cxxflags),
    ]
    case_reports = [run_case(harness, case, args.N, args.skip_audit) for case in cases]
    rows = [summarize(case) for case in case_reports]
    report = {
        "schema": "smooth.lattice.iterative_corrected_highk.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "N": args.N,
            "cases": cases,
            "skip_audit": args.skip_audit,
        },
        "build_runs": build_runs,
        "cases": case_reports,
        "summary_rows": rows,
        "aggregate": aggregate(rows),
    }
    write_outputs(out_dir, report)
    print(f"wrote {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
