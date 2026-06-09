#!/usr/bin/env python3
"""Run a deterministic randomized certified suite for the corrected analytic oracle."""
from __future__ import annotations

import argparse
import csv
import importlib.util
import itertools
import json
import random
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)
DEFAULT_RANKS = (1_000_000_000, 1_000_000_000_000)
DEFAULT_K_VALUES = (5, 6, 8)
DEFAULT_CASES_PER_K_RANK = 2
DEFAULT_SEED = 20260609
RADIUS_BY_K = {5: 250, 6: 1500, 7: 5000, 8: 10000}
DEFAULT_MAX_CANDIDATES = 200_000


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


def exps_csv(value: Any) -> str | None:
    if isinstance(value, list):
        return ",".join(str(part) for part in value)
    if isinstance(value, str) and value:
        stripped = value.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            stripped = stripped[1:-1]
        return stripped.replace(" ", "")
    return None


def parse_csv_ints(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part]


def case_slug(primes: str, n: int, index: int) -> str:
    return f"N{n}_case{index}_P{primes.replace(',', '_')}"


def default_radius(k: int) -> int:
    if k in RADIUS_BY_K:
        return RADIUS_BY_K[k]
    return max(1000, 1000 * (k - 4))


def generate_cases(
    *,
    seed: int,
    ranks: list[int],
    k_values: list[int],
    cases_per_k_rank: int,
    max_candidates: int,
) -> list[dict[str, Any]]:
    cases: list[dict[str, Any]] = []
    for n in ranks:
        for k in k_values:
            combos = list(itertools.combinations(SUPPORTED_PRIMES, k))
            rng = random.Random(seed + 1009 * k + 9176 * len(str(n)) + n % 104729)
            rng.shuffle(combos)
            take = min(cases_per_k_rank, len(combos))
            for local_idx, combo in enumerate(combos[:take]):
                primes = ",".join(str(part) for part in combo)
                cases.append(
                    {
                        "label": f"k{k}_{case_slug(primes, n, local_idx)}",
                        "primes": primes,
                        "N": n,
                        "rank_radius": default_radius(k),
                        "max_candidates": max_candidates,
                    }
                )
    return cases


def case_from_arg(raw: str) -> dict[str, Any]:
    parts = raw.split(":")
    if len(parts) not in {3, 5}:
        raise argparse.ArgumentTypeError(
            "case must be label:primes_csv:N or label:primes_csv:N:rank_radius:max_candidates"
        )
    primes = parts[1]
    n = int(parts[2])
    k = len(parse_csv_ints(primes))
    return {
        "label": parts[0],
        "primes": primes,
        "N": n,
        "rank_radius": int(parts[3]) if len(parts) == 5 else default_radius(k),
        "max_candidates": int(parts[4]) if len(parts) == 5 else DEFAULT_MAX_CANDIDATES,
    }


def build_full_xplusy(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
    (ROOT / "bin").mkdir(exist_ok=True)
    return harness.run_measured(
        "build_xplusy_full_unrank",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/smooth_xplusy_full_unrank.cpp",
            "-o",
            "bin/smooth_xplusy_full_unrank",
        ],
    )


def run_case(harness, case: dict[str, Any], target_gap: int, skip_audit: bool) -> dict[str, Any]:
    label = case["label"]
    primes = case["primes"]
    n = int(case["N"])
    radius = int(case["rank_radius"])
    max_candidates = int(case["max_candidates"])
    runs: list[dict[str, Any]] = []

    corrected = harness.run_measured(
        "corrected_analytic_band",
        [
            "bin/smooth_layer_compressed_general",
            "analytic-band-corrected",
            primes,
            str(n),
            str(radius),
            str(max_candidates),
        ],
    )
    runs.append(corrected)

    corrected_exps = exps_csv(corrected.get("metrics", {}).get("exps", ""))
    if corrected_exps is not None and not skip_audit:
        audit = harness.run_measured(
            "corrected_interval_audit",
            ["bin/smooth_interval_audit_exps", primes, corrected_exps, str(n)],
        )
        runs.append(audit)

    sums = harness.run_measured(
        "sums_only_mitm",
        ["bin/smooth_sums_only_scalable", "nth", primes, str(n)],
    )
    runs.append(sums)

    xplusy = harness.run_measured(
        "xplusy_full_unrank",
        ["bin/smooth_xplusy_full_unrank", "nth", primes, str(n), str(target_gap)],
    )
    runs.append(xplusy)

    return {
        "label": label,
        "primes": primes,
        "N": n,
        "rank_radius": radius,
        "max_candidates": max_candidates,
        "runs": runs,
        "comparison": compare_runs(runs, n),
    }


def compare_runs(runs: list[dict[str, Any]], n: int) -> dict[str, Any]:
    by_name = {run["name"]: run for run in runs}
    corrected = by_name.get("corrected_analytic_band")
    audit = by_name.get("corrected_interval_audit")
    sums = by_name.get("sums_only_mitm")
    xplusy = by_name.get("xplusy_full_unrank")
    cm = (corrected or {}).get("metrics", {})
    am = (audit or {}).get("metrics", {})
    sm = (sums or {}).get("metrics", {})
    xm = (xplusy or {}).get("metrics", {})
    corrected_exps = cm.get("exps", "")
    comparison: dict[str, Any] = {
        "corrected_certified": bool(am.get("rank_certified", False)) and am.get("certified_count_le", "") == n,
        "sums_same_exps": sm.get("exps", None) == corrected_exps,
        "xplusy_same_exps": xm.get("exps", None) == corrected_exps,
    }
    if ok(corrected) and ok(sums):
        cw = float(corrected.get("elapsed_seconds", 0.0))
        sw = float(sums.get("elapsed_seconds", 0.0))
        if cw > 0:
            comparison["sums_wall_over_corrected_wall"] = sw / cw
    if ok(corrected) and ok(xplusy):
        cw = float(corrected.get("elapsed_seconds", 0.0))
        xw = float(xplusy.get("elapsed_seconds", 0.0))
        if cw > 0:
            comparison["xplusy_wall_over_corrected_wall"] = xw / cw
    if corrected and sums and corrected.get("max_rss_kb", "") and sums.get("max_rss_kb", ""):
        cr = float(corrected["max_rss_kb"])
        sr = float(sums["max_rss_kb"])
        if cr > 0:
            comparison["sums_rss_over_corrected_rss"] = sr / cr
    if corrected and xplusy and corrected.get("max_rss_kb", "") and xplusy.get("max_rss_kb", ""):
        cr = float(corrected["max_rss_kb"])
        xr = float(xplusy["max_rss_kb"])
        if cr > 0:
            comparison["xplusy_rss_over_corrected_rss"] = xr / cr
    return comparison


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in case["runs"]}
    corrected = by_name.get("corrected_analytic_band", {})
    audit = by_name.get("corrected_interval_audit", {})
    sums = by_name.get("sums_only_mitm", {})
    xplusy = by_name.get("xplusy_full_unrank", {})
    cm = corrected.get("metrics", {})
    am = audit.get("metrics", {})
    sm = sums.get("metrics", {})
    xm = xplusy.get("metrics", {})
    comparison = case["comparison"]
    return {
        "label": case["label"],
        "primes": case["primes"],
        "k": len(parse_csv_ints(case["primes"])),
        "N": case["N"],
        "rank_radius": case["rank_radius"],
        "max_candidates": case["max_candidates"],
        "corrected_returncode": corrected.get("returncode", ""),
        "audit_returncode": audit.get("returncode", ""),
        "sums_returncode": sums.get("returncode", ""),
        "xplusy_returncode": xplusy.get("returncode", ""),
        "corrected_wall_seconds": corrected.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "sums_wall_seconds": sums.get("elapsed_seconds", ""),
        "xplusy_wall_seconds": xplusy.get("elapsed_seconds", ""),
        "corrected_reported_seconds": cm.get("seconds", ""),
        "sums_reported_seconds": sm.get("seconds", ""),
        "xplusy_reported_seconds": xm.get("seconds", ""),
        "corrected_max_rss_kb": corrected.get("max_rss_kb", ""),
        "sums_max_rss_kb": sums.get("max_rss_kb", ""),
        "xplusy_max_rss_kb": xplusy.get("max_rss_kb", ""),
        "center_rank_error": int(cm.get("center_count", 0)) - int(case["N"]) if cm.get("center_count", "") != "" else "",
        "band_count": cm.get("band_count", ""),
        "target_inside": bool(cm.get("target_inside", False)),
        "enumerated": bool(cm.get("enumerated", False)),
        "recovered": bool(cm.get("recovered", False)),
        "rank_certified": bool(am.get("rank_certified", False)),
        "certified_count_le": am.get("certified_count_le", ""),
        "cert_seconds": am.get("cert_seconds", ""),
        "ambiguous_possible": am.get("ambiguous_possible", ""),
        "ambiguous_resolved": am.get("ambiguous_resolved", ""),
        "sums_same_exps": comparison.get("sums_same_exps", False),
        "xplusy_same_exps": comparison.get("xplusy_same_exps", False),
        "sums_wall_over_corrected_wall": comparison.get("sums_wall_over_corrected_wall", ""),
        "xplusy_wall_over_corrected_wall": comparison.get("xplusy_wall_over_corrected_wall", ""),
        "sums_rss_over_corrected_rss": comparison.get("sums_rss_over_corrected_rss", ""),
        "xplusy_rss_over_corrected_rss": comparison.get("xplusy_rss_over_corrected_rss", ""),
        "exps": cm.get("exps", ""),
        "digits": cm.get("digits", ""),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [
        row for row in rows
        if row["corrected_returncode"] == 0 and row["sums_returncode"] == 0 and row["xplusy_returncode"] == 0
    ]
    certified = [
        row for row in completed
        if row["rank_certified"] and row["certified_count_le"] == row["N"]
    ]
    all_match = [row for row in certified if row["sums_same_exps"] and row["xplusy_same_exps"]]
    corrected_vs_sums = [
        float(row["sums_wall_over_corrected_wall"])
        for row in all_match
        if row["sums_wall_over_corrected_wall"] != ""
    ]
    corrected_vs_xplusy = [
        float(row["xplusy_wall_over_corrected_wall"])
        for row in all_match
        if row["xplusy_wall_over_corrected_wall"] != ""
    ]
    memory_vs_xplusy = [
        float(row["xplusy_rss_over_corrected_rss"])
        for row in all_match
        if row["xplusy_rss_over_corrected_rss"] != ""
    ]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "certified_cases": len(certified),
        "all_methods_same_exps_cases": len(all_match),
        "corrected_wall_wins_vs_sums": sum(1 for ratio in corrected_vs_sums if ratio > 1.0),
        "corrected_wall_wins_vs_xplusy": sum(1 for ratio in corrected_vs_xplusy if ratio > 1.0),
        "corrected_memory_wins_vs_xplusy": sum(1 for ratio in memory_vs_xplusy if ratio > 1.0),
        "mean_sums_wall_over_corrected": sum(corrected_vs_sums) / len(corrected_vs_sums) if corrected_vs_sums else None,
        "mean_xplusy_wall_over_corrected": sum(corrected_vs_xplusy) / len(corrected_vs_xplusy) if corrected_vs_xplusy else None,
        "mean_xplusy_rss_over_corrected": sum(memory_vs_xplusy) / len(memory_vs_xplusy) if memory_vs_xplusy else None,
        "max_band_count": max([int(row["band_count"]) for row in completed if row["band_count"] != ""] or [0]),
    }


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    rows = report["summary_rows"]
    with (out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    agg = report["aggregate"]
    lines = [
        "# Certified Corrected-Oracle Random Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Seed: `{report['config']['seed']}`",
        f"- Cases: `{agg['cases']}`",
        f"- Completed cases: `{agg['completed_cases']}`",
        f"- Certified corrected-oracle cases: `{agg['certified_cases']}`",
        f"- All-method same-exponent cases: `{agg['all_methods_same_exps_cases']}`",
        f"- Corrected wall-time wins vs sums-only: `{agg['corrected_wall_wins_vs_sums']}`",
        f"- Corrected wall-time wins vs full X+Y: `{agg['corrected_wall_wins_vs_xplusy']}`",
        f"- Corrected memory wins vs full X+Y: `{agg['corrected_memory_wins_vs_xplusy']}`",
        f"- Mean wall ratio sums/corrected: `{agg['mean_sums_wall_over_corrected']}`",
        f"- Mean wall ratio full X+Y/corrected: `{agg['mean_xplusy_wall_over_corrected']}`",
        f"- Mean RSS ratio full X+Y/corrected: `{agg['mean_xplusy_rss_over_corrected']}`",
        f"- Max corrected band count: `{agg['max_band_count']}`",
        "",
        "Cases are deterministic pseudo-random subsets of the supported audit prime",
        "universe. The corrected analytic oracle is independently interval-audited;",
        "the sums-only and full materialized X+Y baselines must return the same",
        "exponent vector to count as matching the certified result.",
        "",
        "| label | k | N | band | cert | sums same | X+Y same | corr wall | sums wall | X+Y wall | X+Y/corr RSS | exps |",
        "|---|---:|---:|---:|---|---|---|---:|---:|---:|---:|---|",
    ]
    for row in rows:
        rss_ratio = row["xplusy_rss_over_corrected_rss"]
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['N']} | {row['band_count']} | "
            f"{row['rank_certified']} | {row['sums_same_exps']} | {row['xplusy_same_exps']} | "
            f"{float(row['corrected_wall_seconds']):.6f} | {float(row['sums_wall_seconds']):.6f} | "
            f"{float(row['xplusy_wall_seconds']):.6f} | "
            f"{float(rss_ratio):.3f} | `{row['exps']}` |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--rank", action="append", type=int, dest="ranks")
    parser.add_argument("--k", action="append", type=int, dest="k_values")
    parser.add_argument("--cases-per-k-rank", type=int, default=DEFAULT_CASES_PER_K_RANK)
    parser.add_argument("--case", action="append", type=case_from_arg, dest="cases")
    parser.add_argument("--max-candidates", type=int, default=DEFAULT_MAX_CANDIDATES)
    parser.add_argument("--target-gap", type=int, default=50_000)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    harness = load_harness()
    ranks = args.ranks or list(DEFAULT_RANKS)
    k_values = args.k_values or list(DEFAULT_K_VALUES)
    cases = args.cases or generate_cases(
        seed=args.seed,
        ranks=ranks,
        k_values=k_values,
        cases_per_k_rank=args.cases_per_k_rank,
        max_candidates=args.max_candidates,
    )
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"corrected_oracle_random_suite_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.corrected_oracle_random_suite.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "seed": args.seed,
            "ranks": ranks,
            "k_values": k_values,
            "cases_per_k_rank": args.cases_per_k_rank,
            "cases": cases,
            "target_gap": args.target_gap,
            "max_candidates": args.max_candidates,
            "skip_audit": args.skip_audit,
        },
        "build_runs": [],
        "cases": [],
    }

    build_core = harness.run_measured("build_core", ["./build.sh"])
    report["build_runs"].append(build_core)
    if not ok(build_core):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"build failed; wrote partial report to {out_dir}")
        return 1

    build_xy = build_full_xplusy(harness, args.cxx, args.cxxflags)
    report["build_runs"].append(build_xy)
    if not ok(build_xy):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "xplusy_full_build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"full X+Y build failed; wrote partial report to {out_dir}")
        return 1

    for case in cases:
        report["cases"].append(run_case(harness, case, args.target_gap, args.skip_audit))

    rows = [summarize_case(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["all_methods_same_exps_cases"] == report["aggregate"]["cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
