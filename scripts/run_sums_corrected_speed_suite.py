#!/usr/bin/env python3
"""Compare residual-corrected sums-only MITM against adaptive sums-only MITM."""
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
DEFAULT_MAX_CANDIDATES = 200_000
RADIUS_BY_K = {5: 250, 6: 1500, 7: 5000, 8: 10000}


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


def parse_csv_ints(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part]


def exps_csv(value: Any) -> str | None:
    if isinstance(value, list):
        return ",".join(str(part) for part in value)
    if isinstance(value, str) and value:
        stripped = value.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            stripped = stripped[1:-1]
        return stripped.replace(" ", "")
    return None


def default_radius(k: int) -> int:
    return RADIUS_BY_K.get(k, max(1000, 1000 * (k - 4)))


def case_slug(primes: str, n: int, index: int) -> str:
    return f"N{n}_case{index}_P{primes.replace(',', '_')}"


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
            for local_idx, combo in enumerate(combos[: min(cases_per_k_rank, len(combos))]):
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


def run_case(harness, case: dict[str, Any], skip_audit: bool) -> dict[str, Any]:
    primes = case["primes"]
    n = int(case["N"])
    radius = int(case["rank_radius"])
    max_candidates = int(case["max_candidates"])
    runs: list[dict[str, Any]] = []

    corrected = harness.run_measured(
        "sums_corrected_band",
        [
            "bin/smooth_sums_only_scalable",
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
            "sums_corrected_interval_audit",
            ["bin/smooth_interval_audit_exps", primes, corrected_exps, str(n)],
        )
        runs.append(audit)

    adaptive = harness.run_measured(
        "sums_only_mitm",
        ["bin/smooth_sums_only_scalable", "nth", primes, str(n)],
    )
    runs.append(adaptive)

    return {
        "label": case["label"],
        "primes": primes,
        "N": n,
        "rank_radius": radius,
        "max_candidates": max_candidates,
        "runs": runs,
        "comparison": compare_runs(runs, n),
    }


def compare_runs(runs: list[dict[str, Any]], n: int) -> dict[str, Any]:
    by_name = {run["name"]: run for run in runs}
    corrected = by_name.get("sums_corrected_band")
    audit = by_name.get("sums_corrected_interval_audit")
    adaptive = by_name.get("sums_only_mitm")
    cm = (corrected or {}).get("metrics", {})
    am = (audit or {}).get("metrics", {})
    sm = (adaptive or {}).get("metrics", {})
    comparison: dict[str, Any] = {
        "corrected_certified": bool(am.get("rank_certified", False)) and am.get("certified_count_le", "") == n,
        "adaptive_same_exps": sm.get("exps", None) == cm.get("exps", ""),
    }
    if ok(corrected) and ok(adaptive):
        cw = float(corrected.get("elapsed_seconds", 0.0))
        aw = float(adaptive.get("elapsed_seconds", 0.0))
        if cw > 0:
            comparison["adaptive_wall_over_corrected_wall"] = aw / cw
        cr = corrected.get("max_rss_kb", "")
        ar = adaptive.get("max_rss_kb", "")
        if cr != "" and ar != "" and float(cr) > 0:
            comparison["adaptive_rss_over_corrected_rss"] = float(ar) / float(cr)
        csec = cm.get("seconds", "")
        asec = sm.get("seconds", "")
        if csec != "" and asec != "" and float(csec) > 0:
            comparison["adaptive_reported_over_corrected_reported"] = float(asec) / float(csec)
    return comparison


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in case["runs"]}
    corrected = by_name.get("sums_corrected_band", {})
    audit = by_name.get("sums_corrected_interval_audit", {})
    adaptive = by_name.get("sums_only_mitm", {})
    cm = corrected.get("metrics", {})
    am = audit.get("metrics", {})
    sm = adaptive.get("metrics", {})
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
        "adaptive_returncode": adaptive.get("returncode", ""),
        "corrected_wall_seconds": corrected.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "adaptive_wall_seconds": adaptive.get("elapsed_seconds", ""),
        "corrected_reported_seconds": cm.get("seconds", ""),
        "adaptive_reported_seconds": sm.get("seconds", ""),
        "corrected_max_rss_kb": corrected.get("max_rss_kb", ""),
        "adaptive_max_rss_kb": adaptive.get("max_rss_kb", ""),
        "corrected_build_seconds": cm.get("build", ""),
        "adaptive_build_seconds": sm.get("build", ""),
        "corrected_count_seconds": cm.get("count_phase", ""),
        "adaptive_count_seconds": sm.get("count", ""),
        "corrected_band_seconds": cm.get("band_phase", ""),
        "adaptive_band_seconds": sm.get("band_phase", ""),
        "adaptive_calls": sm.get("calls", ""),
        "adaptive_rank_gap": sm.get("gap", ""),
        "center_rank_error": cm.get("center_rank_error", ""),
        "correction": cm.get("correction", ""),
        "half_width": cm.get("half_width", ""),
        "corrected_band_count": cm.get("band_count", ""),
        "adaptive_band_count": sm.get("band", ""),
        "target_inside": bool(cm.get("target_inside", False)),
        "enumerated": bool(cm.get("enumerated", False)),
        "recovered": bool(cm.get("recovered", False)),
        "rank_certified": bool(am.get("rank_certified", False)),
        "certified_count_le": am.get("certified_count_le", ""),
        "cert_seconds": am.get("cert_seconds", ""),
        "adaptive_same_exps": comparison.get("adaptive_same_exps", False),
        "adaptive_wall_over_corrected_wall": comparison.get("adaptive_wall_over_corrected_wall", ""),
        "adaptive_reported_over_corrected_reported": comparison.get("adaptive_reported_over_corrected_reported", ""),
        "adaptive_rss_over_corrected_rss": comparison.get("adaptive_rss_over_corrected_rss", ""),
        "exps": cm.get("exps", ""),
        "digits": cm.get("digits", ""),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["corrected_returncode"] == 0 and row["adaptive_returncode"] == 0]
    certified = [
        row for row in completed
        if row["rank_certified"] and row["certified_count_le"] == row["N"] and row["adaptive_same_exps"]
    ]
    wall = [float(row["adaptive_wall_over_corrected_wall"]) for row in certified if row["adaptive_wall_over_corrected_wall"] != ""]
    reported = [
        float(row["adaptive_reported_over_corrected_reported"])
        for row in certified
        if row["adaptive_reported_over_corrected_reported"] != ""
    ]
    rss = [float(row["adaptive_rss_over_corrected_rss"]) for row in certified if row["adaptive_rss_over_corrected_rss"] != ""]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "certified_matching_cases": len(certified),
        "corrected_wall_wins": sum(1 for ratio in wall if ratio > 1.0),
        "corrected_reported_wins": sum(1 for ratio in reported if ratio > 1.0),
        "corrected_memory_wins": sum(1 for ratio in rss if ratio > 1.0),
        "mean_wall_ratio_adaptive_over_corrected": sum(wall) / len(wall) if wall else None,
        "mean_reported_ratio_adaptive_over_corrected": sum(reported) / len(reported) if reported else None,
        "mean_rss_ratio_adaptive_over_corrected": sum(rss) / len(rss) if rss else None,
        "min_wall_ratio_adaptive_over_corrected": min(wall) if wall else None,
        "max_wall_ratio_adaptive_over_corrected": max(wall) if wall else None,
        "max_corrected_band_count": max([int(row["corrected_band_count"]) for row in completed if row["corrected_band_count"] != ""] or [0]),
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
        "# Residual-Corrected Sums-Only MITM Speed Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Seed: `{report['config']['seed']}`",
        f"- Cases: `{agg['cases']}`",
        f"- Completed cases: `{agg['completed_cases']}`",
        f"- Certified matching cases: `{agg['certified_matching_cases']}`",
        f"- Corrected wall-time wins: `{agg['corrected_wall_wins']}`",
        f"- Corrected reported-time wins: `{agg['corrected_reported_wins']}`",
        f"- Corrected memory wins: `{agg['corrected_memory_wins']}`",
        f"- Mean wall ratio adaptive/corrected: `{agg['mean_wall_ratio_adaptive_over_corrected']}`",
        f"- Mean reported ratio adaptive/corrected: `{agg['mean_reported_ratio_adaptive_over_corrected']}`",
        f"- Mean RSS ratio adaptive/corrected: `{agg['mean_rss_ratio_adaptive_over_corrected']}`",
        f"- Max corrected band count: `{agg['max_corrected_band_count']}`",
        "",
        "This suite compares two modes of the same sums-only MITM kernel. The",
        "corrected mode uses the analytic residual correction to center a small",
        "rank band and then sorts only that band. The adaptive mode is the prior",
        "interpolation/bisection sums-only unrank path. Corrected outputs are",
        "checked with the independent interval-rank auditor, and adaptive outputs",
        "must match the certified corrected exponent vector.",
        "",
        "| label | k | N | corr wall | adaptive wall | wall ratio | corr band | adaptive band | calls | cert | exps |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['N']} | "
            f"{float(row['corrected_wall_seconds']):.6f} | {float(row['adaptive_wall_seconds']):.6f} | "
            f"{float(row['adaptive_wall_over_corrected_wall']):.6f} | "
            f"{row['corrected_band_count']} | {row['adaptive_band_count']} | {row['adaptive_calls']} | "
            f"{row['rank_certified']} | `{row['exps']}` |"
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
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"sums_corrected_speed_suite_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.sums_corrected_speed_suite.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "seed": args.seed,
            "ranks": ranks,
            "k_values": k_values,
            "cases_per_k_rank": args.cases_per_k_rank,
            "cases": cases,
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

    for case in cases:
        report["cases"].append(run_case(harness, case, args.skip_audit))

    rows = [summarize_case(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["certified_matching_cases"] == report["aggregate"]["cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
