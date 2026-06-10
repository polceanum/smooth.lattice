#!/usr/bin/env python3
"""Run a certified portfolio benchmark across implemented unrank methods.

This harness is intentionally stricter than single-method speed suites: it runs
all applicable full-unrank implementations for each case, requires matching
exponent vectors among successful methods, audits one matching vector, and then
declares the fastest certified method for that case.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import itertools
import json
import random
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
SUPPORTED_PRIMES = (2, 3, 5, 7, 11, 13, 17, 19)
DEFAULT_SEED = 20260609
DEFAULT_RANKS = (1_000_000_000, 1_000_000_000_000)
DEFAULT_K_VALUES = (3, 5, 6, 8)
DEFAULT_CASES_PER_K_RANK = 2
DEFAULT_MAX_CANDIDATES = 200_000
DEFAULT_TARGET_GAP = 50_000
RADIUS_BY_K = {3: 100, 5: 250, 6: 1500, 7: 5000, 8: 10000}


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


def exps_key(value: Any) -> str | None:
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
            if k == 3:
                combos = [(2, 3, 5), (2, 3, 7), (3, 5, 7)]
            else:
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


def build_full_xplusy(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
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


def build_dp_pointer(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
    return harness.run_measured(
        "build_dp_pointer_baseline",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/dp_pointer_baseline.cpp",
            "-o",
            "bin/dp_pointer_baseline",
        ],
    )


def applicable_methods(case: dict[str, Any], target_gap: int, dp_max_n: int) -> list[tuple[str, list[str]]]:
    primes = case["primes"]
    n = int(case["N"])
    k = len(parse_csv_ints(primes))
    radius = str(int(case["rank_radius"]))
    cap = str(int(case["max_candidates"]))
    methods: list[tuple[str, list[str]]] = []
    if k == 3:
        methods.append(("beatty3", ["bin/smooth_3prime_beatty_ranker", "nth", primes, str(n)]))
    if k >= 5:
        methods.extend(
            [
                ("layer_adaptive", ["bin/smooth_layer_compressed_general", "nth", primes, str(n), str(target_gap)]),
                (
                    "layer_corrected",
                    ["bin/smooth_layer_compressed_general", "analytic-band-corrected", primes, str(n), radius, cap],
                ),
            ]
        )
    methods.extend(
        [
            (
                "sums_corrected",
                ["bin/smooth_sums_only_scalable", "analytic-band-corrected", primes, str(n), radius, cap],
            ),
            ("sums_adaptive", ["bin/smooth_sums_only_scalable", "nth", primes, str(n)]),
            (
                "xplusy_corrected",
                ["bin/smooth_xplusy_full_unrank", "analytic-band-corrected", primes, str(n), radius, cap],
            ),
            ("xplusy_adaptive", ["bin/smooth_xplusy_full_unrank", "nth", primes, str(n), str(target_gap)]),
        ]
    )
    if n <= dp_max_n:
        methods.append(("dp_pointer", ["bin/dp_pointer_baseline", primes, str(n)]))
    return methods


def run_case(
    harness,
    case: dict[str, Any],
    *,
    target_gap: int,
    dp_max_n: int,
    skip_audit: bool,
) -> dict[str, Any]:
    runs = [
        harness.run_measured(name, command)
        for name, command in applicable_methods(case, target_gap, dp_max_n)
    ]
    comparison = compare_methods(runs, int(case["N"]))
    audit = None
    if comparison["matching_exps"] is not None and not skip_audit:
        audit = harness.run_measured(
            "portfolio_interval_audit",
            ["bin/smooth_interval_audit_exps", case["primes"], comparison["matching_exps"], str(case["N"])],
        )
        runs.append(audit)
        comparison = compare_methods(runs, int(case["N"]))
    return {
        "label": case["label"],
        "primes": case["primes"],
        "N": case["N"],
        "rank_radius": case["rank_radius"],
        "max_candidates": case["max_candidates"],
        "runs": runs,
        "comparison": comparison,
    }


def compare_methods(runs: list[dict[str, Any]], n: int) -> dict[str, Any]:
    method_runs = [run for run in runs if run["name"] != "portfolio_interval_audit"]
    successful = [run for run in method_runs if ok(run)]
    exps_by_method: dict[str, str] = {}
    for run in successful:
        key = exps_key(run.get("metrics", {}).get("exps", ""))
        if key is not None:
            exps_by_method[run["name"]] = key
    matching_exps = None
    exps_match = False
    if exps_by_method:
        counts = Counter(exps_by_method.values())
        matching_exps, top_count = counts.most_common(1)[0]
        exps_match = top_count == len(exps_by_method)

    candidates = [run for run in successful if run["name"] in exps_by_method]
    fastest_wall = min(candidates, key=lambda run: float(run["elapsed_seconds"])) if candidates else None
    fastest_reported = min(
        [run for run in candidates if run.get("metrics", {}).get("seconds", "") != ""],
        key=lambda run: float(run["metrics"]["seconds"]),
        default=None,
    )
    audit = next((run for run in runs if run["name"] == "portfolio_interval_audit"), None)
    certified = bool(
        audit
        and ok(audit)
        and audit.get("metrics", {}).get("rank_certified", False)
        and audit.get("metrics", {}).get("certified_count_le", "") == n
    )
    audit_stderr = audit.get("stderr", "") if audit else ""
    audit_blocked = bool(audit and not ok(audit) and exps_match and "scaled interval overflow" in audit_stderr)
    return {
        "successful_methods": [run["name"] for run in successful],
        "exps_by_method": exps_by_method,
        "matching_exps": matching_exps if exps_match else None,
        "all_successful_exps_match": exps_match,
        "certified": certified and exps_match,
        "audit_blocked": audit_blocked,
        "audit_stderr": audit_stderr,
        "fastest_wall_method": fastest_wall["name"] if fastest_wall else "",
        "fastest_wall_seconds": fastest_wall["elapsed_seconds"] if fastest_wall else "",
        "fastest_reported_method": fastest_reported["name"] if fastest_reported else "",
        "fastest_reported_seconds": fastest_reported["metrics"]["seconds"] if fastest_reported else "",
        "audit_returncode": audit.get("returncode", "") if audit else "",
        "audit_wall_seconds": audit.get("elapsed_seconds", "") if audit else "",
        "cert_seconds": audit.get("metrics", {}).get("cert_seconds", "") if audit else "",
    }


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    comparison = case["comparison"]
    runs_by_name = {run["name"]: run for run in case["runs"]}
    method_seconds: dict[str, Any] = {}
    method_wall: dict[str, Any] = {}
    method_rss: dict[str, Any] = {}
    for name, run in runs_by_name.items():
        if name == "portfolio_interval_audit":
            continue
        method_wall[f"{name}_wall_seconds"] = run.get("elapsed_seconds", "")
        method_seconds[f"{name}_reported_seconds"] = run.get("metrics", {}).get("seconds", "")
        method_rss[f"{name}_max_rss_kb"] = run.get("max_rss_kb", "")
    return {
        "label": case["label"],
        "primes": case["primes"],
        "k": len(parse_csv_ints(case["primes"])),
        "N": case["N"],
        "rank_radius": case["rank_radius"],
        "max_candidates": case["max_candidates"],
        "completed_methods": ",".join(comparison["successful_methods"]),
        "all_successful_exps_match": comparison["all_successful_exps_match"],
        "certified": comparison["certified"],
        "audit_blocked": comparison["audit_blocked"],
        "audit_stderr": comparison["audit_stderr"],
        "fastest_wall_method": comparison["fastest_wall_method"],
        "fastest_wall_seconds": comparison["fastest_wall_seconds"],
        "fastest_reported_method": comparison["fastest_reported_method"],
        "fastest_reported_seconds": comparison["fastest_reported_seconds"],
        "audit_wall_seconds": comparison["audit_wall_seconds"],
        "cert_seconds": comparison["cert_seconds"],
        "exps": comparison["matching_exps"] or "",
        **method_wall,
        **method_seconds,
        **method_rss,
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    certified = [row for row in rows if row["certified"]]
    audit_blocked = [row for row in rows if row["audit_blocked"]]
    matching = [row for row in rows if row["all_successful_exps_match"]]
    disagreements = [row for row in rows if not row["all_successful_exps_match"]]
    wall_winners = Counter(row["fastest_wall_method"] for row in certified)
    reported_winners = Counter(row["fastest_reported_method"] for row in certified)
    return {
        "cases": len(rows),
        "certified_cases": len(certified),
        "matching_cases": len(matching),
        "audit_blocked_cases": len(audit_blocked),
        "solver_disagreement_cases": len(disagreements),
        "wall_winners": dict(sorted(wall_winners.items())),
        "reported_winners": dict(sorted(reported_winners.items())),
        "all_cases_certified": len(certified) == len(rows),
        "all_cases_resolved": len(disagreements) == 0 and len(certified) + len(audit_blocked) == len(rows),
    }


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    rows = report["summary_rows"]
    with (out_dir / "summary.csv").open("w", newline="") as f:
        fieldnames = sorted({key for row in rows for key in row.keys()})
        preferred = [
            "label",
            "primes",
            "k",
            "N",
            "certified",
            "fastest_wall_method",
            "fastest_wall_seconds",
            "fastest_reported_method",
            "fastest_reported_seconds",
            "exps",
        ]
        fieldnames = preferred + [key for key in fieldnames if key not in preferred]
        writer = csv.DictWriter(f, fieldnames=fieldnames, lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)

    agg = report["aggregate"]
    lines = [
        "# Certified Portfolio Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Cases: `{agg['cases']}`",
        f"- Certified cases: `{agg['certified_cases']}`",
        f"- Solver-agreement cases: `{agg['matching_cases']}`",
        f"- Audit-blocked agreement cases: `{agg['audit_blocked_cases']}`",
        f"- Solver-disagreement cases: `{agg['solver_disagreement_cases']}`",
        f"- Wall-time winners: `{agg['wall_winners']}`",
        f"- Reported-time winners: `{agg['reported_winners']}`",
        "",
        "This suite is a portfolio benchmark over implemented full-unrank methods.",
        "A row is certified only when successful methods agree on the exponent vector",
        "and the independent interval auditor proves `count_le(candidate)=N`.",
        "",
        "| label | k | N | fastest wall | wall sec | fastest reported | reported sec | cert | audit blocked | exps |",
        "|---|---:|---:|---|---:|---|---:|---|---|---|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['N']} | "
            f"`{row['fastest_wall_method']}` | {float(row['fastest_wall_seconds']):.6f} | "
            f"`{row['fastest_reported_method']}` | {float(row['fastest_reported_seconds']):.6f} | "
            f"{row['certified']} | {row['audit_blocked']} | `{row['exps']}` |"
        )
    lines.extend(
        [
            "",
            "## Best-Known Claim Gate",
            "",
            "This artifact is not, by itself, a broad best-known claim. It is the",
            "current implemented portfolio with per-row certification status.",
            "Rows marked `audit blocked` had solver agreement but exceeded the",
            "current interval auditor's numeric range. Remaining comparator gates",
            "before using best-known language include a faithful soft-heap/selection",
            "baseline, clearer Barvinok/fixed-dimensional lattice-count positioning,",
            "and broader independent certification beyond the current audited prime",
            "universe.",
        ]
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
    parser.add_argument("--target-gap", type=int, default=DEFAULT_TARGET_GAP)
    parser.add_argument("--dp-max-n", type=int, default=1_000_000)
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
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"certified_portfolio_suite_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.certified_portfolio_suite.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "seed": args.seed,
            "ranks": ranks,
            "k_values": k_values,
            "cases_per_k_rank": args.cases_per_k_rank,
            "cases": cases,
            "target_gap": args.target_gap,
            "dp_max_n": args.dp_max_n,
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

    build_dp = build_dp_pointer(harness, args.cxx, args.cxxflags)
    report["build_runs"].append(build_dp)
    if not ok(build_dp):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "dp_pointer_build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"DP pointer build failed; wrote partial report to {out_dir}")
        return 1

    for case in cases:
        report["cases"].append(
            run_case(
                harness,
                case,
                target_gap=args.target_gap,
                dp_max_n=args.dp_max_n,
                skip_audit=args.skip_audit,
            )
        )

    rows = [summarize_case(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["all_cases_resolved"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
