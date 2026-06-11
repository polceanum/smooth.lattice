#!/usr/bin/env python3
"""Benchmark a full exponent-vector heap/frontier smooth-number baseline.

The baseline enumerates exponent vectors in value order using a canonical
nondecreasing-prime-index frontier. This is a serious sequential-generation
comparator, not a random-access algorithm: large ranks remain intentionally
capped by the harness.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CASES = (
    ("k3_N1e5", "2,3,5", 100_000),
    ("k5_N1e5", "2,3,5,7,11", 100_000),
    ("k6_N1e5", "2,3,5,7,11,13", 100_000),
    ("k8_N1e5", "2,3,5,7,11,13,17,19", 100_000),
    ("k3_N1e6", "2,3,5", 1_000_000),
    ("k5_N1e6", "2,3,5,7,11", 1_000_000),
    ("k6_N1e6", "2,3,5,7,11,13", 1_000_000),
    ("k8_N1e6", "2,3,5,7,11,13,17,19", 1_000_000),
)


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


def case_from_arg(raw: str) -> dict[str, Any]:
    parts = raw.split(":")
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("case must be label:primes_csv:N")
    return {"label": parts[0], "primes": parts[1], "N": int(parts[2])}


def default_cases() -> list[dict[str, Any]]:
    return [{"label": label, "primes": primes, "N": n} for label, primes, n in DEFAULT_CASES]


def build_heap(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
    return harness.run_measured(
        "build_heap_frontier_baseline",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/heap_frontier_baseline.cpp",
            "-o",
            "bin/heap_frontier_baseline",
        ],
    )


def build_dp(harness, cxx: str, cxxflags: str) -> dict[str, Any]:
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


def current_solver_command(primes: str, n: int) -> tuple[str, list[str]]:
    k = len(parse_csv_ints(primes))
    if k == 3:
        return ("current_beatty3", ["bin/smooth_3prime_beatty_ranker", "nth", primes, str(n)])
    return ("current_sums_adaptive", ["bin/smooth_sums_only_scalable", "nth", primes, str(n)])


def compare(runs: list[dict[str, Any]], n: int) -> dict[str, Any]:
    method_runs = [run for run in runs if run["name"] != "interval_audit"]
    successful = [run for run in method_runs if ok(run)]
    exps_by_method: dict[str, str] = {}
    for run in successful:
        key = exps_key(run.get("metrics", {}).get("exps", ""))
        if key is not None:
            exps_by_method[run["name"]] = key
    matching_exps = None
    all_match = False
    if exps_by_method:
        counts = Counter(exps_by_method.values())
        matching_exps, top_count = counts.most_common(1)[0]
        all_match = top_count == len(exps_by_method)
    audit = next((run for run in runs if run["name"] == "interval_audit"), None)
    certified = bool(
        audit
        and ok(audit)
        and audit.get("metrics", {}).get("rank_certified", False)
        and audit.get("metrics", {}).get("certified_count_le", "") == n
        and all_match
    )
    candidates = [run for run in successful if run["name"] in exps_by_method]
    fastest = min(candidates, key=lambda run: float(run["elapsed_seconds"])) if candidates else None
    return {
        "successful_methods": [run["name"] for run in successful],
        "exps_by_method": exps_by_method,
        "all_successful_exps_match": all_match,
        "matching_exps": matching_exps if all_match else None,
        "certified": certified,
        "fastest_wall_method": fastest["name"] if fastest else "",
        "fastest_wall_seconds": fastest["elapsed_seconds"] if fastest else "",
        "audit_returncode": audit.get("returncode", "") if audit else "",
        "audit_wall_seconds": audit.get("elapsed_seconds", "") if audit else "",
        "audit_stderr": audit.get("stderr", "") if audit else "",
    }


def run_case(harness, case: dict[str, Any], *, dp_max_n: int, skip_audit: bool) -> dict[str, Any]:
    primes = case["primes"]
    n = int(case["N"])
    runs: list[dict[str, Any]] = []
    runs.append(harness.run_measured("heap_frontier", ["bin/heap_frontier_baseline", primes, str(n)]))
    if n <= dp_max_n:
        runs.append(harness.run_measured("dp_pointer", ["bin/dp_pointer_baseline", primes, str(n)]))
    current_name, current_cmd = current_solver_command(primes, n)
    runs.append(harness.run_measured(current_name, current_cmd))
    comparison = compare(runs, n)
    if comparison["matching_exps"] is not None and not skip_audit:
        runs.append(
            harness.run_measured(
                "interval_audit",
                ["bin/smooth_interval_audit_exps", primes, comparison["matching_exps"], str(n)],
            )
        )
        comparison = compare(runs, n)
    return {**case, "runs": runs, "comparison": comparison}


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    runs = {run["name"]: run for run in case["runs"]}
    comparison = case["comparison"]
    heap = runs.get("heap_frontier", {})
    dp = runs.get("dp_pointer", {})
    current = runs.get("current_beatty3") or runs.get("current_sums_adaptive") or {}
    heap_wall = float(heap.get("elapsed_seconds", 0) or 0)
    dp_wall = float(dp.get("elapsed_seconds", 0) or 0)
    current_wall = float(current.get("elapsed_seconds", 0) or 0)
    return {
        "label": case["label"],
        "primes": case["primes"],
        "k": len(parse_csv_ints(case["primes"])),
        "N": case["N"],
        "certified": comparison["certified"],
        "all_successful_exps_match": comparison["all_successful_exps_match"],
        "fastest_wall_method": comparison["fastest_wall_method"],
        "fastest_wall_seconds": comparison["fastest_wall_seconds"],
        "heap_wall_seconds": heap.get("elapsed_seconds", ""),
        "heap_reported_seconds": heap.get("metrics", {}).get("seconds", ""),
        "heap_max_rss_kb": heap.get("max_rss_kb", ""),
        "heap_pushes": heap.get("metrics", {}).get("heap_pushes", ""),
        "heap_max_heap": heap.get("metrics", {}).get("max_heap", ""),
        "dp_wall_seconds": dp.get("elapsed_seconds", ""),
        "dp_reported_seconds": dp.get("metrics", {}).get("seconds", ""),
        "dp_max_rss_kb": dp.get("max_rss_kb", ""),
        "current_method": current.get("name", ""),
        "current_wall_seconds": current.get("elapsed_seconds", ""),
        "current_reported_seconds": current.get("metrics", {}).get("seconds", ""),
        "current_max_rss_kb": current.get("max_rss_kb", ""),
        "heap_over_current_wall_ratio": heap_wall / current_wall if heap_wall and current_wall else "",
        "heap_over_dp_wall_ratio": heap_wall / dp_wall if heap_wall and dp_wall else "",
        "audit_wall_seconds": comparison["audit_wall_seconds"],
        "audit_returncode": comparison["audit_returncode"],
        "audit_stderr": comparison["audit_stderr"],
        "exps": comparison["matching_exps"] or "",
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    certified = [row for row in rows if row["certified"]]
    matching = [row for row in rows if row["all_successful_exps_match"]]
    winners = Counter(row["fastest_wall_method"] for row in certified)
    ratios = [float(row["heap_over_current_wall_ratio"]) for row in certified if row["heap_over_current_wall_ratio"] != ""]
    dp_ratios = [float(row["heap_over_dp_wall_ratio"]) for row in certified if row["heap_over_dp_wall_ratio"] != ""]
    return {
        "cases": len(rows),
        "certified_cases": len(certified),
        "matching_cases": len(matching),
        "solver_disagreement_cases": len(rows) - len(matching),
        "wall_winners": dict(sorted(winners.items())),
        "mean_heap_over_current_wall_ratio": sum(ratios) / len(ratios) if ratios else "",
        "mean_heap_over_dp_wall_ratio": sum(dp_ratios) / len(dp_ratios) if dp_ratios else "",
        "all_cases_certified": len(certified) == len(rows),
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
            "heap_wall_seconds",
            "dp_wall_seconds",
            "current_method",
            "current_wall_seconds",
            "heap_over_current_wall_ratio",
            "exps",
        ]
        writer = csv.DictWriter(
            f,
            fieldnames=preferred + [key for key in fieldnames if key not in preferred],
            lineterminator="\n",
            extrasaction="ignore",
        )
        writer.writeheader()
        writer.writerows(rows)

    agg = report["aggregate"]
    lines = [
        "# Heap Frontier Baseline Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Cases: `{agg['cases']}`",
        f"- Certified cases: `{agg['certified_cases']}`",
        f"- Solver-agreement cases: `{agg['matching_cases']}`",
        f"- Solver-disagreement cases: `{agg['solver_disagreement_cases']}`",
        f"- Wall-time winners: `{agg['wall_winners']}`",
        f"- Mean heap/current wall ratio: `{agg['mean_heap_over_current_wall_ratio']}`",
        f"- Mean heap/DP wall ratio: `{agg['mean_heap_over_dp_wall_ratio']}`",
        "",
        "The heap row is a canonical frontier generator that returns the full",
        "exponent vector. A row is certified only when the heap, DP baseline when",
        "enabled, and current solver agree on the exponent vector and the",
        "independent interval auditor proves `count_le(candidate)=N`.",
        "",
        "| label | k | N | fastest | heap s | DP s | current s | heap/current | cert | exps |",
        "|---|---:|---:|---|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        dp_s = "" if row["dp_wall_seconds"] == "" else f"{float(row['dp_wall_seconds']):.6f}"
        ratio = "" if row["heap_over_current_wall_ratio"] == "" else f"{float(row['heap_over_current_wall_ratio']):.4f}"
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['N']} | `{row['fastest_wall_method']}` | "
            f"{float(row['heap_wall_seconds']):.6f} | {dp_s} | "
            f"{float(row['current_wall_seconds']):.6f} | {ratio} | {row['certified']} | `{row['exps']}` |"
        )
    lines.extend(
        [
            "",
            "## Claim Status",
            "",
            "This suite discharges the repository's heap/frontier baseline obligation",
            "for the tested small-to-medium ranks. It is not a broad best-known",
            "claim, and it does not replace the still-open soft-heap, full",
            "Frederickson-Johnson, or Barvinok-style comparator obligations.",
        ]
    )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=case_from_arg, dest="cases")
    parser.add_argument("--dp-max-n", type=int, default=1_000_000)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    harness = load_harness()
    cases = args.cases or default_cases()
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"heap_frontier_baseline_suite_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.heap_frontier_baseline_suite.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "cases": cases,
            "dp_max_n": args.dp_max_n,
            "skip_audit": args.skip_audit,
        },
        "build_runs": [],
        "cases": [],
    }

    for build in (
        harness.run_measured("build_core", ["./build.sh"]),
        build_heap(harness, args.cxx, args.cxxflags),
        build_dp(harness, args.cxx, args.cxxflags),
    ):
        report["build_runs"].append(build)
        if not ok(build):
            report["summary_rows"] = []
            report["aggregate"] = {"status": build["name"] + "_failed"}
            write_outputs(out_dir, report)
            print(f"{build['name']} failed; wrote partial report to {out_dir}")
            return 1

    for case in cases:
        report["cases"].append(run_case(harness, case, dp_max_n=args.dp_max_n, skip_audit=args.skip_audit))

    report["summary_rows"] = [summarize_case(case) for case in report["cases"]]
    report["aggregate"] = aggregate(report["summary_rows"])
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["all_cases_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
