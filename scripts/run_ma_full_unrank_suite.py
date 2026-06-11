#!/usr/bin/env python3
"""Run the full Mirzaian-Arjomandi X+Y unrank comparator.

This suite compares a faithful Mirzaian-Arjomandi sorted-matrix value selector,
wrapped with exact exponent reconstruction, against the current
analytic-corrected X+Y unranker.  It is intentionally narrow: both rows return a
full exponent vector and every matching row is independently rank-audited.
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


def result_ok(run: dict[str, Any] | None) -> bool:
    return run is not None and int(run.get("returncode", 1)) == 0


def exps_csv(exps: Any) -> str | None:
    if isinstance(exps, list) and all(isinstance(item, int) for item in exps):
        return ",".join(str(item) for item in exps)
    if isinstance(exps, str) and exps:
        stripped = exps.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            return stripped[1:-1].replace(" ", "")
    return None


def case_slug(primes: str) -> str:
    return "P" + primes.replace(",", "_")


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


def audit_exps(single, primes: str, exps: Any, n: int) -> dict[str, Any] | None:
    exps_text = exps_csv(exps)
    if exps_text is None:
        return None
    return single.run_measured(
        "interval_audit",
        ["bin/smooth_interval_audit_exps", primes, exps_text, str(n)],
    )


def run_case(
    single,
    primes: str,
    n: int,
    rank_radius: int,
    max_candidates: int,
    max_n: int,
    max_middle: int,
    skip_audit: bool,
) -> dict[str, Any]:
    runs = []
    corrected = single.run_measured(
        "xplusy_corrected",
        [
            "bin/smooth_xplusy_full_unrank",
            "analytic-band-corrected",
            primes,
            str(n),
            str(rank_radius),
            str(max_candidates),
        ],
    )
    runs.append(corrected)
    ma_full = single.run_measured(
        "xplusy_ma_full",
        [
            "bin/smooth_xplusy_full_unrank",
            "ma-full",
            primes,
            str(n),
            str(rank_radius),
            str(max_candidates),
            str(max_n),
            str(max_middle),
        ],
    )
    runs.append(ma_full)

    audit = None
    if not skip_audit and result_ok(corrected) and result_ok(ma_full):
        if corrected.get("metrics", {}).get("exps") == ma_full.get("metrics", {}).get("exps"):
            audit = audit_exps(single, primes, corrected.get("metrics", {}).get("exps"), n)
            if audit is not None:
                runs.append(audit)

    return {
        "primes": primes,
        "N": n,
        "slug": case_slug(primes),
        "runs": runs,
        "comparison": compute_comparison(runs, n),
    }


def compute_comparison(runs: list[dict[str, Any]], n: int) -> dict[str, Any]:
    by_name = {run["name"]: run for run in runs}
    corrected = by_name.get("xplusy_corrected")
    ma_full = by_name.get("xplusy_ma_full")
    audit = by_name.get("interval_audit")
    comparison: dict[str, Any] = {
        "same_exps": False,
        "certified": False,
        "ma_wall_over_corrected_wall": "",
        "ma_reported_over_corrected_reported": "",
    }
    if result_ok(corrected) and result_ok(ma_full):
        corrected_exps = corrected.get("metrics", {}).get("exps")
        ma_exps = ma_full.get("metrics", {}).get("exps")
        comparison["same_exps"] = corrected_exps == ma_exps
        corrected_wall = float(corrected["elapsed_seconds"])
        ma_wall = float(ma_full["elapsed_seconds"])
        if corrected_wall > 0:
            comparison["ma_wall_over_corrected_wall"] = ma_wall / corrected_wall
        corrected_sec = corrected.get("metrics", {}).get("seconds", "")
        ma_sec = ma_full.get("metrics", {}).get("seconds", "")
        if corrected_sec != "" and ma_sec != "" and float(corrected_sec) > 0:
            comparison["ma_reported_over_corrected_reported"] = float(ma_sec) / float(corrected_sec)
    comparison["certified"] = bool(
        audit
        and result_ok(audit)
        and audit.get("metrics", {}).get("rank_certified", False)
        and audit.get("metrics", {}).get("certified_count_le", "") == n
        and comparison["same_exps"]
    )
    return comparison


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in case["runs"]}
    corrected = by_name.get("xplusy_corrected", {})
    ma_full = by_name.get("xplusy_ma_full", {})
    audit = by_name.get("interval_audit", {})
    comparison = case["comparison"]
    return {
        "primes": case["primes"],
        "N": case["N"],
        "corrected_returncode": corrected.get("returncode", ""),
        "ma_returncode": ma_full.get("returncode", ""),
        "audit_returncode": audit.get("returncode", ""),
        "same_exps": comparison["same_exps"],
        "certified": comparison["certified"],
        "corrected_wall_seconds": corrected.get("elapsed_seconds", ""),
        "ma_wall_seconds": ma_full.get("elapsed_seconds", ""),
        "audit_wall_seconds": audit.get("elapsed_seconds", ""),
        "corrected_reported_seconds": corrected.get("metrics", {}).get("seconds", ""),
        "ma_reported_seconds": ma_full.get("metrics", {}).get("seconds", ""),
        "ma_phase_seconds": ma_full.get("metrics", {}).get("ma_phase", ""),
        "corrected_max_rss_kb": corrected.get("max_rss_kb", ""),
        "ma_max_rss_kb": ma_full.get("max_rss_kb", ""),
        "ma_wall_over_corrected_wall": comparison["ma_wall_over_corrected_wall"],
        "ma_reported_over_corrected_reported": comparison["ma_reported_over_corrected_reported"],
        "ma_A": ma_full.get("metrics", {}).get("A", ""),
        "ma_B": ma_full.get("metrics", {}).get("B", ""),
        "ma_n_square": ma_full.get("metrics", {}).get("n_square", ""),
        "ma_padded_a": ma_full.get("metrics", {}).get("padded_a", ""),
        "ma_padded_b": ma_full.get("metrics", {}).get("padded_b", ""),
        "ma_band_count": ma_full.get("metrics", {}).get("band_count", ""),
        "ma_recovered": ma_full.get("metrics", {}).get("recovered", ""),
        "ma_skipped": ma_full.get("metrics", {}).get("ma_skipped", ""),
        "corrected_band_count": corrected.get("metrics", {}).get("band_count", ""),
        "exps": corrected.get("metrics", {}).get("exps", ""),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["corrected_returncode"] == 0 and row["ma_returncode"] == 0]
    certified = [row for row in completed if row["same_exps"] and row["certified"]]
    ratios = [
        float(row["ma_wall_over_corrected_wall"])
        for row in certified
        if row["ma_wall_over_corrected_wall"] != ""
    ]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "same_exps_certified_cases": len(certified),
        "ma_wall_wins": sum(1 for row in certified if float(row["ma_wall_over_corrected_wall"]) < 1.0),
        "all_completed": len(completed) == len(rows),
        "all_same_exps_certified": len(certified) == len(rows),
        "min_ma_wall_over_corrected_wall": min(ratios) if ratios else None,
        "max_ma_wall_over_corrected_wall": max(ratios) if ratios else None,
        "mean_ma_wall_over_corrected_wall": sum(ratios) / len(ratios) if ratios else None,
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
        "# MA Full X+Y Unrank Suite",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- N: `{report['config']['N']}`",
        f"- Cases: `{report['aggregate']['cases']}`",
        f"- Completed cases: `{report['aggregate']['completed_cases']}`",
        f"- Same-exponent certified cases: `{report['aggregate']['same_exps_certified_cases']}`",
        f"- MA wall-time wins: `{report['aggregate']['ma_wall_wins']}`",
        f"- Mean MA/corrected wall ratio: `{report['aggregate']['mean_ma_wall_over_corrected_wall']}`",
        "",
        "Both rows return full exponent vectors. `xplusy_ma_full` uses the",
        "Mirzaian-Arjomandi sorted-matrix selector to choose the X+Y log value,",
        "then reconstructs and exact-sorts a narrow exponent band. Matching rows",
        "are checked by the independent interval-rank auditor.",
        "",
        "| P | corrected wall s | MA wall s | ratio | MA phase s | corrected RSS KB | MA RSS KB | MA n | MA band | certified | exps |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|",
    ]
    def fmt_float(value: Any) -> str:
        if value == "" or value is None:
            return ""
        return f"{float(value):.6f}"

    for row in rows:
        ratio = row["ma_wall_over_corrected_wall"]
        lines.append(
            f"| `{row['primes']}` | {fmt_float(row['corrected_wall_seconds'])} | "
            f"{fmt_float(row['ma_wall_seconds'])} | {fmt_float(ratio)} | "
            f"{fmt_float(row['ma_phase_seconds'])} | {row['corrected_max_rss_kb']} | "
            f"{row['ma_max_rss_kb']} | {row['ma_n_square']} | {row['ma_band_count']} | "
            f"{row['certified']} | `{row['exps']}` |"
        )
    lines.extend(
        [
            "",
            "## Claim Status",
            "",
            "This is a published sorted-matrix/X+Y selector comparison wrapped into",
            "full exponent-vector unranking and independent rank auditing. It is a",
            "narrow best-known-style comparator for the sorted-matrix value-selection",
            "subproblem, not a claim against soft-heap top-k/output algorithms, full",
            "Frederickson-Johnson ranking/selection, or Barvinok-style lattice",
            "counting.",
        ]
    )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--prime-set", action="append", dest="prime_sets")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--rank-radius", type=int, default=250)
    parser.add_argument("--max-candidates", type=int, default=200_000)
    parser.add_argument("--ma-max-n", type=int, default=30_000_000)
    parser.add_argument("--ma-max-middle", type=int, default=200_000_000)
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    single = load_single_harness()
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"ma_full_unrank_suite_{utc_stamp()}")
    prime_sets = args.prime_sets or DEFAULT_PRIME_SETS
    build = build_full_xplusy(single, args.cxx, args.cxxflags)
    cases = [
        run_case(
            single,
            primes,
            args.N,
            args.rank_radius,
            args.max_candidates,
            args.ma_max_n,
            args.ma_max_middle,
            args.skip_audit,
        )
        for primes in prime_sets
    ]
    rows = [summarize_case(case) for case in cases]
    report = {
        "schema": "smooth.lattice.ma_full_unrank_suite.v1",
        "metadata": single.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "N": args.N,
            "prime_sets": prime_sets,
            "rank_radius": args.rank_radius,
            "max_candidates": args.max_candidates,
            "ma_max_n": args.ma_max_n,
            "ma_max_middle": args.ma_max_middle,
            "skip_audit": args.skip_audit,
        },
        "build_runs": [build],
        "cases": cases,
        "summary_rows": rows,
        "aggregate": aggregate(rows),
    }
    write_outputs(out_dir, report)
    print(f"wrote {out_dir}")
    return 0 if report["aggregate"]["all_same_exps_certified"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
