#!/usr/bin/env python3
"""Run sorted-matrix/FJ-style and LOH probes for X+Y baselines.

This is deliberately a workbench, not a headline baseline. The range-pruned
block counter is a practical sorted-matrix probe, not a faithful implementation
of Frederickson-Johnson selection. The LOH row is an output-style top-k probe
and is only run at a capped rank.
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


def result_ok(result: dict[str, Any]) -> bool:
    return int(result.get("returncode", 1)) == 0


def build_workbench(single, cxx: str, cxxflags: str) -> dict[str, Any]:
    single.BIN.mkdir(exist_ok=True)
    return single.run_measured(
        "build_sorted_matrix_workbench",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/smooth_xplusy_fj_loh_workbench.cpp",
            "-o",
            "bin/smooth_xplusy_fj_loh_workbench",
        ],
    )


def parse_workbench_stdout(single, stdout: str) -> list[dict[str, Any]]:
    rows = []
    for line in stdout.splitlines():
        line = line.strip()
        if not line:
            continue
        row_type = line.split(None, 1)[0]
        metrics = single.parse_kv(line)
        metrics["row_type"] = row_type
        rows.append(metrics)
    return rows


def run_case(single, primes: str, n: int) -> dict[str, Any]:
    run = single.run_measured(
        "sorted_matrix_workbench",
        ["bin/smooth_xplusy_fj_loh_workbench", primes, str(n)],
    )
    rows = parse_workbench_stdout(single, run.get("stdout", ""))
    return {
        "primes": primes,
        "N": n,
        "run": run,
        "rows": rows,
        "summary": summarize_rows(rows),
    }


def summarize_rows(rows: list[dict[str, Any]]) -> dict[str, Any]:
    block_rows = [row for row in rows if row.get("row_type") == "block_rank"]
    loh_rows = [row for row in rows if row.get("row_type") == "loh_topk_probe"]
    summary: dict[str, Any] = {
        "block_rows": len(block_rows),
        "loh_rows": len(loh_rows),
    }
    if block_rows:
        best_block = min(block_rows, key=lambda row: float(row.get("block_total", float("inf"))))
        linear_total = float(best_block["linear_total"])
        block_total = float(best_block["block_total"])
        summary.update(
            {
                "best_leaf": best_block.get("leaf"),
                "linear_total": linear_total,
                "best_block_total": block_total,
                "block_over_linear": block_total / linear_total if linear_total > 0 else None,
                "range_pruning_wins": block_total < linear_total,
                "block_log_delta": best_block.get("block_log_delta"),
                "A": best_block.get("A"),
                "B": best_block.get("B"),
            }
        )
    if loh_rows:
        loh = loh_rows[0]
        summary.update(
            {
                "loh_N_probe": loh.get("N_probe"),
                "loh_skipped": loh.get("skipped"),
                "loh_sec": loh.get("sec"),
                "loh_cand": loh.get("cand"),
                "loh_pairs": loh.get("pairs"),
                "loh_reason": loh.get("reason", ""),
            }
        )
    return summary


def summary_row(case: dict[str, Any]) -> dict[str, Any]:
    run = case["run"]
    summary = case["summary"]
    return {
        "primes": case["primes"],
        "N": case["N"],
        "returncode": run.get("returncode"),
        "wall_seconds": run.get("elapsed_seconds", ""),
        "max_rss_kb": run.get("max_rss_kb", ""),
        "A": summary.get("A", ""),
        "B": summary.get("B", ""),
        "linear_total": summary.get("linear_total", ""),
        "best_block_total": summary.get("best_block_total", ""),
        "block_over_linear": summary.get("block_over_linear", ""),
        "best_leaf": summary.get("best_leaf", ""),
        "range_pruning_wins": summary.get("range_pruning_wins", ""),
        "block_log_delta": summary.get("block_log_delta", ""),
        "loh_N_probe": summary.get("loh_N_probe", ""),
        "loh_skipped": summary.get("loh_skipped", ""),
        "loh_sec": summary.get("loh_sec", ""),
        "loh_cand": summary.get("loh_cand", ""),
        "loh_pairs": summary.get("loh_pairs", ""),
        "loh_reason": summary.get("loh_reason", ""),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["returncode"] == 0]
    comparable = [row for row in completed if row["block_over_linear"] != ""]
    wins = [row for row in comparable if row["range_pruning_wins"] is True]
    ratios = [float(row["block_over_linear"]) for row in comparable]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "comparable_cases": len(comparable),
        "range_pruning_wins": len(wins),
        "all_completed": len(completed) == len(rows),
        "min_block_over_linear": min(ratios) if ratios else None,
        "max_block_over_linear": max(ratios) if ratios else None,
        "mean_block_over_linear": sum(ratios) / len(ratios) if ratios else None,
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
        "# Sorted-Matrix / LOH X+Y Workbench",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- N: `{report['config']['N']}`",
        f"- Cases: `{report['aggregate']['cases']}`",
        f"- Completed cases: `{report['aggregate']['completed_cases']}`",
        f"- Range-pruning wins over linear saddleback count: `{report['aggregate']['range_pruning_wins']}`",
        f"- Mean block/linear internal time ratio: `{report['aggregate']['mean_block_over_linear']}`",
        "",
        "The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful",
        "Frederickson-Johnson implementation. The `loh_topk_probe` row uses a capped",
        "output-style top-k rank (`N_probe`) and is not evidence for random access at",
        "the full rank when `N_probe != N`.",
        "",
        "| P | linear total s | best block s | block/linear | best leaf | wall s | RSS KB | LOH probe N | LOH s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['primes']}` | {float(row['linear_total']):.6f} | "
            f"{float(row['best_block_total']):.6f} | "
            f"{float(row['block_over_linear']):.6f} | "
            f"{row['best_leaf']} | {float(row['wall_seconds']):.6f} | "
            f"{row['max_rss_kb']} | {row['loh_N_probe']} | {float(row['loh_sec']):.6f} |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--prime-set", action="append", dest="prime_sets")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    args = parser.parse_args()

    single = load_single_harness()
    prime_sets = args.prime_sets or DEFAULT_PRIME_SETS
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"sorted_matrix_workbench_{utc_stamp()}")

    report: dict[str, Any] = {
        "schema": "smooth.lattice.sorted_matrix_workbench.v1",
        "metadata": single.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "N": args.N,
            "prime_sets": prime_sets,
        },
        "build_runs": [],
        "cases": [],
    }

    build = build_workbench(single, args.cxx, args.cxxflags)
    report["build_runs"].append(build)
    if not result_ok(build):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"build failed; wrote partial report to {out_dir}")
        return 1

    for primes in prime_sets:
        report["cases"].append(run_case(single, primes, args.N))

    rows = [summary_row(case) for case in report["cases"]]
    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["all_completed"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
