#!/usr/bin/env python3
"""Run sorted-matrix/FJ-style and LOH probes for X+Y baselines.

This is deliberately a workbench, not a headline baseline. The range-pruned
block counter is a practical sorted-matrix probe, not a faithful implementation
of Frederickson-Johnson selection. The Mirzaian-Arjomandi row is a value
selection probe, not full exponent-vector unranking. The `matselect2_heap` row
implements the Kaplan/Frederickson-Johnson-style exponential-block row-sorted
selector using an exact binary heap for the Mat-Select1 primitive; it is not a
soft-heap time-bound implementation. The LOH row is an output-style top-k probe
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
    ma_rows = [row for row in rows if row.get("row_type") == "ma_select_probe"]
    mat_rows = [row for row in rows if row.get("row_type") == "matselect2_heap_probe"]
    soft_rows = [row for row in rows if row.get("row_type") == "matselect2_soft_probe"]
    loh_rows = [row for row in rows if row.get("row_type") == "loh_topk_probe"]
    summary: dict[str, Any] = {
        "block_rows": len(block_rows),
        "ma_rows": len(ma_rows),
        "matselect2_rows": len(mat_rows),
        "matselect2_soft_rows": len(soft_rows),
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
    if ma_rows:
        ma = ma_rows[0]
        ma_sec = ma.get("sec", "")
        linear_total = summary.get("linear_total")
        ma_over_linear = ""
        if ma_sec != "" and linear_total not in {None, "", 0}:
            ma_over_linear = float(ma_sec) / float(linear_total)
        summary.update(
            {
                "ma_skipped": ma.get("skipped"),
                "ma_sec": ma_sec,
                "ma_over_linear": ma_over_linear,
                "ma_log_delta": ma.get("log_delta"),
                "ma_n_square": ma.get("n_square"),
                "ma_padded_a": ma.get("padded_a"),
                "ma_padded_b": ma.get("padded_b"),
                "ma_reason": ma.get("reason", ""),
                "ma_wins_linear": bool(ma_over_linear != "" and ma_over_linear < 1.0),
            }
        )
    if mat_rows:
        mat = mat_rows[0]
        mat_sec = mat.get("sec", "")
        linear_total = summary.get("linear_total")
        mat_over_linear = ""
        if mat_sec != "" and linear_total not in {None, "", 0} and mat.get("skipped") is False:
            mat_over_linear = float(mat_sec) / float(linear_total)
        summary.update(
            {
                "matselect2_skipped": mat.get("skipped"),
                "matselect2_sec": mat_sec,
                "matselect2_over_linear": mat_over_linear,
                "matselect2_log_delta": mat.get("log_delta"),
                "matselect2_iterations": mat.get("iterations"),
                "matselect2_rep_pops": mat.get("rep_pops"),
                "matselect2_base_pops": mat.get("base_pops"),
                "matselect2_removed": mat.get("removed"),
                "matselect2_max_active_rows": mat.get("max_active_rows"),
                "matselect2_reason": mat.get("reason", ""),
                "matselect2_wins_linear": bool(mat_over_linear != "" and mat_over_linear < 1.0),
            }
        )
    if soft_rows:
        soft = soft_rows[0]
        soft_sec = soft.get("sec", "")
        linear_total = summary.get("linear_total")
        soft_over_linear = ""
        if soft_sec != "" and linear_total not in {None, "", 0} and soft.get("skipped") is False:
            soft_over_linear = float(soft_sec) / float(linear_total)
        summary.update(
            {
                "matselect2_soft_skipped": soft.get("skipped"),
                "matselect2_soft_sec": soft_sec,
                "matselect2_soft_over_linear": soft_over_linear,
                "matselect2_soft_log_delta": soft.get("log_delta"),
                "matselect2_soft_iterations": soft.get("iterations"),
                "matselect2_soft_extracts": soft.get("extracts"),
                "matselect2_soft_inserted": soft.get("inserted"),
                "matselect2_soft_candidates": soft.get("candidates"),
                "matselect2_soft_corrupt_returns": soft.get("corrupt_returns"),
                "matselect2_soft_newly_corrupt": soft.get("newly_corrupt"),
                "matselect2_soft_removed": soft.get("removed"),
                "matselect2_soft_max_active_rows": soft.get("max_active_rows"),
                "matselect2_soft_reason": soft.get("reason", ""),
                "matselect2_soft_wins_linear": bool(soft_over_linear != "" and soft_over_linear < 1.0),
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
        "ma_skipped": summary.get("ma_skipped", ""),
        "ma_sec": summary.get("ma_sec", ""),
        "ma_over_linear": summary.get("ma_over_linear", ""),
        "ma_log_delta": summary.get("ma_log_delta", ""),
        "ma_n_square": summary.get("ma_n_square", ""),
        "ma_padded_a": summary.get("ma_padded_a", ""),
        "ma_padded_b": summary.get("ma_padded_b", ""),
        "ma_wins_linear": summary.get("ma_wins_linear", ""),
        "ma_reason": summary.get("ma_reason", ""),
        "matselect2_skipped": summary.get("matselect2_skipped", ""),
        "matselect2_sec": summary.get("matselect2_sec", ""),
        "matselect2_over_linear": summary.get("matselect2_over_linear", ""),
        "matselect2_log_delta": summary.get("matselect2_log_delta", ""),
        "matselect2_iterations": summary.get("matselect2_iterations", ""),
        "matselect2_rep_pops": summary.get("matselect2_rep_pops", ""),
        "matselect2_base_pops": summary.get("matselect2_base_pops", ""),
        "matselect2_removed": summary.get("matselect2_removed", ""),
        "matselect2_max_active_rows": summary.get("matselect2_max_active_rows", ""),
        "matselect2_wins_linear": summary.get("matselect2_wins_linear", ""),
        "matselect2_reason": summary.get("matselect2_reason", ""),
        "matselect2_soft_skipped": summary.get("matselect2_soft_skipped", ""),
        "matselect2_soft_sec": summary.get("matselect2_soft_sec", ""),
        "matselect2_soft_over_linear": summary.get("matselect2_soft_over_linear", ""),
        "matselect2_soft_log_delta": summary.get("matselect2_soft_log_delta", ""),
        "matselect2_soft_iterations": summary.get("matselect2_soft_iterations", ""),
        "matselect2_soft_extracts": summary.get("matselect2_soft_extracts", ""),
        "matselect2_soft_inserted": summary.get("matselect2_soft_inserted", ""),
        "matselect2_soft_candidates": summary.get("matselect2_soft_candidates", ""),
        "matselect2_soft_corrupt_returns": summary.get("matselect2_soft_corrupt_returns", ""),
        "matselect2_soft_newly_corrupt": summary.get("matselect2_soft_newly_corrupt", ""),
        "matselect2_soft_removed": summary.get("matselect2_soft_removed", ""),
        "matselect2_soft_max_active_rows": summary.get("matselect2_soft_max_active_rows", ""),
        "matselect2_soft_wins_linear": summary.get("matselect2_soft_wins_linear", ""),
        "matselect2_soft_reason": summary.get("matselect2_soft_reason", ""),
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
    ma_comparable = [row for row in completed if row["ma_over_linear"] != "" and row["ma_skipped"] is False]
    ma_wins = [row for row in ma_comparable if row["ma_wins_linear"] is True]
    ma_ratios = [float(row["ma_over_linear"]) for row in ma_comparable]
    ma_correct = [row for row in ma_comparable if float(row["ma_log_delta"]) == 0.0]
    mat_comparable = [row for row in completed if row["matselect2_over_linear"] != "" and row["matselect2_skipped"] is False]
    mat_wins = [row for row in mat_comparable if row["matselect2_wins_linear"] is True]
    mat_ratios = [float(row["matselect2_over_linear"]) for row in mat_comparable]
    mat_correct = [row for row in mat_comparable if float(row["matselect2_log_delta"]) == 0.0]
    soft_comparable = [row for row in completed if row["matselect2_soft_over_linear"] != "" and row["matselect2_soft_skipped"] is False]
    soft_wins = [row for row in soft_comparable if row["matselect2_soft_wins_linear"] is True]
    soft_ratios = [float(row["matselect2_soft_over_linear"]) for row in soft_comparable]
    soft_correct = [row for row in soft_comparable if float(row["matselect2_soft_log_delta"]) == 0.0]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "comparable_cases": len(comparable),
        "range_pruning_wins": len(wins),
        "ma_comparable_cases": len(ma_comparable),
        "ma_wins": len(ma_wins),
        "ma_exact_log_matches": len(ma_correct),
        "matselect2_comparable_cases": len(mat_comparable),
        "matselect2_wins": len(mat_wins),
        "matselect2_exact_log_matches": len(mat_correct),
        "matselect2_soft_comparable_cases": len(soft_comparable),
        "matselect2_soft_wins": len(soft_wins),
        "matselect2_soft_exact_log_matches": len(soft_correct),
        "all_completed": len(completed) == len(rows),
        "min_block_over_linear": min(ratios) if ratios else None,
        "max_block_over_linear": max(ratios) if ratios else None,
        "mean_block_over_linear": sum(ratios) / len(ratios) if ratios else None,
        "min_ma_over_linear": min(ma_ratios) if ma_ratios else None,
        "max_ma_over_linear": max(ma_ratios) if ma_ratios else None,
        "mean_ma_over_linear": sum(ma_ratios) / len(ma_ratios) if ma_ratios else None,
        "min_matselect2_over_linear": min(mat_ratios) if mat_ratios else None,
        "max_matselect2_over_linear": max(mat_ratios) if mat_ratios else None,
        "mean_matselect2_over_linear": sum(mat_ratios) / len(mat_ratios) if mat_ratios else None,
        "min_matselect2_soft_over_linear": min(soft_ratios) if soft_ratios else None,
        "max_matselect2_soft_over_linear": max(soft_ratios) if soft_ratios else None,
        "mean_matselect2_soft_over_linear": sum(soft_ratios) / len(soft_ratios) if soft_ratios else None,
    }


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    rows = report["summary_rows"]
    with (out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    validation = next((run for run in report.get("build_runs", []) if run.get("name") == "validate_ma_selector"), {})
    validation_metrics = validation.get("metrics", {})
    mat_validation = next((run for run in report.get("build_runs", []) if run.get("name") == "validate_matselect2_heap"), {})
    mat_validation_metrics = mat_validation.get("metrics", {})
    soft_validation = next((run for run in report.get("build_runs", []) if run.get("name") == "validate_matselect2_soft"), {})
    soft_validation_metrics = soft_validation.get("metrics", {})
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
        f"- Mirzaian-Arjomandi probe wins over linear saddleback count: `{report['aggregate']['ma_wins']}`",
        f"- Mean Mirzaian-Arjomandi/linear internal time ratio: `{report['aggregate']['mean_ma_over_linear']}`",
        f"- Mat-Select2 heap-primitive wins over linear saddleback count: `{report['aggregate']['matselect2_wins']}`",
        f"- Mean Mat-Select2 heap-primitive/linear internal time ratio: `{report['aggregate']['mean_matselect2_over_linear']}`",
        f"- Mat-Select2 soft-selector wins over linear saddleback count: `{report['aggregate']['matselect2_soft_wins']}`",
        f"- Mean Mat-Select2 soft-selector/linear internal time ratio: `{report['aggregate']['mean_matselect2_soft_over_linear']}`",
        f"- Mirzaian-Arjomandi exhaustive validation cases: `{validation_metrics.get('cases', '')}`",
        f"- Mirzaian-Arjomandi exhaustive validation failures: `{validation_metrics.get('failures', '')}`",
        f"- Mirzaian-Arjomandi exhaustive validation max delta: `{validation_metrics.get('max_delta', '')}`",
        f"- Mat-Select2 heap-primitive exhaustive validation cases: `{mat_validation_metrics.get('cases', '')}`",
        f"- Mat-Select2 heap-primitive exhaustive validation failures: `{mat_validation_metrics.get('failures', '')}`",
        f"- Mat-Select2 heap-primitive exhaustive validation max delta: `{mat_validation_metrics.get('max_delta', '')}`",
        f"- Mat-Select2 soft-selector exhaustive validation cases: `{soft_validation_metrics.get('cases', '')}`",
        f"- Mat-Select2 soft-selector exhaustive validation failures: `{soft_validation_metrics.get('failures', '')}`",
        f"- Mat-Select2 soft-selector exhaustive validation max delta: `{soft_validation_metrics.get('max_delta', '')}`",
        "",
        "The `block_rank` rows are sorted-matrix range-pruning probes, not a faithful",
        "Frederickson-Johnson implementation. The `ma_select_probe` row adapts the",
        "Mirzaian-Arjomandi square sorted-matrix selector by padding the shorter MITM",
        "side and selects only a log value. The `matselect2_heap_probe` row implements",
        "the exponential-block selector from the Kaplan/Frederickson-Johnson line of",
        "algorithms with an exact binary heap for the primitive that is soft-heap based",
        "in the asymptotic paper algorithm. The `matselect2_soft_probe` row uses the",
        "same exponential-block structure but replaces the row-list primitive with the",
        "validated soft-sequence-heap Soft-Select bridge; it remains a practical probe,",
        "not a proof of the published time bound.",
        "The `loh_topk_probe` row uses a capped",
        "output-style top-k rank (`N_probe`) and is not evidence for random access at",
        "the full rank when `N_probe != N`.",
        "",
        "| P | linear total s | best block s | block/linear | MA s | MA/linear | Mat2 heap s | Mat2 heap/linear | Mat2 soft s | Mat2 soft/linear | Mat2 soft delta | wall s | RSS KB | LOH probe N | LOH s |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        mat_sec = row["matselect2_sec"]
        mat_ratio = row["matselect2_over_linear"]
        mat_delta = row["matselect2_log_delta"]
        mat_sec_text = f"{float(mat_sec):.6f}" if mat_sec != "" else ""
        mat_ratio_text = f"{float(mat_ratio):.6f}" if mat_ratio != "" else ""
        soft_sec = row["matselect2_soft_sec"]
        soft_ratio = row["matselect2_soft_over_linear"]
        soft_delta = row["matselect2_soft_log_delta"]
        soft_sec_text = f"{float(soft_sec):.6f}" if soft_sec != "" else ""
        soft_ratio_text = f"{float(soft_ratio):.6f}" if soft_ratio != "" else ""
        soft_delta_text = f"{float(soft_delta):.6g}" if soft_delta != "" else ""
        lines.append(
            f"| `{row['primes']}` | {float(row['linear_total']):.6f} | "
            f"{float(row['best_block_total']):.6f} | "
            f"{float(row['block_over_linear']):.6f} | "
            f"{float(row['ma_sec']):.6f} | "
            f"{float(row['ma_over_linear']):.6f} | "
            f"{mat_sec_text} | {mat_ratio_text} | "
            f"{soft_sec_text} | {soft_ratio_text} | {soft_delta_text} | "
            f"{float(row['wall_seconds']):.6f} | "
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
    validate = single.run_measured(
        "validate_ma_selector",
        ["bin/smooth_xplusy_fj_loh_workbench", "validate-ma"],
    )
    report["build_runs"].append(validate)
    if not result_ok(validate):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "ma_validation_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"MA validation failed; wrote partial report to {out_dir}")
        return 1
    validate_mat = single.run_measured(
        "validate_matselect2_heap",
        ["bin/smooth_xplusy_fj_loh_workbench", "validate-matselect2"],
    )
    report["build_runs"].append(validate_mat)
    if not result_ok(validate_mat):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "matselect2_validation_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"Mat-Select2 validation failed; wrote partial report to {out_dir}")
        return 1
    validate_soft = single.run_measured(
        "validate_matselect2_soft",
        ["bin/smooth_xplusy_fj_loh_workbench", "validate-matselect2-soft"],
    )
    report["build_runs"].append(validate_soft)
    if not result_ok(validate_soft):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "matselect2_soft_validation_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"Mat-Select2 soft validation failed; wrote partial report to {out_dir}")
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
