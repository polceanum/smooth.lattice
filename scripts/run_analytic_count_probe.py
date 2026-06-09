#!/usr/bin/env python3
"""Run analytic lattice-count residual probes at certified target vectors."""
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
    ("k5_P2_3_5_7_11", "2,3,5,7,11", "1052,26,33,53,4", DEFAULT_N),
    ("k5_P2_3_5_7_13", "2,3,5,7,13", "205,279,119,131,16", DEFAULT_N),
    ("k5_P2_3_5_11_13", "2,3,5,11,13", "254,220,258,4,52", DEFAULT_N),
    ("k5_P2_3_7_11_13", "2,3,7,11,13", "291,331,90,84,28", DEFAULT_N),
    ("k5_P2_5_7_11_13", "2,5,7,11,13", "106,306,164,53,32", DEFAULT_N),
    ("k5_P3_5_7_11_13", "3,5,7,11,13", "273,219,5,219,5", DEFAULT_N),
    ("k6_first6", "2,3,5,7,11,13", "55,126,27,54,2,52", DEFAULT_N),
    ("k8_first8", "2,3,5,7,11,13,17,19", "75,28,9,16,3,22,5,1", DEFAULT_N),
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


def ok(run: dict[str, Any]) -> bool:
    return int(run.get("returncode", 1)) == 0


def case_from_arg(raw: str) -> tuple[str, str, str, int]:
    parts = raw.split(":")
    if len(parts) != 4:
        raise argparse.ArgumentTypeError("case must be label:primes_csv:exps_csv:N")
    return (parts[0], parts[1], parts[2], int(parts[3]))


def summarize_case(label: str, primes: str, exps: str, expected_n: int, run: dict[str, Any]) -> dict[str, Any]:
    metrics = run.get("metrics", {})
    return {
        "label": label,
        "primes": primes,
        "N": expected_n,
        "returncode": run.get("returncode"),
        "wall_seconds": run.get("elapsed_seconds", ""),
        "reported_seconds": metrics.get("seconds", ""),
        "build_seconds": metrics.get("build", ""),
        "count_phase_seconds": metrics.get("count_phase", ""),
        "k": metrics.get("k", ""),
        "T": metrics.get("T", ""),
        "analytic_count": metrics.get("analytic_count", ""),
        "analytic_derivative": metrics.get("analytic_derivative", ""),
        "layer_count": metrics.get("layer_count", ""),
        "expected_N": metrics.get("expected_N", expected_n),
        "residual_expected_minus_analytic": metrics.get("residual_expected_minus_analytic", ""),
        "residual_layer_minus_analytic": metrics.get("residual_layer_minus_analytic", ""),
        "residual_layer_minus_expected": metrics.get("residual_layer_minus_expected", ""),
        "relative_residual": metrics.get("relative_residual", ""),
        "A": metrics.get("A", ""),
        "Base": metrics.get("Base", ""),
        "splitA": metrics.get("splitA", ""),
        "splitBase": metrics.get("splitBase", ""),
        "outer": metrics.get("outer", ""),
        "digits": metrics.get("digits", ""),
        "exps": metrics.get("exps", exps),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["returncode"] == 0]
    abs_residuals = [
        abs(float(row["residual_expected_minus_analytic"]))
        for row in completed
        if row["residual_expected_minus_analytic"] != ""
    ]
    rel_residuals = [
        abs(float(row["relative_residual"]))
        for row in completed
        if row["relative_residual"] != ""
    ]
    layer_matches = [
        row for row in completed
        if row["residual_layer_minus_expected"] != "" and abs(float(row["residual_layer_minus_expected"])) < 0.5
    ]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "layer_count_matches_expected": len(layer_matches),
        "max_abs_expected_residual": max(abs_residuals) if abs_residuals else None,
        "mean_abs_expected_residual": sum(abs_residuals) / len(abs_residuals) if abs_residuals else None,
        "max_abs_relative_residual": max(rel_residuals) if rel_residuals else None,
        "mean_abs_relative_residual": sum(rel_residuals) / len(rel_residuals) if rel_residuals else None,
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
        "# Analytic Lattice Count Probe",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Cases: `{agg['cases']}`",
        f"- Completed cases: `{agg['completed_cases']}`",
        f"- Layer-count matches expected rank: `{agg['layer_count_matches_expected']}`",
        f"- Max |expected - analytic|: `{agg['max_abs_expected_residual']}`",
        f"- Mean |expected - analytic|: `{agg['mean_abs_expected_residual']}`",
        f"- Max absolute relative residual: `{agg['max_abs_relative_residual']}`",
        "",
        "These probes evaluate the asymptotic analytic count approximation at",
        "previously certified exponent vectors. `expected_N` is the certified rank",
        "for the vector; `layer_count` is the current floating-log layer counter at",
        "the same log height. This is diagnostic evidence only, not a proof of an",
        "analytic error bound.",
        "",
        "| label | k | residual expected-analytic | relative residual | derivative | layer-expected | wall s | A | Base |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['k']} | "
            f"{float(row['residual_expected_minus_analytic']):.6f} | "
            f"{float(row['relative_residual']):.9f} | "
            f"{float(row['analytic_derivative']):.6f} | "
            f"{float(row['residual_layer_minus_expected']):.0f} | "
            f"{float(row['wall_seconds']):.6f} | "
            f"{row['A']} | {row['Base']} |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=case_from_arg, dest="cases")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    args = parser.parse_args()

    harness = load_harness()
    cases = args.cases or DEFAULT_CASES
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"analytic_count_probe_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.analytic_count_probe.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "cases": [
                {"label": label, "primes": primes, "exps": exps, "N": expected_n}
                for label, primes, exps, expected_n in cases
            ],
        },
        "build_runs": [],
        "runs": [],
    }

    build = harness.run_measured("build_core", ["./build.sh"])
    report["build_runs"].append(build)
    if not ok(build):
        report["summary_rows"] = []
        report["aggregate"] = {"status": "build_failed"}
        out_dir.mkdir(parents=True, exist_ok=True)
        (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
        print(f"build failed; wrote partial report to {out_dir}")
        return 1

    rows = []
    for label, primes, exps, expected_n in cases:
        run = harness.run_measured(
            label,
            [
                "bin/smooth_layer_compressed_general",
                "count-probe-exps",
                primes,
                exps,
                str(expected_n),
            ],
        )
        report["runs"].append(run)
        rows.append(summarize_case(label, primes, exps, expected_n, run))

    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["completed_cases"] == report["aggregate"]["cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
