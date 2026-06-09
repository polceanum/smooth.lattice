#!/usr/bin/env python3
"""Run analytic boundary-correction band probes at certified target ranks."""
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
    ("k5_P2_3_5_7_11", "2,3,5,7,11", DEFAULT_N, "1052,26,33,53,4", 1_000, 200_000),
    ("k5_P2_3_5_7_13", "2,3,5,7,13", DEFAULT_N, "205,279,119,131,16", 1_000, 200_000),
    ("k5_P2_3_5_11_13", "2,3,5,11,13", DEFAULT_N, "254,220,258,4,52", 1_000, 200_000),
    ("k5_P2_3_7_11_13", "2,3,7,11,13", DEFAULT_N, "291,331,90,84,28", 1_000, 200_000),
    ("k5_P2_5_7_11_13", "2,5,7,11,13", DEFAULT_N, "106,306,164,53,32", 1_000, 200_000),
    ("k5_P3_5_7_11_13", "3,5,7,11,13", DEFAULT_N, "273,219,5,219,5", 1_000, 200_000),
    ("k6_first6", "2,3,5,7,11,13", DEFAULT_N, "55,126,27,54,2,52", 2_000, 200_000),
    ("k8_first8", "2,3,5,7,11,13,17,19", DEFAULT_N, "75,28,9,16,3,22,5,1", 2_500_000, 200_000),
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


def parse_exps(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part]


def case_from_arg(raw: str) -> tuple[str, str, int, str, int, int]:
    parts = raw.split(":")
    if len(parts) != 6:
        raise argparse.ArgumentTypeError("case must be label:primes_csv:N:exps_csv:rank_radius:max_candidates")
    return (parts[0], parts[1], int(parts[2]), parts[3], int(parts[4]), int(parts[5]))


def summarize_case(
    label: str,
    primes: str,
    n: int,
    expected_exps: str,
    rank_radius: int,
    max_candidates: int,
    run: dict[str, Any],
) -> dict[str, Any]:
    metrics = run.get("metrics", {})
    got_exps = metrics.get("exps", "")
    expected_list = parse_exps(expected_exps)
    return {
        "label": label,
        "primes": primes,
        "N": n,
        "expected_exps": expected_list,
        "rank_radius": rank_radius,
        "max_candidates": max_candidates,
        "returncode": run.get("returncode"),
        "wall_seconds": run.get("elapsed_seconds", ""),
        "reported_seconds": metrics.get("seconds", ""),
        "build_seconds": metrics.get("build", ""),
        "count_phase_seconds": metrics.get("count_phase", ""),
        "band_phase_seconds": metrics.get("band_phase", ""),
        "exact_seconds": metrics.get("exact", ""),
        "k": metrics.get("k", ""),
        "T": metrics.get("T", ""),
        "derivative": metrics.get("derivative", ""),
        "half_width": metrics.get("half_width", ""),
        "below": metrics.get("below", ""),
        "above": metrics.get("above", ""),
        "band_count": metrics.get("band_count", ""),
        "target_inside": bool(metrics.get("target_inside", False)),
        "enumerated": bool(metrics.get("enumerated", False)),
        "recovered": bool(metrics.get("recovered", False)),
        "recovered_expected_exps": got_exps == expected_list,
        "cands": metrics.get("cands", ""),
        "A": metrics.get("A", ""),
        "Base": metrics.get("Base", ""),
        "splitA": metrics.get("splitA", ""),
        "splitBase": metrics.get("splitBase", ""),
        "outer": metrics.get("outer", ""),
        "digits": metrics.get("digits", ""),
        "exps": got_exps,
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    completed = [row for row in rows if row["returncode"] == 0]
    inside = [row for row in completed if row["target_inside"]]
    enumerated = [row for row in completed if row["enumerated"]]
    recovered = [row for row in completed if row["recovered"]]
    expected = [row for row in recovered if row["recovered_expected_exps"]]
    band_counts = [int(row["band_count"]) for row in completed if row["band_count"] != ""]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "target_inside_cases": len(inside),
        "enumerated_cases": len(enumerated),
        "recovered_cases": len(recovered),
        "recovered_expected_cases": len(expected),
        "max_band_count": max(band_counts) if band_counts else None,
        "mean_band_count": sum(band_counts) / len(band_counts) if band_counts else None,
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
        "# Analytic Boundary Band Probe",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Cases: `{agg['cases']}`",
        f"- Completed cases: `{agg['completed_cases']}`",
        f"- Target inside band: `{agg['target_inside_cases']}`",
        f"- Enumerated cases: `{agg['enumerated_cases']}`",
        f"- Recovered expected vectors: `{agg['recovered_expected_cases']}`",
        f"- Max band count: `{agg['max_band_count']}`",
        f"- Mean band count: `{agg['mean_band_count']}`",
        "",
        "The band center is the analytic solution for rank N. The half-width is",
        "`rank_radius / analytic_derivative`; endpoints are checked with the current",
        "layer counter. Cases are enumerated only when the endpoint band count is at",
        "or below `max_candidates`. This is an experimental correction-band probe,",
        "not a certified analytic oracle.",
        "",
        "| label | k | radius | band count | inside | enumerated | recovered | wall s | cands |",
        "|---|---:|---:|---:|---|---|---|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['rank_radius']} | "
            f"{row['band_count']} | {row['target_inside']} | {row['enumerated']} | "
            f"{row['recovered_expected_exps']} | {float(row['wall_seconds']):.6f} | "
            f"{row['cands']} |"
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
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"analytic_band_probe_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.analytic_band_probe.v1",
        "metadata": harness.machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "cases": [
                {
                    "label": label,
                    "primes": primes,
                    "N": n,
                    "expected_exps": exps,
                    "rank_radius": radius,
                    "max_candidates": max_candidates,
                }
                for label, primes, n, exps, radius, max_candidates in cases
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
    for label, primes, n, exps, radius, max_candidates in cases:
        run = harness.run_measured(
            label,
            [
                "bin/smooth_layer_compressed_general",
                "analytic-band",
                primes,
                str(n),
                str(radius),
                str(max_candidates),
            ],
        )
        report["runs"].append(run)
        rows.append(summarize_case(label, primes, n, exps, radius, max_candidates, run))

    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if report["aggregate"]["completed_cases"] == report["aggregate"]["cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
