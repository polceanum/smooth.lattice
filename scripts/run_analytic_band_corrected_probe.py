#!/usr/bin/env python3
"""Run residual-corrected analytic boundary-band probes."""
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
    ("k5_P2_3_5_7_11", "2,3,5,7,11", DEFAULT_N, "1052,26,33,53,4", 100, 200_000),
    ("k5_P2_3_5_7_13", "2,3,5,7,13", DEFAULT_N, "205,279,119,131,16", 100, 200_000),
    ("k5_P2_3_5_11_13", "2,3,5,11,13", DEFAULT_N, "254,220,258,4,52", 100, 200_000),
    ("k5_P2_3_7_11_13", "2,3,7,11,13", DEFAULT_N, "291,331,90,84,28", 100, 200_000),
    ("k5_P2_5_7_11_13", "2,5,7,11,13", DEFAULT_N, "106,306,164,53,32", 100, 200_000),
    ("k5_P3_5_7_11_13", "3,5,7,11,13", DEFAULT_N, "273,219,5,219,5", 100, 200_000),
    ("k6_first6", "2,3,5,7,11,13", DEFAULT_N, "55,126,27,54,2,52", 500, 200_000),
    ("k8_first8", "2,3,5,7,11,13,17,19", DEFAULT_N, "75,28,9,16,3,22,5,1", 10_000, 200_000),
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


def exps_csv(value: Any) -> str | None:
    if isinstance(value, list):
        return ",".join(str(part) for part in value)
    if isinstance(value, str) and value:
        stripped = value.strip()
        if stripped.startswith("[") and stripped.endswith("]"):
            stripped = stripped[1:-1]
        return stripped.replace(" ", "")
    return None


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
    audit: dict[str, Any] | None = None,
) -> dict[str, Any]:
    metrics = run.get("metrics", {})
    audit_metrics = (audit or {}).get("metrics", {})
    got_exps = metrics.get("exps", "")
    expected_list = parse_exps(expected_exps)
    center_count = metrics.get("center_count", "")
    center_error = ""
    if center_count != "":
        center_error = int(center_count) - n
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
        "raw_T": metrics.get("raw_T", ""),
        "raw_derivative": metrics.get("raw_derivative", ""),
        "center_count": center_count,
        "center_rank_error": center_error,
        "reported_center_rank_error": metrics.get("center_rank_error", ""),
        "correction": metrics.get("correction", ""),
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
        "audit_returncode": (audit or {}).get("returncode", ""),
        "audit_wall_seconds": (audit or {}).get("elapsed_seconds", ""),
        "rank_certified": bool(audit_metrics.get("rank_certified", False)),
        "certified_count_le": audit_metrics.get("certified_count_le", ""),
        "cert_seconds": audit_metrics.get("cert_seconds", ""),
        "groupA": audit_metrics.get("groupA", ""),
        "groupB": audit_metrics.get("groupB", ""),
        "ambiguous_possible": audit_metrics.get("ambiguous_possible", ""),
        "ambiguous_resolved": audit_metrics.get("ambiguous_resolved", ""),
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
    certified = [row for row in expected if row["rank_certified"] and row["certified_count_le"] == row["N"]]
    band_counts = [int(row["band_count"]) for row in completed if row["band_count"] != ""]
    return {
        "cases": len(rows),
        "completed_cases": len(completed),
        "target_inside_cases": len(inside),
        "enumerated_cases": len(enumerated),
        "recovered_cases": len(recovered),
        "recovered_expected_cases": len(expected),
        "rank_certified_cases": len(certified),
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
        "# Residual-Corrected Analytic Boundary Band Probe",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Cases: `{agg['cases']}`",
        f"- Completed cases: `{agg['completed_cases']}`",
        f"- Target inside band: `{agg['target_inside_cases']}`",
        f"- Enumerated cases: `{agg['enumerated_cases']}`",
        f"- Recovered expected vectors: `{agg['recovered_expected_cases']}`",
        f"- Independently rank-certified vectors: `{agg['rank_certified_cases']}`",
        f"- Max band count: `{agg['max_band_count']}`",
        f"- Mean band count: `{agg['mean_band_count']}`",
        "",
        "The band center starts at the analytic solution for rank N. The harness",
        "then performs one exact layer count at that center, shifts the center by",
        "`(N - count(center)) / analytic_derivative(center)`, and checks a smaller",
        "rank-radius band around the corrected center. Endpoint counts and band",
        "enumeration still use the floating-log layer counter. Each recovered vector",
        "is then passed to the independent interval-log rank auditor; only rows with",
        "`rank_certified=True` should be used as certified correctness evidence.",
        "",
        "| label | k | radius | center error | band count | inside | recovered | certified | wall s | audit s |",
        "|---|---:|---:|---:|---:|---|---|---|---:|---:|",
    ]
    for row in rows:
        audit_s = "" if row["audit_wall_seconds"] == "" else f"{float(row['audit_wall_seconds']):.6f}"
        lines.append(
            f"| `{row['label']}` | {row['k']} | {row['rank_radius']} | "
            f"{row['center_rank_error']} | {row['band_count']} | "
            f"{row['target_inside']} | {row['recovered_expected_exps']} | "
            f"{row['rank_certified']} | {float(row['wall_seconds']):.6f} | {audit_s} |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", action="append", type=case_from_arg, dest="cases")
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default="g++")
    parser.add_argument("--cxxflags", default="-O3 -std=c++17")
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    harness = load_harness()
    cases = args.cases or DEFAULT_CASES
    out_dir = args.out_dir or (ROOT / "results" / "local" / f"analytic_band_corrected_probe_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.analytic_band_corrected_probe.v1",
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
            "skip_audit": args.skip_audit,
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
                "analytic-band-corrected",
                primes,
                str(n),
                str(radius),
                str(max_candidates),
            ],
        )
        report["runs"].append(run)
        audit = None
        recovered_exps = exps_csv(run.get("metrics", {}).get("exps", ""))
        if recovered_exps is not None and not args.skip_audit:
            audit = harness.run_measured(
                f"{label}_interval_audit",
                ["bin/smooth_interval_audit_exps", primes, recovered_exps, str(n)],
            )
            report["runs"].append(audit)
        rows.append(summarize_case(label, primes, n, exps, radius, max_candidates, run, audit))

    report["summary_rows"] = rows
    report["aggregate"] = aggregate(rows)
    write_outputs(out_dir, report)
    print(out_dir)
    all_completed = report["aggregate"]["completed_cases"] == report["aggregate"]["cases"]
    all_certified = args.skip_audit or report["aggregate"]["rank_certified_cases"] == report["aggregate"]["cases"]
    return 0 if all_completed and all_certified else 1


if __name__ == "__main__":
    raise SystemExit(main())
