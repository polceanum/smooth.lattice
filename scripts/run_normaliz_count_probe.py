#!/usr/bin/env python3
"""Run bounded PyNormaliz lattice-count probes for smooth-number simplexes.

The target polytope is the rational simplex

    { e in R^k : e_i >= 0, sum_i a_i e_i <= B },

where ``a_i`` and ``B`` are integer rationalizations of ``log(p_i)`` and the
certified target value.  Normaliz receives this simplex through homogeneous
rational vertices ``0`` and ``(B/a_i) e_i``.

This is a comparator probe, not a correctness oracle for the irrational-log
problem.  The existing interval auditor remains the correctness certificate.
"""
from __future__ import annotations

import argparse
import csv
import json
import multiprocessing as mp
import subprocess
import time
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SCALE_DIGITS = (1, 2, 3)


def git_metadata() -> dict[str, Any]:
    def run(command: list[str]) -> str:
        proc = subprocess.run(command, cwd=ROOT, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return proc.stdout.strip() if proc.returncode == 0 else ""

    status = run(["git", "status", "--short"])
    return {
        "git_commit": run(["git", "rev-parse", "HEAD"]),
        "git_dirty": bool(status),
    }


def parse_primes(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part]


def rationalized_vertices(
    primes: list[int],
    exps: list[int],
    scale_digits: int,
    mode: str,
) -> tuple[list[list[int]], dict[str, Any]]:
    getcontext().prec = max(80, scale_digits + 30)
    scale = Decimal(10) ** scale_digits
    logs = [Decimal(p).ln() for p in primes]
    target = sum(Decimal(e) * logp for e, logp in zip(exps, logs))

    if mode == "inner":
        coeffs = [int((logp * scale).to_integral_value(rounding=ROUND_CEILING)) for logp in logs]
        bound = int((target * scale).to_integral_value(rounding=ROUND_FLOOR))
    elif mode == "outer":
        coeffs = [int((logp * scale).to_integral_value(rounding=ROUND_FLOOR)) for logp in logs]
        bound = int((target * scale).to_integral_value(rounding=ROUND_CEILING))
    elif mode == "rounded":
        coeffs = [int((logp * scale).to_integral_value()) for logp in logs]
        bound = int((target * scale).to_integral_value())
    else:
        raise ValueError(f"unknown mode: {mode}")

    vertices = [[0] * len(primes) + [1]]
    for i, coeff in enumerate(coeffs):
        vertex = [0] * len(primes) + [coeff]
        vertex[i] = bound
        vertices.append(vertex)

    meta = {
        "scale_digits": scale_digits,
        "scale": str(scale),
        "mode": mode,
        "coeffs": coeffs,
        "bound": bound,
        "vertex_format": "Normaliz rational vertices, coordinates followed by denominator",
    }
    return vertices, meta


def normaliz_count_worker(vertices: list[list[int]], queue: mp.Queue) -> None:
    started = time.perf_counter()
    try:
        from PyNormaliz import Cone

        cone = Cone(vertices=vertices)
        count = cone.NumberLatticePoints()
        queue.put(
            {
                "status": "ok",
                "count": int(count),
                "seconds": time.perf_counter() - started,
            }
        )
    except Exception as exc:  # pragma: no cover - exercised through external tool behavior
        queue.put(
            {
                "status": "error",
                "error_type": type(exc).__name__,
                "error": str(exc),
                "seconds": time.perf_counter() - started,
            }
        )


def timed_normaliz_count(vertices: list[list[int]], timeout_seconds: float) -> dict[str, Any]:
    ctx = mp.get_context("spawn")
    queue: mp.Queue = ctx.Queue()
    proc = ctx.Process(target=normaliz_count_worker, args=(vertices, queue))
    started = time.perf_counter()
    proc.start()
    proc.join(timeout_seconds)
    elapsed = time.perf_counter() - started
    if proc.is_alive():
        proc.terminate()
        proc.join(5)
        return {
            "status": "timeout",
            "seconds": elapsed,
            "timeout_seconds": timeout_seconds,
            "exitcode": proc.exitcode,
        }
    if not queue.empty():
        result = queue.get()
        result["exitcode"] = proc.exitcode
        return result
    return {
        "status": "error",
        "error_type": "NoWorkerResult",
        "error": "worker exited without putting a result on the queue",
        "seconds": elapsed,
        "exitcode": proc.exitcode,
    }


def toy_cases(timeout_seconds: float) -> list[dict[str, Any]]:
    cases = [
        {
            "name": "integer_2simplex_x_plus_y_le_3",
            "vertices": [[0, 0, 1], [3, 0, 1], [0, 3, 1]],
            "expected_count": 10,
        },
        {
            "name": "rational_2simplex_x_plus_y_le_3_over_2",
            "vertices": [[0, 0, 1], [3, 0, 2], [0, 3, 2]],
            "expected_count": 3,
        },
    ]
    rows = []
    for case in cases:
        result = timed_normaliz_count(case["vertices"], timeout_seconds)
        count = result.get("count")
        rows.append(
            {
                **case,
                **result,
                "passed": result.get("status") == "ok" and count == case["expected_count"],
            }
        )
    return rows


def target_cases(
    ma_report: Path,
    out_dir: Path,
    scale_digits: list[int],
    modes: list[str],
    timeout_seconds: float,
    max_cases: int | None,
) -> list[dict[str, Any]]:
    report = json.loads(ma_report.read_text())
    summary_rows = report.get("summary_rows", [])
    if max_cases is not None:
        summary_rows = summary_rows[:max_cases]

    vertices_dir = out_dir / "vertices"
    vertices_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, Any]] = []
    for row_index, source in enumerate(summary_rows):
        primes = parse_primes(source["primes"])
        exps = [int(x) for x in source["exps"]]
        expected_n = int(source["N"])
        slug = "P" + "_".join(str(p) for p in primes)
        for digits in scale_digits:
            for mode in modes:
                vertices, meta = rationalized_vertices(primes, exps, digits, mode)
                vertex_path = vertices_dir / f"{slug}_d{digits}_{mode}.json"
                vertex_path.write_text(json.dumps({"vertices": vertices, **meta}, indent=2, sort_keys=True) + "\n")
                result = timed_normaliz_count(vertices, timeout_seconds)
                count = result.get("count")
                delta = int(count) - expected_n if isinstance(count, int) else None
                rows.append(
                    {
                        "row_index": row_index,
                        "primes": source["primes"],
                        "k": len(primes),
                        "N": expected_n,
                        "exps": exps,
                        "scale_digits": digits,
                        "mode": mode,
                        "vertex_file": str(vertex_path.relative_to(out_dir)),
                        "coeffs": meta["coeffs"],
                        "bound": meta["bound"],
                        **result,
                        "delta_from_certified_N": delta,
                    }
                )
    return rows


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    target_rows = report["target_counts"]
    with (out_dir / "target_counts.csv").open("w", newline="") as f:
        fieldnames = [
            "row_index",
            "primes",
            "k",
            "N",
            "scale_digits",
            "mode",
            "status",
            "count",
            "delta_from_certified_N",
            "seconds",
            "timeout_seconds",
            "exitcode",
            "vertex_file",
            "error_type",
            "error",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(target_rows)

    aggregate = report["aggregate"]
    lines = [
        "# Normaliz Count Probe",
        "",
        f"- Tool: `PyNormaliz`",
        f"- Toy cases passed: `{aggregate['toy_passed']}/{aggregate['toy_cases']}`",
        f"- Target cases completed: `{aggregate['target_ok']}/{aggregate['target_cases']}`",
        f"- Target timeouts: `{aggregate['target_timeouts']}`",
        f"- Target errors: `{aggregate['target_errors']}`",
        f"- Timeout per count: `{report['config']['timeout_seconds']}` seconds",
        "",
        "The target counts use rationalized log simplexes and are comparator",
        "measurements, not independent certificates for the original irrational-log",
        "rank. Correctness certification still comes from the interval auditor.",
        "",
        "| P | d | mode | status | seconds | count | delta from certified N |",
        "|---|---:|---|---|---:|---:|---:|",
    ]
    for row in target_rows:
        seconds = row.get("seconds", "")
        seconds_text = f"{float(seconds):.6f}" if isinstance(seconds, (int, float)) else ""
        count = row.get("count", "")
        delta = row.get("delta_from_certified_N", "")
        lines.append(
            f"| `{row['primes']}` | {row['scale_digits']} | {row['mode']} | "
            f"{row['status']} | {seconds_text} | {count} | {delta} |"
        )
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ma-report", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--scale-digits", type=int, nargs="+", default=list(DEFAULT_SCALE_DIGITS))
    parser.add_argument("--mode", choices=("rounded", "inner", "outer"), nargs="+", default=["rounded"])
    parser.add_argument("--timeout-seconds", type=float, default=10.0)
    parser.add_argument("--max-cases", type=int)
    args = parser.parse_args()

    out_dir = args.out_dir if args.out_dir.is_absolute() else ROOT / args.out_dir
    ma_report = args.ma_report if args.ma_report.is_absolute() else ROOT / args.ma_report

    toys = toy_cases(args.timeout_seconds)
    targets = target_cases(ma_report, out_dir, args.scale_digits, args.mode, args.timeout_seconds, args.max_cases)
    aggregate = {
        "toy_cases": len(toys),
        "toy_passed": sum(1 for row in toys if row["passed"]),
        "target_cases": len(targets),
        "target_ok": sum(1 for row in targets if row.get("status") == "ok"),
        "target_timeouts": sum(1 for row in targets if row.get("status") == "timeout"),
        "target_errors": sum(1 for row in targets if row.get("status") == "error"),
    }
    report = {
        "schema": "smooth.lattice.normaliz_count_probe.v1",
        "metadata": git_metadata(),
        "config": {
            "ma_report": str(ma_report.relative_to(ROOT) if ma_report.is_relative_to(ROOT) else ma_report),
            "scale_digits": args.scale_digits,
            "modes": args.mode,
            "timeout_seconds": args.timeout_seconds,
            "max_cases": args.max_cases,
        },
        "aggregate": aggregate,
        "toy_counts": toys,
        "target_counts": targets,
    }
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if aggregate["toy_passed"] == aggregate["toy_cases"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
