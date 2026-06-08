#!/usr/bin/env python3
"""Reproducible layer-compressed versus X+Y benchmark harness.

The default target is the first five primes at N=10^12:

    P = 2,3,5,7,11
    N = 1000000000000

The X+Y comparator currently measures adaptive Cartesian-sum value selection,
not exponent reconstruction. That makes the comparison conservative for the
layer-compressed solver, which returns the exponent vector.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
BIN = ROOT / "bin"
DEFAULT_PRIMES = "2,3,5,7,11"
DEFAULT_N = 1_000_000_000_000
KV_RE = re.compile(r"(\w+)=((?:\[[^\]]*\])|(?:\([^\)]*\))|(?:\"[^\"]*\")|(?:\S+))")


def utc_stamp() -> str:
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")


def run_text(args: list[str], *, check: bool = False) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        args,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=check,
    )


def optional_cmd(args: list[str]) -> str | None:
    try:
        proc = run_text(args)
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def parse_value(raw: str) -> Any:
    value = raw.strip().strip('"')
    if value in {"true", "false"}:
        return value == "true"
    if value.startswith("[") and value.endswith("]"):
        body = value[1:-1].strip()
        if not body:
            return []
        parts = [part.strip() for part in body.split(",")]
        try:
            return [int(part) for part in parts]
        except ValueError:
            return parts
    if value.startswith("(") and value.endswith(")"):
        body = value[1:-1].strip()
        if not body:
            return []
        parts = [part.strip() for part in body.split(",")]
        try:
            return [int(part) for part in parts]
        except ValueError:
            return parts
    try:
        return int(value)
    except ValueError:
        pass
    try:
        return float(value)
    except ValueError:
        return value


def parse_kv(stdout: str) -> dict[str, Any]:
    return {match.group(1): parse_value(match.group(2)) for match in KV_RE.finditer(stdout)}


def time_prefix(time_file: Path) -> list[str] | None:
    time_bin = Path("/usr/bin/time")
    if not time_bin.exists():
        return None
    system = platform.system()
    if system == "Darwin":
        return [str(time_bin), "-l", "-o", str(time_file)]
    if system == "Linux":
        return [str(time_bin), "-v", "-o", str(time_file)]
    return None


def parse_time_file(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {}
    text = path.read_text(errors="replace")
    out: dict[str, Any] = {"time_raw": text}
    for line in text.splitlines():
        stripped = line.strip()
        if "Maximum resident set size" in stripped:
            if ":" in stripped:
                _, value = stripped.rsplit(":", 1)
                try:
                    out["max_rss_kb"] = int(value.strip())
                except ValueError:
                    pass
            else:
                parts = stripped.split()
                if parts and parts[0].isdigit():
                    rss = int(parts[0])
                    # Darwin reports bytes for "maximum resident set size".
                    out["max_rss_kb"] = rss // 1024
        elif stripped.endswith("maximum resident set size"):
            parts = stripped.split()
            if parts and parts[0].isdigit():
                out["max_rss_kb"] = int(parts[0]) // 1024
    return out


def run_measured(name: str, command: list[str]) -> dict[str, Any]:
    with tempfile.TemporaryDirectory(prefix="smooth_bench_time_") as tmp:
        time_file = Path(tmp) / "time.txt"
        prefix = time_prefix(time_file)
        full_command = (prefix or []) + command
        start = time.perf_counter()
        proc = subprocess.run(
            full_command,
            cwd=ROOT,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
        elapsed = time.perf_counter() - start
        result = {
            "name": name,
            "command": command,
            "returncode": proc.returncode,
            "elapsed_seconds": elapsed,
            "stdout": proc.stdout.strip(),
            "stderr": proc.stderr.strip(),
            "metrics": parse_kv(proc.stdout),
        }
        result.update(parse_time_file(time_file))
        return result


def machine_metadata(cxx: str, cxxflags: str) -> dict[str, Any]:
    compiler_path = shutil.which(cxx) or cxx
    compiler_version = optional_cmd([cxx, "--version"])
    cpu_model = None
    mem_bytes = None
    if platform.system() == "Darwin":
        cpu_model = optional_cmd(["sysctl", "-n", "machdep.cpu.brand_string"])
        mem_raw = optional_cmd(["sysctl", "-n", "hw.memsize"])
        if mem_raw and mem_raw.isdigit():
            mem_bytes = int(mem_raw)
    elif platform.system() == "Linux":
        cpuinfo = Path("/proc/cpuinfo")
        meminfo = Path("/proc/meminfo")
        if cpuinfo.exists():
            for line in cpuinfo.read_text(errors="replace").splitlines():
                if line.startswith("model name"):
                    cpu_model = line.split(":", 1)[1].strip()
                    break
        if meminfo.exists():
            for line in meminfo.read_text(errors="replace").splitlines():
                if line.startswith("MemTotal:"):
                    parts = line.split()
                    if len(parts) >= 2 and parts[1].isdigit():
                        mem_bytes = int(parts[1]) * 1024
                    break
    return {
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "platform": platform.platform(),
        "system": platform.system(),
        "python": sys.version.split()[0],
        "git_commit": optional_cmd(["git", "rev-parse", "HEAD"]),
        "git_dirty": bool(optional_cmd(["git", "status", "--porcelain"])),
        "compiler": {
            "cxx": cxx,
            "path": compiler_path,
            "version": compiler_version,
            "cxxflags": cxxflags,
        },
        "cpu_model": cpu_model,
        "memory_bytes": mem_bytes,
    }


def build_baseline(cxx: str, cxxflags: str) -> dict[str, Any]:
    BIN.mkdir(exist_ok=True)
    return run_measured(
        "build_xplusy_baseline",
        [
            cxx,
            *cxxflags.split(),
            "benchmarks/smooth_xplusy_baseline.cpp",
            "-o",
            "bin/smooth_xplusy_baseline",
        ],
    )


def result_ok(result: dict[str, Any]) -> bool:
    return int(result.get("returncode", 1)) == 0


def exps_csv(exps: Any) -> str | None:
    if isinstance(exps, list) and all(isinstance(item, int) for item in exps):
        return ",".join(str(item) for item in exps)
    return None


def compute_comparison(runs: list[dict[str, Any]]) -> dict[str, Any]:
    by_name = {run["name"]: run for run in runs}
    layer = by_name.get("layer_compressed")
    xplusy = by_name.get("xplusy_adaptive")
    comparison: dict[str, Any] = {}
    if layer and xplusy and result_ok(layer) and result_ok(xplusy):
        layer_wall = float(layer["elapsed_seconds"])
        xplusy_wall = float(xplusy["elapsed_seconds"])
        if layer_wall > 0:
            comparison["xplusy_wall_over_layer_wall"] = xplusy_wall / layer_wall
        layer_sec = layer["metrics"].get("seconds")
        xplusy_metrics = xplusy["metrics"]
        xplusy_total = None
        if "build" in xplusy_metrics and "adaptive_seconds" in xplusy_metrics:
            xplusy_total = float(xplusy_metrics["build"]) + float(xplusy_metrics["adaptive_seconds"])
        if layer_sec and xplusy_total:
            comparison["xplusy_reported_total_over_layer_reported"] = xplusy_total / float(layer_sec)
        comparison["same_rank"] = layer["metrics"].get("N") == xplusy["metrics"].get("N")
        comparison["same_primes"] = layer["metrics"].get("P") == xplusy["metrics"].get("P")
    audit = by_name.get("interval_audit")
    if audit:
        comparison["rank_certified"] = bool(audit.get("metrics", {}).get("rank_certified", False))
        comparison["certified_count_le"] = audit.get("metrics", {}).get("certified_count_le")
    return comparison


def write_outputs(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")

    rows = []
    for run in report["runs"]:
        metrics = run.get("metrics", {})
        rows.append(
            {
                "name": run["name"],
                "returncode": run["returncode"],
                "elapsed_seconds": f"{run['elapsed_seconds']:.6f}",
                "max_rss_kb": run.get("max_rss_kb", ""),
                "method": metrics.get("method", run["name"]),
                "P": metrics.get("P", report["config"]["primes"]),
                "N": metrics.get("N", report["config"]["N"]),
                "reported_seconds": metrics.get("seconds", ""),
                "build_seconds": metrics.get("build", ""),
                "adaptive_seconds": metrics.get("adaptive_seconds", ""),
                "band": metrics.get("band", metrics.get("adaptive_band_n", "")),
                "exps": metrics.get("exps", ""),
                "rank_certified": metrics.get("rank_certified", ""),
                "certified_count_le": metrics.get("certified_count_le", ""),
            }
        )
    with (out_dir / "summary.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else ["name"])
        writer.writeheader()
        writer.writerows(rows)

    comparison = report.get("comparison", {})
    lines = [
        "# X+Y vs Layer-Compressed Benchmark",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        f"- Primes: `{report['config']['primes']}`",
        f"- N: `{report['config']['N']}`",
        f"- CXX: `{report['metadata']['compiler']['cxx']}`",
        f"- CXXFLAGS: `{report['metadata']['compiler']['cxxflags']}`",
        "",
        "The `xplusy_adaptive` row measures adaptive Cartesian-sum value selection only.",
        "The `layer_compressed` row measures full exponent-vector unranking.",
        "",
        "| method | return | wall s | max RSS KB | reported s | exps | certified |",
        "|---|---:|---:|---:|---:|---|---|",
    ]
    for row in rows:
        reported = row["reported_seconds"] or row["adaptive_seconds"]
        lines.append(
            f"| {row['name']} | {row['returncode']} | {row['elapsed_seconds']} | "
            f"{row['max_rss_kb']} | {reported} | `{row['exps']}` | {row['rank_certified']} |"
        )
    if comparison:
        lines.extend(["", "## Comparison", ""])
        for key, value in comparison.items():
            lines.append(f"- `{key}`: `{value}`")
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primes", default=DEFAULT_PRIMES)
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--out-dir", type=Path)
    parser.add_argument("--cxx", default=os.environ.get("CXX", "g++"))
    parser.add_argument("--cxxflags", default=os.environ.get("CXXFLAGS", "-O3 -std=c++17"))
    parser.add_argument("--target-gap", type=int, default=50_000)
    parser.add_argument("--heap-limit", type=int, default=0)
    parser.add_argument("--skip-audit", action="store_true")
    args = parser.parse_args()

    out_dir = args.out_dir or (ROOT / "results" / "local" / f"xplusy_vs_layer5_{utc_stamp()}")
    report: dict[str, Any] = {
        "schema": "smooth.lattice.xplusy_vs_layer5.v1",
        "metadata": machine_metadata(args.cxx, args.cxxflags),
        "config": {
            "primes": args.primes,
            "N": args.N,
            "target_gap": args.target_gap,
            "heap_limit": args.heap_limit,
            "skip_audit": args.skip_audit,
        },
        "runs": [],
    }

    build_core = run_measured("build_core", ["./build.sh"])
    report["runs"].append(build_core)
    if not result_ok(build_core):
        report["comparison"] = {"status": "build_failed"}
        write_outputs(out_dir, report)
        print(f"build failed; wrote partial report to {out_dir}", file=sys.stderr)
        return 1

    build_xy = build_baseline(args.cxx, args.cxxflags)
    report["runs"].append(build_xy)
    if not result_ok(build_xy):
        report["comparison"] = {"status": "xplusy_build_failed"}
        write_outputs(out_dir, report)
        print(f"X+Y baseline build failed; wrote partial report to {out_dir}", file=sys.stderr)
        return 1

    xplusy = run_measured(
        "xplusy_adaptive",
        ["bin/smooth_xplusy_baseline", "bench", args.primes, str(args.N), str(args.heap_limit)],
    )
    report["runs"].append(xplusy)

    layer = run_measured(
        "layer_compressed",
        [
            "bin/smooth_layer_compressed_general",
            "nth",
            args.primes,
            str(args.N),
            str(args.target_gap),
        ],
    )
    report["runs"].append(layer)

    layer_exps = exps_csv(layer.get("metrics", {}).get("exps"))
    if layer_exps and not args.skip_audit:
        audit = run_measured(
            "interval_audit",
            ["bin/smooth_interval_audit_exps_k6", args.primes, layer_exps, str(args.N)],
        )
        report["runs"].append(audit)

    report["comparison"] = compute_comparison(report["runs"])
    write_outputs(out_dir, report)
    print(out_dir)
    return 0 if all(result_ok(run) for run in report["runs"]) else 1


if __name__ == "__main__":
    raise SystemExit(main())
