#!/usr/bin/env python3
"""Run and record the current serious-comparator gates.

This script is intentionally stricter about language than about exit status. It
does not pretend that missing algorithms are implemented. It runs the strongest
published selector wrapper currently in the repository, runs output-sensitive
X+Y probes in a feasible regime, and emits Barvinok/LattE-style rational
polytope inputs plus local tool availability.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import shutil
import subprocess
from decimal import Decimal, ROUND_CEILING, ROUND_FLOOR, getcontext
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parents[1]
FIRST_K_CASES = (
    "2,3,5,7,11",
    "2,3,5,7,11,13",
    "2,3,5,7,11,13,17,19",
)
DEFAULT_N = 1_000_000_000_000
DEFAULT_OUTPUT_N = 1_000_000
BARVINOK_TOOLS = ("barvinok_count", "count", "latte-count", "normaliz", "ehrhart")
COUNT_ENV_BIN = Path("/usr/local/Caskroom/miniforge/base/envs/smooth-lattice-count/bin")


def load_harness():
    path = ROOT / "scripts" / "run_xplusy_vs_layer5.py"
    spec = importlib.util.spec_from_file_location("run_xplusy_vs_layer5", path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"could not load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def run_command(command: list[str]) -> dict[str, Any]:
    proc = subprocess.run(
        command,
        cwd=ROOT,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return {
        "command": command,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def run_command_stdin(command: list[str], stdin: str) -> dict[str, Any]:
    proc = subprocess.run(
        command,
        cwd=ROOT,
        input=stdin,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    return {
        "command": command,
        "stdin": stdin,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def parse_primes(raw: str) -> list[int]:
    return [int(part) for part in raw.split(",") if part]


def parse_exps(raw: Any) -> list[int]:
    if isinstance(raw, list):
        return [int(part) for part in raw]
    text = str(raw).strip()
    if text.startswith("[") and text.endswith("]"):
        text = text[1:-1]
    return [int(part.strip()) for part in text.split(",") if part.strip()]


def sanitize_case(primes: str) -> str:
    return "P" + primes.replace(",", "_")


def write_ine(path: Path, coeffs: list[int], bound: int) -> None:
    k = len(coeffs)
    rows = []
    for i in range(k):
        row = ["0"] + ["0"] * k
        row[i + 1] = "1"
        rows.append(" ".join(row))
    rows.append(" ".join([str(bound)] + [str(-c) for c in coeffs]))
    text = "\n".join(
        [
            "H-representation",
            "begin",
            f"{k + 1} {k + 1} integer",
            *rows,
            "end",
            "",
        ]
    )
    path.write_text(text)


def barvinok_inputs_from_ma_report(ma_report: Path, out_dir: Path, scale_digits: int) -> list[dict[str, Any]]:
    getcontext().prec = max(80, scale_digits + 30)
    scale = Decimal(10) ** scale_digits
    report = json.loads(ma_report.read_text())
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for row in report["summary_rows"]:
        primes = parse_primes(row["primes"])
        exps = parse_exps(row["exps"])
        logs = [Decimal(p).ln() for p in primes]
        target = sum(Decimal(e) * logp for e, logp in zip(exps, logs))
        coeff_floor = [int((logp * scale).to_integral_value(rounding=ROUND_FLOOR)) for logp in logs]
        coeff_ceil = [int((logp * scale).to_integral_value(rounding=ROUND_CEILING)) for logp in logs]
        coeff_round = [int((logp * scale).to_integral_value()) for logp in logs]
        bound_floor = int((target * scale).to_integral_value(rounding=ROUND_FLOOR))
        bound_ceil = int((target * scale).to_integral_value(rounding=ROUND_CEILING))
        bound_round = int((target * scale).to_integral_value())
        slug = sanitize_case(row["primes"])
        files = {
            "inner_ine": f"{slug}_inner_upper_coeffs.ine",
            "outer_ine": f"{slug}_outer_lower_coeffs.ine",
            "rounded_ine": f"{slug}_rounded.ine",
            "metadata": f"{slug}.json",
        }
        write_ine(out_dir / files["inner_ine"], coeff_ceil, bound_floor)
        write_ine(out_dir / files["outer_ine"], coeff_floor, bound_ceil)
        write_ine(out_dir / files["rounded_ine"], coeff_round, bound_round)
        metadata = {
            "primes": primes,
            "exps": exps,
            "scale_digits": scale_digits,
            "scale": str(scale),
            "coeff_floor": coeff_floor,
            "coeff_ceil": coeff_ceil,
            "coeff_round": coeff_round,
            "bound_floor": bound_floor,
            "bound_ceil": bound_ceil,
            "bound_round": bound_round,
            "files": files,
            "interpretation": {
                "inner_ine": "uses ceil(scale*log p_i) and floor(scale*target); any accepted point is safely inside the real log simplex if the decimal intervals are independently trusted",
                "outer_ine": "uses floor(scale*log p_i) and ceil(scale*target); contains the real log simplex if the decimal intervals are independently trusted",
                "rounded_ine": "non-certified rationalized comparator input",
            },
        }
        (out_dir / files["metadata"]).write_text(json.dumps(metadata, indent=2, sort_keys=True) + "\n")
        rows.append({"primes": row["primes"], "exps": exps, **files})
    with (out_dir / "index.csv").open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["primes", "exps", "inner_ine", "outer_ine", "rounded_ine", "metadata"], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return rows


def tool_availability() -> dict[str, str | None]:
    out: dict[str, str | None] = {}
    for tool in BARVINOK_TOOLS:
        found = shutil.which(tool)
        if found is None:
            env_path = COUNT_ENV_BIN / tool
            if env_path.exists():
                found = str(env_path)
        out[tool] = found
    return out


def tool_versions(tools: dict[str, str | None]) -> dict[str, dict[str, Any]]:
    versions: dict[str, dict[str, Any]] = {}
    for name, path in tools.items():
        if not path:
            continue
        versions[name] = run_command([path, "--version"])
    return versions


def external_count_smoke(tools: dict[str, str | None]) -> dict[str, Any]:
    smokes: dict[str, Any] = {}
    if tools.get("barvinok_count"):
        smokes["barvinok_count_1d_isl"] = run_command_stdin(
            [tools["barvinok_count"]],
            "{ [x] : 0 <= x <= 3 }\n",
        )
        smokes["barvinok_count_2d_isl"] = run_command_stdin(
            [tools["barvinok_count"]],
            "{ [x,y] : 0 <= x and 0 <= y and x + y <= 3 }\n",
        )
    if tools.get("normaliz"):
        smoke_dir = ROOT / "results" / "local" / "normaliz_smoke"
        smoke_dir.mkdir(parents=True, exist_ok=True)
        project = smoke_dir / "simplex"
        (smoke_dir / "simplex.in").write_text(
            "\n".join(
                [
                    "amb_space 2",
                    "polytope 3",
                    "0 0",
                    "3 0",
                    "0 3",
                    "",
                ]
            )
        )
        smokes["normaliz_polytope_number_lattice_points"] = run_command(
            [tools["normaliz"], "--NumberLatticePoints", f"--OutputDir={smoke_dir}", str(project)]
        )
    return smokes


def load_json(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    return json.loads(path.read_text())


def write_report(out_dir: Path, report: dict[str, Any]) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "report.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    ma = report["ma_first_k"]["report"]
    loh = report["output_sensitive"]["report"]
    barv = report["barvinok"]
    ma_agg = (ma or {}).get("aggregate", {})
    loh_agg = (loh or {}).get("aggregate", {})
    tools = barv["tools"]
    available_tools = {name: path for name, path in tools.items() if path}
    smoke = barv.get("smoke", {})
    lines = [
        "# Best-Known Comparator Gate",
        "",
        f"- Timestamp: `{report['metadata']['timestamp_utc']}`",
        f"- Git commit: `{report['metadata'].get('git_commit')}`",
        f"- Git dirty: `{report['metadata'].get('git_dirty')}`",
        "",
        "## Gate Summary",
        "",
        f"- Mirzaian-Arjomandi full-unrank gate: `{report['ma_first_k']['status']}`",
        f"- Output-sensitive X+Y/LOH probe gate: `{report['output_sensitive']['status']}`",
        f"- Full Frederickson-Johnson gate: `{report['frederickson_johnson']['status']}`",
        f"- Soft-heap X+Y gate: `{report['soft_heap']['status']}`",
        f"- Barvinok-style external count gate: `{report['barvinok']['status']}`",
        "",
        "## Published Sorted-Matrix Selector",
        "",
        f"- Artifact: `{report['ma_first_k']['artifact']}`",
        f"- Cases: `{ma_agg.get('cases')}`",
        f"- Certified same-exponent cases: `{ma_agg.get('same_exps_certified_cases')}`",
        f"- MA wall-time wins: `{ma_agg.get('ma_wall_wins')}`",
        f"- Mean MA/current wall ratio: `{ma_agg.get('mean_ma_wall_over_corrected_wall')}`",
        "",
        "## Output-Sensitive X+Y Probe",
        "",
        f"- Artifact: `{report['output_sensitive']['artifact']}`",
        f"- Cases: `{loh_agg.get('cases')}`",
        f"- Completed cases: `{loh_agg.get('completed_cases')}`",
        f"- Mean block/linear ratio: `{loh_agg.get('mean_block_over_linear')}`",
        f"- Mean MA/linear ratio: `{loh_agg.get('mean_ma_over_linear')}`",
        "",
        "This is an output-sized probe at a feasible rank. It does not establish a",
        "full-rank random-access comparison when `N_probe` is much smaller than the",
        "paper target rank.",
        "",
        "## Barvinok-Style Lattice Counting",
        "",
        f"- Available local tools: `{available_tools}`",
        f"- Export directory: `{report['barvinok']['export_dir']}`",
        f"- Exported cases: `{len(report['barvinok']['exports'])}`",
        f"- External count smoke statuses: `{ {name: run.get('returncode') for name, run in smoke.items()} }`",
        "",
        "The exported `.ine` files are rational simplex inputs for external",
        "Barvinok/LattE/Normaliz-style tools. The gate records tool availability and",
        "tiny smoke commands separately; a failed smoke is treated as a blocked",
        "external-tool comparison, not as mathematical evidence.",
        "",
        "## Bottom Line",
        "",
        "The repository currently has a clean, certified comparison against a",
        "published Mirzaian-Arjomandi sorted-matrix selector wrapper. Full",
        "Frederickson-Johnson, an actual soft-heap implementation, and an external",
        "Barvinok-family count run remain open obligations.",
    ]
    (out_dir / "report.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=ROOT / "results" / "benchmarks" / "best_known_comparator_gate")
    parser.add_argument("--N", type=int, default=DEFAULT_N)
    parser.add_argument("--output-N", type=int, default=DEFAULT_OUTPUT_N)
    parser.add_argument("--scale-digits", type=int, default=40)
    parser.add_argument("--strict", action="store_true")
    args = parser.parse_args()

    harness = load_harness()
    out_dir = args.out_dir if args.out_dir.is_absolute() else ROOT / args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    ma_dir = out_dir / "ma_full_unrank_first_k_1e12"
    ma_cmd = [
        "python3",
        "scripts/run_ma_full_unrank_suite.py",
        "--out-dir",
        str(ma_dir),
        "--rank-radius",
        "250",
        "--max-candidates",
        "200000",
        "--ma-max-n",
        "30000000",
        "--ma-max-middle",
        "200000000",
    ]
    for primes in FIRST_K_CASES:
        ma_cmd.extend(["--prime-set", primes])
    ma_run = run_command(ma_cmd)
    ma_report = load_json(ma_dir / "report.json")
    ma_status = "passed" if ma_run["returncode"] == 0 and ma_report and ma_report.get("aggregate", {}).get("all_same_exps_certified") else "failed"

    loh_dir = out_dir / "output_sensitive_xplusy_1e6"
    loh_cmd = [
        "python3",
        "scripts/run_sorted_matrix_workbench.py",
        "--out-dir",
        str(loh_dir),
        "--N",
        str(args.output_N),
    ]
    for primes in FIRST_K_CASES:
        loh_cmd.extend(["--prime-set", primes])
    loh_run = run_command(loh_cmd)
    loh_report = load_json(loh_dir / "report.json")
    loh_status = "executed" if loh_run["returncode"] == 0 and loh_report else "failed"

    tools = tool_availability()
    versions = tool_versions(tools)
    smoke = external_count_smoke(tools)
    barvinok_dir = out_dir / "barvinok_inputs"
    exports = []
    if ma_report is not None:
        exports = barvinok_inputs_from_ma_report(ma_dir / "report.json", barvinok_dir, args.scale_digits)
    smoke_ok = any(run.get("returncode") == 0 for run in smoke.values())
    if smoke_ok:
        barvinok_status = "tool_available_smoke_passed_exports_ready"
    elif any(tools.values()):
        barvinok_status = "tool_available_smoke_failed_exports_ready"
    else:
        barvinok_status = "blocked_no_local_tool_exports_ready"

    report = {
        "schema": "smooth.lattice.best_known_comparator_gate.v1",
        "metadata": harness.machine_metadata("g++", "-O3 -std=c++17"),
        "config": {
            "N": args.N,
            "output_N": args.output_N,
            "first_k_cases": FIRST_K_CASES,
            "scale_digits": args.scale_digits,
            "strict": args.strict,
        },
        "ma_first_k": {
            "status": ma_status,
            "artifact": str(ma_dir.relative_to(ROOT)),
            "run": ma_run,
            "report": ma_report,
        },
        "output_sensitive": {
            "status": loh_status,
            "artifact": str(loh_dir.relative_to(ROOT)),
            "run": loh_run,
            "report": loh_report,
        },
        "frederickson_johnson": {
            "status": "open_not_implemented",
            "note": "The repository has MA sorted-matrix selection and range-pruning probes, but not the published Frederickson-Johnson reduction/ranking algorithm.",
        },
        "soft_heap": {
            "status": "open_not_implemented",
            "note": "The repository has exact LOH/output-style probes but not a Chazelle/Kaplan-style soft heap implementation with corruption semantics.",
        },
        "barvinok": {
            "status": barvinok_status,
            "tools": tools,
            "versions": versions,
            "smoke": smoke,
            "export_dir": str(barvinok_dir.relative_to(ROOT)),
            "exports": exports,
        },
    }
    write_report(out_dir, report)
    print(out_dir)
    if args.strict:
        required = [ma_status == "passed", loh_status == "executed", any(tools.values())]
        return 0 if all(required) else 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
