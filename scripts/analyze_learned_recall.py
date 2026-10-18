#!/usr/bin/env python3
"""Dependency-free audit of saved learned-recall trials and matching controls."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import random
import statistics
import importlib.util
from itertools import product

ROOT = Path(__file__).resolve().parents[1]


def verify_implementation(directory):
    """Optional PyTorch checks of the actual two-dimensional trained encoder."""
    import torch
    from torch.nn import functional as functional
    import sys
    spec = importlib.util.spec_from_file_location("recall_verification", ROOT / "scripts/run_learned_recall.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    fixtures = []
    for kind in ("structured_kda", "structured_gdn2"):
        cfg = module.Config()
        model = module.make_model(kind, cfg, 71).double()
        x = torch.randn(2, 5, cfg.key_dim, dtype=torch.double)
        phi, hidden = model.features(x)
        latent = model.feature_head(hidden).tanh()
        theta = torch.cat((latent.new_ones(1), functional.softplus(model.raw_theta) + 1e-6))
        support = list(product(range(cfg.degree + 1), repeat=cfg.latent_dim))
        monomials = torch.stack([torch.prod(latent ** latent.new_tensor(a), dim=-1) for a in support], dim=-1)
        expanded = functional.normalize(monomials * theta.sqrt(), dim=-1)
        direct = (expanded @ expanded.transpose(-1, -2)).square()
        collected = phi @ phi.transpose(-1, -2)
        assert torch.allclose(collected, direct, atol=5e-13, rtol=5e-12)
        weights = torch.randn_like(direct)
        active = list(model.trunk.parameters()) + list(model.feature_head.parameters()) + [model.raw_theta]
        grads_collected = torch.autograd.grad((collected * weights).sum(), active, retain_graph=True)
        grads_direct = torch.autograd.grad((direct * weights).sum(), active)
        for actual, expected in zip(grads_collected, grads_direct):
            assert torch.allclose(actual, expected, atol=5e-11, rtol=5e-9)
        k, v, q, _ = module.episode(cfg, torch.Generator().manual_seed(23), 3, 1, 1, batch=2)
        k, v, q = k.double(), v.double(), q.double()
        loss = model(k, v, q).square().sum()
        loss.backward()
        directions = {name: torch.randn_like(p) / math.sqrt(p.numel()) for name, p in model.named_parameters()}
        analytic = sum((p.grad * directions[name]).sum().item() for name, p in model.named_parameters())
        original = {name: p.detach().clone() for name, p in model.named_parameters()}
        endpoints = []
        with torch.no_grad():
            for sign in (-1, 1):
                for name, p in model.named_parameters():
                    p.copy_(original[name] + sign * 1e-6 * directions[name])
                endpoints.append(model(k, v, q).square().sum().item())
        numeric = (endpoints[1] - endpoints[0]) / 2e-6
        assert math.isclose(analytic, numeric, rel_tol=2e-5, abs_tol=2e-6)
        fixtures.append({"kind": kind, "latent_dimension": 2, "degree": 7, "feature_width": 225,
                         "direct_kernel_max_abs_error": (direct - collected).abs().max().item(),
                         "independent_kernel_parameter_gradients_passed": True,
                         "full_network_directional_gradient_passed": True})
    (directory / "implementation_checks.json").write_text(json.dumps({"actual_training_configuration_checked": True, "fixtures": fixtures}, indent=2) + "\n")


def training_counts(config, seed, steps):
    rng = random.Random(200000 + seed)
    queries = writes = 0
    for _ in range(steps):
        n = rng.choice(config["train_pairs"])
        overwrite = n // 2 if rng.random() < 0.25 else 0
        distract = n // 2 if rng.random() < 0.25 else 0
        queries += config["batch_size"] * n
        writes += config["batch_size"] * (n + overwrite + distract)
    return queries, writes


def analyze(directory):
    protocol = json.loads((directory / "protocol.json").read_text())
    summary = json.loads((directory / "summary.json").read_text())
    assert hashlib.sha256((ROOT / "scripts/run_learned_recall.py").read_bytes()).hexdigest() == protocol["script_sha256"]
    rows = [json.loads(p.read_text()) for p in sorted(directory.glob("run_*.json"))]
    expected = {(d, k, s) for d in protocol["key_dimensions"] for k in protocol["models"] for s in protocol["seeds"]}
    actual = {(r["config"]["key_dim"], r["kind"], r["seed"]) for r in rows}
    assert actual == expected and len(rows) == len(expected)
    for row in rows:
        assert row["state_scalars"] == 3600 and row["state_bytes_float32"] == 14400
        assert row["parameters"] <= protocol["config"]["parameter_target"]
        assert row["persistent_buffer_bytes"] == 0
        assert row["parameter_bytes_float32"] == 4 * row["parameters"]
        assert protocol["training_cpu_seconds_per_run"] <= row["training_cpu_seconds"] < protocol["training_cpu_seconds_per_run"] + 1
        assert (row["training_queries"], row["training_writes"]) == training_counts(row["config"], row["seed"], row["steps"])
        for name, metric in row["test"].items():
            n, _, _ = protocol["task"]["test"][name]
            assert metric["query_count"] == 32 * protocol["evaluation_batches_per_case"] * n
            assert 0 <= metric["correct"] <= metric["query_count"]
            assert metric["accuracy"] == metric["correct"] / metric["query_count"]
            assert math.isfinite(metric["cross_entropy"])
        assert row["test_primary"] == statistics.mean(row["test"][k]["accuracy"] for k in ("recall_16", "overwrite_16"))
    comparisons = []
    for comparison in summary["paired_comparisons"]:
        d, family = comparison["key_dim"], comparison["recurrence"]
        differences = []
        for seed in protocol["seeds"]:
            ours = next(r for r in rows if (r["config"]["key_dim"], r["kind"], r["seed"]) == (d, f"structured_{family}", seed))
            native = next(r for r in rows if (r["config"]["key_dim"], r["kind"], r["seed"]) == (d, f"learned_{family}", seed))
            differences.append(ours["test_primary"] - native["test_primary"])
        assert comparison["paired_differences"] == differences
        assert comparison["mean"] == statistics.mean(differences)
        interval = comparison["descriptive_95pct_t_interval"]
        if len(differences) == 3:
            half = 4.302652729911275 * statistics.stdev(differences) / math.sqrt(3)
            assert all(math.isclose(x, y, abs_tol=1e-15) for x, y in
                       zip(interval, (statistics.mean(differences) - half, statistics.mean(differences) + half)))
        passes = comparison["all_seeds_positive"] and interval is not None and interval[0] > 0
        comparisons.append({"key_dim": d, "recurrence": family, "structured_minus_native_percentage_points": 100 * comparison["mean"],
                            "descriptive_interval_percentage_points": [100 * x for x in interval] if interval else None,
                            "predeclared_pilot_superiority_gate_passed": passes})
    common = summary["common_step_controls"]
    common_comparisons = []
    for d in protocol["key_dimensions"]:
        for family in ("kda", "gdn2"):
            deltas = []
            steps = []
            for seed in protocol["seeds"]:
                subset = [r for r in common if r["key_dim"] == d and r["seed"] == seed]
                assert len(subset) == 4
                assert len({r["step"] for r in subset}) == 1
                ours = next(r for r in subset if r["kind"] == f"structured_{family}")
                native = next(r for r in subset if r["kind"] == f"learned_{family}")
                deltas.append(ours["test_primary"] - native["test_primary"])
                steps.append(ours["step"])
            common_comparisons.append({"key_dim": d, "recurrence": family,
                                       "seed_differences_percentage_points": [100 * x for x in deltas],
                                       "common_steps_by_seed": steps,
                                       "mean_difference_percentage_points": 100 * statistics.mean(deltas)})
    table = []
    for cell in summary["cells"]:
        cell_rows = [r for r in rows if r["config"]["key_dim"] == cell["key_dim"] and r["kind"] == cell["kind"]]
        assert cell["mean_primary"] == statistics.mean(r["test_primary"] for r in cell_rows)
        table.append({"key_dim": cell["key_dim"], "kind": cell["kind"],
                      "primary_accuracy_percent": 100 * cell["mean_primary"],
                      "recall_16_percent": 100 * cell["test"]["recall_16"]["mean"],
                      "overwrite_16_percent": 100 * cell["test"]["overwrite_16"]["mean"],
                      "recall_64_percent": 100 * cell["test"]["recall_64"]["mean"],
                      "distractors_16_percent": 100 * cell["test"]["distractors_16"]["mean"],
                      "train_cpu_ms_per_step": cell["mean_training_cpu_ms_per_step"],
                      "episode_16_cpu_ms": cell["median_episode_16_cpu_ms"],
                      "parameters": cell["parameters"], "state_bytes": cell["state_bytes_float32"]})
    report = {"schema_version": 1, "matching_audits_passed": True, "complete_main_runs": len(rows),
              "host_logical_cpu_count": os.cpu_count(),
              "analysis_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "decision": "At least one narrow pilot comparison passes; inspect every condition and equal-exposure control before any advantage claim." if any(r["predeclared_pilot_superiority_gate_passed"] for r in comparisons) else "No tested structured variant passes the predeclared pilot superiority gate.",
              "cost_matched_comparisons": comparisons, "equal_exposure_comparisons": common_comparisons,
              "memory": {"state_bytes_exact": 14400, "parameter_byte_cap": 131072,
                         "actual_parameter_byte_min": min(r["parameter_bytes_float32"] for r in rows),
                         "actual_parameter_byte_max": max(r["parameter_bytes_float32"] for r in rows),
                         "training_activation_peak_matched": False},
              "table": table,
              "statistics_limit": "Three training seeds, fixed shared test episodes, descriptive t intervals; no multiple-comparison-adjusted confirmatory significance or general architecture superiority claim.",
              "sources_sha256": {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in sorted(directory.glob("*.json")) if p.name != "analysis.json"}}
    (directory / "analysis.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    with (directory / "comparison.csv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(table[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(table)
    print(json.dumps({k: report[k] for k in ("matching_audits_passed", "complete_main_runs", "decision", "cost_matched_comparisons", "equal_exposure_comparisons")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", nargs="?", default="results/benchmarks/learned_recall_v2")
    parser.add_argument("--check-model", action="store_true", help="Also check actual 2D kernels and backward using optional PyTorch")
    args = parser.parse_args()
    if args.check_model:
        verify_implementation(ROOT / args.directory)
    analyze(ROOT / args.directory)
