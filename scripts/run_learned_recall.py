#!/usr/bin/env python3
"""Opt-in CPU experiment; requires PyTorch and NumPy, not toolkit dependencies.

Reference memory-block recurrences, not complete Kimi/GDN-2 architectures.
See --help and the saved protocol for all matching and task limitations.
"""
from __future__ import annotations

import argparse
import copy
from dataclasses import asdict, dataclass
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import platform
import random
import statistics
import sys
import time

import numpy as np
import torch
from torch import nn
from torch.nn import functional as F


ROOT = Path(__file__).resolve().parents[1]
MODELS = ("learned_kda", "structured_kda", "learned_gdn2", "structured_gdn2")


@dataclass(frozen=True)
class Config:
    key_dim: int = 8
    latent_dim: int = 2
    degree: int = 7
    value_dim: int = 16
    parameter_target: int = 32768
    batch_size: int = 24
    train_pairs: tuple[int, ...] = (4, 8, 16)

    @property
    def width(self):
        return (2 * self.degree + 1) ** self.latent_dim


def seed_all(seed):
    random.seed(seed)
    np.random.seed(seed)
    torch.manual_seed(seed)


def coefficient_square(theta, degree, dimension):
    # Full polynomial convolution, differentiating through both factors.
    shape = (1, 1) + (degree + 1,) * dimension
    kernel = theta.reshape(shape)
    convolution = {1: F.conv1d, 2: F.conv2d, 3: F.conv3d}[dimension]
    return convolution(kernel, kernel.flip(tuple(range(2, 2 + dimension))),
                       padding=degree).reshape(-1)


class Memory(nn.Module):
    def __init__(self, kind, config, hidden):
        super().__init__()
        self.kind, self.config, self.hidden = kind, config, hidden
        self.structured = kind.startswith("structured")
        self.decoupled = kind.endswith("gdn2")
        m, v = config.width, config.value_dim
        self.trunk = nn.Linear(config.key_dim, hidden)
        self.feature_head = nn.Linear(hidden, config.latent_dim if self.structured else m)
        gate_width = 2 * m + v if self.decoupled else m + 1
        self.gate_head = nn.Linear(hidden + v, gate_width)
        nn.init.zeros_(self.gate_head.weight)
        with torch.no_grad():
            self.gate_head.bias[:m].fill_(6.0)
            self.gate_head.bias[m:].fill_(2.0)
        self.log_scale = nn.Parameter(torch.tensor(math.log(5.0)))
        self.output_bias = nn.Parameter(torch.zeros(v))
        if self.structured:
            count = (config.degree + 1) ** config.latent_dim
            self.raw_theta = nn.Parameter(torch.full((count - 1,), math.log(math.expm1(1.0))))

    def coefficients(self):
        theta = torch.cat((self.raw_theta.new_ones(1), F.softplus(self.raw_theta) + 1e-6))
        return coefficient_square(theta, self.config.degree, self.config.latent_dim)

    def features(self, keys, coeffs=None):
        hidden = torch.tanh(self.trunk(keys))
        raw = self.feature_head(hidden)
        if self.structured:
            latent = torch.tanh(raw)
            coeffs = self.coefficients() if coeffs is None else coeffs
            # Compute a table of powers by multiplication, including at zero.
            powers = [torch.ones_like(latent)]
            for _ in range(2 * self.config.degree):
                powers.append(powers[-1] * latent)
            table = torch.stack(powers, dim=-1)
            raw = torch.ones((*keys.shape[:-1], self.config.width), dtype=keys.dtype, device=keys.device)
            raw = raw[..., :1]
            for axis in range(self.config.latent_dim):
                raw = (raw[..., :, None] * table[..., axis, None, :]).flatten(-2, -1)
            raw = raw * torch.sqrt(coeffs)
        else:
            # A final linear head would artificially restrict all queries to
            # hidden+1 dimensions. SiLU avoids that baseline rank bottleneck.
            raw = F.silu(raw)
        return F.normalize(raw, dim=-1, eps=1e-8), hidden

    def forward(self, keys, values, queries):
        coefficients = self.coefficients() if self.structured else None
        k, h = self.features(keys, coefficients)
        q, _ = self.features(queries, coefficients)
        gates = self.gate_head(torch.cat((h, values), dim=-1)).sigmoid()
        m, v = self.config.width, self.config.value_dim
        state = keys.new_zeros(keys.shape[0], m, v)
        for t in range(keys.shape[1]):
            decayed = gates[:, t, :m, None] * state
            key = k[:, t]
            if self.decoupled:
                erase = gates[:, t, m:2 * m] * key
                write = gates[:, t, 2 * m:] * values[:, t]
                read = (erase[:, :, None] * decayed).sum(dim=1)
                state = decayed + key[:, :, None] * (write - read)[:, None, :]
            else:
                read = (key[:, :, None] * decayed).sum(dim=1)
                correction = gates[:, t, m, None] * (values[:, t] - read)
                state = decayed + key[:, :, None] * correction[:, None, :]
        reads = torch.bmm(q, state)
        return self.log_scale.clamp(-4, 6).exp() * reads + self.output_bias


def make_model(kind, cfg, seed):
    # Every parameter participates in the model; no unused padding parameters.
    def count(hidden):
        m, v = cfg.width, cfg.value_dim
        out = cfg.latent_dim if kind.startswith("structured") else m
        gates = 2 * m + v if kind.endswith("gdn2") else m + 1
        theta = (cfg.degree + 1) ** cfg.latent_dim - 1 if kind.startswith("structured") else 0
        return (cfg.key_dim + 1) * hidden + (hidden + 1) * out + (hidden + v + 1) * gates + theta + 1 + v
    hidden = max(h for h in range(4, 513) if count(h) <= cfg.parameter_target)
    seed_all(seed)
    model = Memory(kind, cfg, hidden)
    assert sum(p.numel() for p in model.parameters()) == count(hidden)
    assert abs(count(hidden) - cfg.parameter_target) / cfg.parameter_target < 0.025
    return model


def episode(cfg, generator, pairs, overwrites=0, distractors=0, batch=None):
    """Explicit key/value writes and read-only queries; all keys are fresh.

    This isolates memory and learned feature encoding, not token parsing.
    Later overwrites must determine the answer, never the original value.
    """
    b = cfg.batch_size if batch is None else batch
    keys = torch.randn(b, pairs, cfg.key_dim, generator=generator)
    keys = F.normalize(keys, dim=-1)
    labels = torch.randint(cfg.value_dim, (b, pairs), generator=generator)
    latest = labels.clone()
    writes, values = [keys], [F.one_hot(labels, cfg.value_dim).float()]
    batch_index = torch.arange(b)
    for _ in range(overwrites):
        idx = torch.randint(pairs, (b,), generator=generator)
        val = torch.randint(cfg.value_dim, (b,), generator=generator)
        writes.append(keys[batch_index, idx, None])
        values.append(F.one_hot(val, cfg.value_dim).float()[:, None])
        latest[batch_index, idx] = val
    if distractors:
        extra = F.normalize(torch.randn(b, distractors, cfg.key_dim, generator=generator), dim=-1)
        val = torch.randint(cfg.value_dim, (b, distractors), generator=generator)
        writes.append(extra)
        values.append(F.one_hot(val, cfg.value_dim).float())
    order = torch.rand(b, pairs, generator=generator).argsort(dim=-1)
    queries = keys.gather(1, order[:, :, None].expand(-1, -1, cfg.key_dim))
    targets = latest.gather(1, order)
    return torch.cat(writes, dim=1), torch.cat(values, dim=1), queries, targets


CASES = {
    "recall_8": (8, 0, 0), "recall_16": (16, 0, 0),
    "recall_32": (32, 0, 0), "recall_64": (64, 0, 0),
    "overwrite_16": (16, 8, 0), "overwrite_32": (32, 16, 0),
    "distractors_16": (16, 0, 32),
}


@torch.no_grad()
def evaluate(model, cfg, batches=8, seed=91001):
    model.eval()
    result = {}
    for name, (pairs, overwrites, distractors) in CASES.items():
        generator = torch.Generator().manual_seed(seed + pairs + 101 * overwrites + 1009 * distractors)
        correct = total = 0
        loss = 0.0
        for _ in range(batches):
            k, v, q, y = episode(cfg, generator, pairs, overwrites, distractors, batch=32)
            logits = model(k, v, q)
            correct += (logits.argmax(-1) == y).sum().item()
            total += y.numel()
            loss += F.cross_entropy(logits.reshape(-1, cfg.value_dim), y.reshape(-1), reduction="sum").item()
        result[name] = {"accuracy": correct / total, "cross_entropy": loss / total,
                        "correct": correct, "query_count": total}
    model.train()
    return result


def primary_score(metrics):
    return statistics.mean(metrics[name]["accuracy"] for name in ("recall_16", "overwrite_16"))


def train(kind, cfg, seed, lr, seconds, *, fixed_steps=None, progress=None):
    model = make_model(kind, cfg, seed)
    optimizer = torch.optim.AdamW(model.parameters(), lr=lr, weight_decay=1e-4)
    generator = torch.Generator().manual_seed(100000 + seed)
    data_rng = random.Random(200000 + seed)
    cpu = wall = 0.0
    steps = query_count = writes = 0
    observations = []
    # Fixed exposures are a secondary control. CPU time is the primary budget.
    snapshots = {}
    common_steps = (100, 250, 500, 1000)
    last_loss = math.nan
    while (steps < fixed_steps if fixed_steps is not None else cpu < seconds):
        c0, w0 = time.process_time(), time.perf_counter()
        pairs = data_rng.choice(cfg.train_pairs)
        overwrites = pairs // 2 if data_rng.random() < 0.25 else 0
        distractors = pairs // 2 if data_rng.random() < 0.25 else 0
        k, v, q, y = episode(cfg, generator, pairs, overwrites, distractors)
        optimizer.zero_grad(set_to_none=True)
        logits = model(k, v, q)
        loss = F.cross_entropy(logits.reshape(-1, cfg.value_dim), y.reshape(-1))
        if not torch.isfinite(loss):
            raise RuntimeError(f"Nonfinite loss in {kind}, seed={seed}, step={steps}")
        loss.backward()
        grad_norm = nn.utils.clip_grad_norm_(model.parameters(), 1.0, error_if_nonfinite=True)
        optimizer.step()
        cpu += time.process_time() - c0
        wall += time.perf_counter() - w0
        steps += 1
        query_count += y.numel()
        writes += k.shape[0] * k.shape[1]
        last_loss = loss.item()
        if steps in common_steps:
            snapshots[str(steps)] = {k: t.detach().clone() for k, t in model.state_dict().items()}
        if steps % 100 == 0:
            row = {"step": steps, "training_cpu_seconds": cpu, "training_wall_seconds": wall,
                   "loss": last_loss, "gradient_norm_before_clip": grad_norm.item(),
                   "training_queries": query_count, "training_writes": writes}
            observations.append(row)
            if progress:
                progress(row)
    return model, {"kind": kind, "config": asdict(cfg), "seed": seed, "learning_rate": lr,
                   "parameters": sum(p.numel() for p in model.parameters()), "hidden_width": model.hidden,
                   "parameter_bytes_float32": 4 * sum(p.numel() for p in model.parameters()),
                   "persistent_buffer_bytes": sum(b.numel() * b.element_size() for b in model.buffers()),
                   "state_scalars": cfg.width * cfg.value_dim, "state_bytes_float32": 4 * cfg.width * cfg.value_dim,
                   "training_cpu_seconds": cpu, "training_wall_seconds": wall, "steps": steps,
                   "training_queries": query_count, "training_writes": writes, "final_batch_loss": last_loss,
                   "history": observations}, snapshots


@torch.no_grad()
def benchmark_episode(model, cfg):
    """Complete CPU reference forward, not a fused decode kernel benchmark."""
    result = {}
    for pairs in (16, 64):
        k, v, q, _ = episode(cfg, torch.Generator().manual_seed(33001), pairs, batch=1)
        for _ in range(3):
            model(k, v, q)
        cpu, wall = [], []
        for _ in range(15):
            c0, w0 = time.process_time(), time.perf_counter()
            model(k, v, q)
            cpu.append(1000 * (time.process_time() - c0))
            wall.append(1000 * (time.perf_counter() - w0))
        result[str(pairs)] = {"batch": 1, "writes": pairs, "queries": pairs, "repeats": 15,
                              "median_cpu_ms": statistics.median(cpu), "median_wall_ms": statistics.median(wall)}
    return result


def self_check():
    seed_all(7)
    cfg = Config(key_dim=4, latent_dim=1, degree=2, parameter_target=2048, batch_size=2)
    for kind in MODELS:
        model = Memory(kind, cfg, 7).double()
        k, v, q, _ = episode(cfg, torch.Generator().manual_seed(14), 3, 1, 1)
        k, v, q = k.double(), v.double(), q.double()
        prediction = model(k, v, q)
        # Independent explicit transition matrices, rather than rank-one loops.
        features, hidden = model.features(k)
        queries, _ = model.features(q)
        gates = model.gate_head(torch.cat((hidden, v), -1)).sigmoid()
        state = torch.zeros(2, cfg.width, cfg.value_dim, dtype=torch.double)
        identity = torch.eye(cfg.width, dtype=torch.double)[None].expand(2, -1, -1)
        for t in range(k.shape[1]):
            key = features[:, t]
            decay = torch.diag_embed(gates[:, t, :cfg.width])
            if model.decoupled:
                erase = gates[:, t, cfg.width:2 * cfg.width] * key
                value = gates[:, t, 2 * cfg.width:] * v[:, t]
                transition = identity - key[:, :, None] * erase[:, None, :]
            else:
                beta = gates[:, t, cfg.width]
                transition = identity - beta[:, None, None] * key[:, :, None] * key[:, None, :]
                value = beta[:, None] * v[:, t]
            state = torch.bmm(torch.bmm(transition, decay), state) + key[:, :, None] * value[:, None, :]
        expected = model.log_scale.exp() * torch.bmm(queries, state) + model.output_bias
        assert torch.allclose(prediction, expected, atol=1e-11, rtol=1e-11)
        # A full-network directional gradient check includes coefficients/gates.
        prediction.square().sum().backward()
        directions = {name: torch.randn_like(p) for name, p in model.named_parameters()}
        analytic = sum((p.grad * directions[name]).sum().item() for name, p in model.named_parameters())
        original = {name: p.detach().clone() for name, p in model.named_parameters()}
        endpoints = []
        with torch.no_grad():
            for sign in (-1, 1):
                for name, p in model.named_parameters():
                    p.copy_(original[name] + sign * 1e-6 * directions[name])
                endpoints.append(model(k, v, q).square().sum().item())
            for name, p in model.named_parameters():
                p.copy_(original[name])
        numeric = (endpoints[1] - endpoints[0]) / 2e-6
        assert math.isclose(analytic, numeric, rel_tol=2e-5, abs_tol=2e-6), (kind, analytic, numeric)
        if model.structured:
            phi, _ = model.features(k)
            latent = torch.tanh(model.feature_head(torch.tanh(model.trunk(k))))
            theta = torch.cat((latent.new_ones(1), F.softplus(model.raw_theta) + 1e-6))
            support = list(product(range(cfg.degree + 1), repeat=cfg.latent_dim))
            z = torch.stack([math.sqrt(1) * torch.prod(latent ** latent.new_tensor(a), dim=-1)
                             for a in support], dim=-1) * torch.sqrt(theta)
            z = F.normalize(z, dim=-1)
            assert torch.allclose(phi @ phi.transpose(-1, -2), (z @ z.transpose(-1, -2)).square(), atol=1e-11)
    # Latest-label oracle and positive/negative task controls.
    native = Memory("learned_kda", cfg, 2).double()
    probes = torch.randn(1, 9, cfg.key_dim, dtype=torch.double)
    native_features, _ = native.features(probes)
    assert torch.linalg.matrix_rank(native_features[0]).item() > native.hidden + 1
    c = Config(key_dim=8)
    k, v, q, labels = episode(c, torch.Generator().manual_seed(3), 16, 8, 8, batch=16)
    matches = (q[:, :, None, :] - k[:, None, :, :]).abs().amax(-1) < 1e-7
    positions = torch.arange(k.shape[1])[None, None, :].expand_as(matches)
    latest = positions.masked_fill(~matches, -1).amax(-1)
    recovered = v.argmax(-1).gather(1, latest)
    assert torch.equal(recovered, labels)
    # The fixed-state orthogonal-address control stores all 16 classes exactly.
    generator = torch.Generator().manual_seed(9)
    values = F.one_hot(torch.randint(16, (8, 16), generator=generator), 16).float()
    state = torch.zeros(8, 225, 16)
    for t in range(16):
        key = F.one_hot(torch.tensor(t), 225).float()
        state = state + key[None, :, None] * (values[:, t] - state[:, t])[:, None, :]
    assert torch.equal(state[:, :16].argmax(-1), values.argmax(-1))
    return {"independent_dense_recurrence_checks": 4, "full_network_directional_gradient_checks": 4,
            "structured_direct_kernel_checks": 2, "latest_overwrite_label_oracle": True,
            "learned_features_exceed_linear_hidden_rank_limit": True,
            "orthogonal_address_positive_control_accuracy": 1.0}


def save_json(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + "\n")


def protocol(args):
    return {"schema_version": 1,
            "command": (f"python3 scripts/run_learned_recall.py run --output {args.output} "
                        f"--seconds {args.seconds:g} --tune-seconds {args.tune_seconds:g} "
                        f"--key-dims {' '.join(map(str, args.key_dims))} "
                        f"--seeds {' '.join(map(str, args.seeds))} "
                        f"--learning-rates {' '.join(map(str, args.learning_rates))} "
                        f"--eval-batches {args.eval_batches}" +
                        (f" --ram-bytes {args.ram_bytes}" if args.ram_bytes is not None else "")),
            "scope": "CPU float32, one-layer explicit-write/read memory blocks. Native published recurrences with a common small encoder; not full Kimi Linear/GDN-2 models, official Zoology MQAR, fused kernels or GPU throughput.",
            "models": list(MODELS), "key_dimensions": args.key_dims, "seeds": args.seeds,
            "config": asdict(Config()), "learning_rates": args.learning_rates,
            "feature_encoders": "Shared query/key tanh trunk. Native: SiLU on the final M-feature projection followed by L2 normalization, avoiding a linear hidden-rank bottleneck. Structured: bounded tanh latent projection, learned positive polynomial coefficients, exact degree-two collected features and L2 normalization. No persistent index/basis buffers.",
            "preflight_correction": "Discarded a linear-head native baseline before any main test run; restarted all tuning with the nonlinear head. Final test data were not used to make that correction.",
            "tuning_seed": 101, "tuning_cpu_seconds_per_trial": args.tune_seconds,
            "training_cpu_seconds_per_run": args.seconds, "evaluation_batches_per_case": args.eval_batches,
            "matching": {"recurrent_state": "Exactly 225*16 float32 scalars = 14400 bytes per episode for all four blocks; encoder/gates/optimizer/activations reported separately.",
                         "parameters": "Common cap of 32768 active float32 parameters (128 KiB), largest fitting hidden width, within 2.5% of cap; no persistent basis/index buffers or unused padding. Optimizer state scales with actual parameter count.",
                         "training_cost": "Primary: equal process CPU seconds, one intra-op/inter-op thread. Includes synthetic batch generation, feature/coefficient preparation, all forward/backward, clipping and AdamW; excludes evaluation/checkpoint I/O. One final step may overshoot; actual times reported.",
                         "secondary": "Same-seed data-stream prefix at the largest predeclared common step milestone reached by every block. Equal examples, not a compute-equivalence claim. Peak training activation/RSS memory is not matched or claimed equal.",
                         "tuning": "Equal two-learning-rate CPU budget per architecture and key dimension, separate seed and validation keys; select on mean recall_16/overwrite_16 accuracy. Never select on final test data."},
            "task": {"keys": "Independent Gaussian vectors normalized to unit length, fresh each episode. Query keys equal stored keys. Values are independent uniform 16-way classes encoded one-hot.",
                     "training": "Uniform 4/8/16 target associations; independently 25% episodes add pairs/2 overwrites and 25% add pairs/2 unrelated writes. All target keys queried in shuffled order after writes; query events do not alter memory.",
                     "test": CASES, "chance_accuracy": 1 / 16,
                     "primary_test_score": "Unweighted mean of recall_16 and overwrite_16 accuracy, same fixed held-out episodes across models and seeds.",
                     "limitations": "Explicit episode boundaries, write/read masks, one-hot values and shared query/key encoder simplify the task. Synthetic diagnostic only; no natural-language generalization claim."},
            "decision_rule": {"advantage": "Require positive paired primary-score differences in every training seed and a positive descriptive paired 95% t interval; also inspect overwrite/longer-context behavior and measured end-to-end cost. Three seeds only supports a pilot conclusion.",
                              "negative": "If matched-state learned features win reliably, record no advantage for this tested structured family, preserve results, and do not justify custom kernels from state counts."},
            "primary_recurrence_sources": ["https://arxiv.org/html/2510.26692v2#S3", "https://arxiv.org/html/2605.22791v1#S3.SS1"],
            "runtime": {"python": platform.python_version(), "pytorch": torch.__version__, "numpy": np.__version__,
                        "system": platform.system(), "machine": platform.machine(), "device": "cpu", "threads": 1,
                        "compiler_flags": "Prebuilt PyTorch CPU wheel; no custom compiled kernel", "ram_bytes": args.ram_bytes},
            "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}


def run(args):
    output = ROOT / args.output
    checkpoint_dir = ROOT / "results/local/learned_recall_checkpoints" / output.name
    output.mkdir(parents=True, exist_ok=True)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    p = protocol(args)
    if (output / "protocol.json").exists():
        assert json.loads((output / "protocol.json").read_text()) == p, "Use a fresh output directory for a changed protocol."
    else:
        save_json(output / "protocol.json", p)
    save_json(output / "self_checks.json", self_check())
    tuning = {}
    for key_dim in args.key_dims:
        cfg = Config(key_dim=key_dim)
        for kind in MODELS:
            path = output / f"tune_d{key_dim}_{kind}.json"
            if path.exists():
                tuning[(key_dim, kind)] = json.loads(path.read_text())
                continue
            trials = []
            for lr in args.learning_rates:
                model, row, _ = train(kind, cfg, 101, lr, args.tune_seconds)
                metrics = evaluate(model, cfg, batches=2, seed=81001)
                row["validation"] = metrics
                row["validation_primary"] = primary_score(metrics)
                trials.append(row)
                print(json.dumps({"stage": "tuning", "key_dim": key_dim, "kind": kind,
                                  "lr": lr, "steps": row["steps"], "validation": row["validation_primary"]}), flush=True)
            best = max(trials, key=lambda row: row["validation_primary"])
            chosen = {"learning_rate": best["learning_rate"], "trials": trials}
            save_json(path, chosen)
            tuning[(key_dim, kind)] = chosen
    for key_dim in args.key_dims:
        cfg = Config(key_dim=key_dim)
        for seed in args.seeds:
            # Rotate architecture order to reduce systematic runtime drift.
            offset = args.seeds.index(seed) % len(MODELS)
            order = MODELS[offset:] + MODELS[:offset]
            for kind in order:
                path = output / f"run_d{key_dim}_{kind}_s{seed}.json"
                if path.exists():
                    continue
                lr = tuning[(key_dim, kind)]["learning_rate"]
                model, row, snapshots = train(kind, cfg, seed, lr, args.seconds)
                row["test"] = evaluate(model, cfg, batches=args.eval_batches)
                row["test_primary"] = primary_score(row["test"])
                row["complete_episode_cpu_reference_timing"] = benchmark_episode(model, cfg)
                row["equal_step_milestones"] = sorted(int(step) for step in snapshots)
                # Private checkpoints allow independent equal-exposure and timing audits.
                torch.save({"final": model.state_dict(), "snapshots": snapshots},
                           checkpoint_dir / f"d{key_dim}_{kind}_s{seed}.pt")
                save_json(path, row)
                print(json.dumps({"stage": "main", "key_dim": key_dim, "kind": kind, "seed": seed,
                                  "steps": row["steps"], "cpu_seconds": row["training_cpu_seconds"],
                                  "test_primary": row["test_primary"], "recall_64": row["test"]["recall_64"]["accuracy"]}), flush=True)
    summarize(args, output, checkpoint_dir)


def paired_stats(values):
    mean = statistics.mean(values)
    sd = statistics.stdev(values) if len(values) > 1 else 0.0
    # Predeclared three-seed pilot: two-sided 95% Student-t, df=2.
    critical = {3: 4.302652729911275}.get(len(values))
    half = critical * sd / math.sqrt(len(values)) if critical is not None else None
    return {"paired_differences": values, "mean": mean, "sample_sd": sd,
            "descriptive_95pct_t_interval": [mean - half, mean + half] if half is not None else None,
            "all_seeds_positive": all(x > 0 for x in values)}


def summarize(args, output, checkpoint_dir):
    rows = [json.loads(p.read_text()) for p in sorted(output.glob("run_*.json"))]
    summary = {"scope": "Predeclared small CPU synthetic memory-block comparison; no GPU or full-model claim.",
               "cells": [], "paired_comparisons": [], "common_step_controls": []}
    for key_dim in args.key_dims:
        cfg = Config(key_dim=key_dim)
        for kind in MODELS:
            cell = [r for r in rows if r["kind"] == kind and r["config"]["key_dim"] == key_dim]
            assert len(cell) == len(args.seeds)
            cases = {name: {"mean": statistics.mean(r["test"][name]["accuracy"] for r in cell),
                            "seed_values": [r["test"][name]["accuracy"] for r in cell]} for name in CASES}
            summary["cells"].append({"key_dim": key_dim, "kind": kind, "test": cases,
                                      "mean_primary": statistics.mean(r["test_primary"] for r in cell),
                                      "parameters": cell[0]["parameters"], "state_bytes_float32": cell[0]["state_bytes_float32"],
                                      "parameter_bytes_float32": cell[0]["parameter_bytes_float32"],
                                      "persistent_buffer_bytes": cell[0]["persistent_buffer_bytes"],
                                      "steps": [r["steps"] for r in cell],
                                      "mean_training_cpu_ms_per_step": statistics.mean(1000 * r["training_cpu_seconds"] / r["steps"] for r in cell),
                                      "median_episode_16_cpu_ms": statistics.median(r["complete_episode_cpu_reference_timing"]["16"]["median_cpu_ms"] for r in cell),
                                      "training_cpu_seconds": [r["training_cpu_seconds"] for r in cell]})
        for suffix in ("kda", "gdn2"):
            differences = []
            for seed in args.seeds:
                native = next(r for r in rows if r["kind"] == f"learned_{suffix}" and r["seed"] == seed and r["config"]["key_dim"] == key_dim)
                ours = next(r for r in rows if r["kind"] == f"structured_{suffix}" and r["seed"] == seed and r["config"]["key_dim"] == key_dim)
                differences.append(ours["test_primary"] - native["test_primary"])
            summary["paired_comparisons"].append({"key_dim": key_dim, "recurrence": suffix, **paired_stats(differences)})
        for seed in args.seeds:
            cell = [r for r in rows if r["seed"] == seed and r["config"]["key_dim"] == key_dim]
            common = set.intersection(*(set(r["equal_step_milestones"]) for r in cell))
            if not common:
                summary["common_step_controls"].append({"key_dim": key_dim, "seed": seed, "available": False})
                continue
            step = max(common)
            for kind in MODELS:
                lr = next(r["learning_rate"] for r in cell if r["kind"] == kind)
                model = make_model(kind, cfg, seed)
                states = torch.load(checkpoint_dir / f"d{key_dim}_{kind}_s{seed}.pt", map_location="cpu", weights_only=True)
                model.load_state_dict(states["snapshots"][str(step)])
                metrics = evaluate(model, cfg, batches=args.eval_batches)
                summary["common_step_controls"].append({"key_dim": key_dim, "seed": seed, "kind": kind,
                                                        "step": step, "learning_rate": lr,
                                                        "test_primary": primary_score(metrics), "test": metrics})
    save_json(output / "summary.json", summary)
    print(json.dumps({"stage": "complete", "paired_comparisons": summary["paired_comparisons"]}), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("self-check", "pilot", "run"))
    parser.add_argument("--output", default="results/benchmarks/learned_recall_v2")
    parser.add_argument("--key-dims", nargs="+", type=int, default=[2, 8])
    parser.add_argument("--seeds", nargs="+", type=int, default=[11, 23, 37])
    parser.add_argument("--learning-rates", nargs="+", type=float, default=[0.001, 0.003])
    parser.add_argument("--seconds", type=float, default=45.0)
    parser.add_argument("--tune-seconds", type=float, default=10.0)
    parser.add_argument("--eval-batches", type=int, default=8)
    parser.add_argument("--ram-bytes", type=int, default=None)
    args = parser.parse_args()
    torch.set_num_threads(1)
    torch.set_num_interop_threads(1)
    torch.use_deterministic_algorithms(True)
    if args.action == "self-check":
        print(json.dumps(self_check(), indent=2))
    elif args.action == "pilot":
        print(json.dumps(self_check(), indent=2))
        for kind in MODELS:
            cfg = Config()
            model, row, _ = train(kind, cfg, 101, 0.003, 5.0)
            row["validation"] = evaluate(model, cfg, batches=1, seed=81001)
            print(json.dumps({k: v for k, v in row.items() if k != "history"}), flush=True)
    else:
        run(args)


if __name__ == "__main__":
    main()
