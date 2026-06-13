"""
appendix_diagnostics.py
=======================
One-shot script to produce the two missing numbers for Appendix A of the
reptation-partitions paper, for the UNIFORM-measure sampler:

  (1) acceptance rate at each production L
  (2) burn-in convergence diagnostics (running mean, Gelman-Rubin R-hat,
      starting-point sensitivity) at L = 100, 200, 500

Place in the same directory as mc_uniform_fast.py and partition_numpy.py
(e.g. the mc_uniform/ production directory) and run:

    python3 appendix_diagnostics.py

Pure standard library + numpy; uses the existing sampler so the dynamics
are identical to the production runs. Runtime ~a few minutes total.
"""

import math
import random
import time

import numpy as np

from partition_numpy import init_arrays, partition_to_arrays
from mc_uniform_fast import metropolis_step_uniform_fast


def start_state(L, mode):
    """mode in {'low','default','high'} -> (vals, mults, k)."""
    if mode == "default":
        v, m = init_arrays(L)                 # ~1.28 sqrt(L) parts
    elif mode == "low":
        v, m = partition_to_arrays([L])       # single part, k=1 (under-dispersed)
    elif mode == "high":
        v, m = partition_to_arrays([1] * L)   # all ones, k=L (over-dispersed)
    else:
        raise ValueError(mode)
    return v, m, int(m.sum())


def acceptance_rates(L_vals, n_steps=200000, seed=42):
    print("=" * 60)
    print("(1) Acceptance rate (uniform sampler, production move mix)")
    print("=" * 60)
    print(f"  {'L':>6}  {'acc_rate':>9}  {'steps':>9}")
    for L in L_vals:
        rng = random.Random(seed)
        vals, mults, k = start_state(L, "default")
        for _ in range(max(20000, 100 * L)):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        acc = 0
        for _ in range(n_steps):
            vals, mults, k, a = metropolis_step_uniform_fast(vals, mults, k, rng)
            acc += a
        print(f"  {L:>6}  {acc / n_steps:>9.4f}  {n_steps:>9}")
    print()


def burnin_diagnostics(L, burnin=None, n_post=None, n_chains=4, seed=42):
    if burnin is None:
        burnin = 1000 * L
    if n_post is None:
        n_post = burnin

    print("=" * 60)
    print(f"(2) Burn-in diagnostics  L={L}  burnin={burnin}  chains={n_chains}")
    print("=" * 60)

    # running mean from the default start
    rng = random.Random(seed)
    vals, mults, k = start_state(L, "default")
    k0 = k
    ks = []
    for _ in range(burnin):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        ks.append(k)
    ks = np.array(ks)
    target = ks[len(ks) // 2:].mean()
    run = np.cumsum(ks) / np.arange(1, len(ks) + 1)
    band = 0.01 * target
    enter = next((i for i, r in enumerate(run) if abs(r - target) < band), None)
    print(f"  initial k:              {k0}")
    print(f"  2nd-half mean k:        {target:.2f}")
    if enter is not None:
        print(f"  running mean in +/-1% band at step: {enter} "
              f"(fraction {enter / burnin:.3f})")
    else:
        print(f"  running mean NOT in +/-1% band within burn-in")

    # Gelman-Rubin across n_chains, default start, different seeds
    chain_means, chain_vars = [], []
    for c in range(n_chains):
        rng = random.Random(1000 + c)
        vals, mults, k = start_state(L, "default")
        for _ in range(burnin):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        post = []
        for _ in range(n_post):                       # collect post-burn-in k
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
            post.append(k)
        post = np.array(post)
        chain_means.append(post.mean())
        chain_vars.append(post.var(ddof=1))
        print(f"  chain {c+1}: mean k = {post.mean():.2f}  std = {post.std():.2f}")
    m = burnin
    W = np.mean(chain_vars)
    B = m * np.var(chain_means, ddof=1)
    var_hat = (1 - 1 / m) * W + B / m
    Rhat = math.sqrt(var_hat / W) if W > 0 else float("nan")
    print(f"  W = {W:.3f}  B = {B:.3f}  R-hat = {Rhat:.5f}  "
          f"({'PASS' if Rhat < 1.01 else 'CHECK'})")

    # starting-point sensitivity
    finals = {}
    for mode in ("low", "default", "high"):
        rng = random.Random(seed)
        vals, mults, k = start_state(L, mode)
        for _ in range(burnin):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        tail = []
        for _ in range(burnin // 2):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
            tail.append(k)
        finals[mode] = np.mean(tail)
        print(f"  start {mode:>7} (k0={int(start_state(L,mode)[2])}): "
              f"2nd-half mean k = {finals[mode]:.2f}")
    spread = max(finals.values()) - min(finals.values())
    rel = spread / np.mean(list(finals.values())) * 100
    print(f"  max spread: {spread:.2f} ({rel:.2f}%)  "
          f"({'PASS' if rel < 2 else 'CHECK'})")
    print()


def main():
    t0 = time.time()
    for L in (200, 500):
        burnin_diagnostics(L)            # burnin = n_post = 1000*L
    print(f"Total elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
