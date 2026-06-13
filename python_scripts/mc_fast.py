"""
mc_fast.py
==========
Numpy-accelerated Metropolis-Hastings sampler over integer partitions
of L with weight  w_tilde(lambda) = (k! / prod m_j!) * 2^k.

Replaces mc.py with the same algorithms but using the numpy partition
representation from partition_numpy.py for O(sqrt(L)) vectorised
inner operations, giving ~10-50x speedup over the pure-Python version.

Can be run in two modes:

  1. Standalone:  python3 mc_fast.py
     Runs validation, run-length distributions, and scaling study
     sequentially (same output as mc.py).

  2. SLURM job-array task:  python3 mc_fast.py --task <i> --outdir <dir>
     Runs only the i-th L value from SCALING_L_VALS and writes a JSON
     result file to <dir>/scaling_<L>.json.
     Used by mc_slurm.sh.

See also: mc_uniform_fast.py  (uniform-measure sampler)
          aggregate.py         (combines SLURM job-array output)
"""

import argparse
import json
import math
import os
import random
import sys
import time
from collections import Counter

import numpy as np

from partition_numpy import (
    init_arrays, arrays_to_tuple,
    apply_transfer_np, apply_split_np, apply_merge_np,
    all_transfers_np, all_splits_np, all_merges_np,
    n_transfers_np, n_splits_np, n_merges_np,
    log_w_np,
    autocorrelation_np, tau_int_np,
)

# L values used by the scaling study (indexed by SLURM_ARRAY_TASK_ID)
SCALING_L_VALS = [20, 50, 100, 200, 300, 500]


# ============================================================
# MH step (weighted measure)
# ============================================================

def metropolis_step_fast(vals, mults, k, rng, p_transfer=0.5):
    """
    One MH step under the weighted measure w_tilde.
    Acceptance: min(1, exp(Delta_log_w) * N_fwd / N_rev).
    Returns (new_vals, new_mults, new_k, accepted).
    """
    u = rng.random()

    if u < p_transfer:
        moves = all_transfers_np(vals, mults)
        if not moves:
            return vals, mults, k, False
        a, b      = moves[rng.randrange(len(moves))]
        nv, nm, nk = apply_transfer_np(vals, mults, k, a, b)
        N_fwd     = len(moves)
        N_rev     = n_transfers_np(nv, nm)

    elif u < p_transfer + (1 - p_transfer) / 2:
        moves = all_splits_np(vals)
        if not moves:
            return vals, mults, k, False
        a, c      = moves[rng.randrange(len(moves))]
        nv, nm, nk = apply_split_np(vals, mults, k, a, c)
        N_fwd     = len(moves)
        N_rev     = n_merges_np(nv, nm)

    else:
        moves = all_merges_np(vals, mults)
        if not moves:
            return vals, mults, k, False
        c, d      = moves[rng.randrange(len(moves))]
        nv, nm, nk = apply_merge_np(vals, mults, k, c, d)
        N_fwd     = len(moves)
        N_rev     = n_splits_np(nv)

    if N_rev == 0:
        return vals, mults, k, False

    lw_old    = log_w_np(mults, k)
    lw_new    = log_w_np(nm, nk)
    log_alpha = (lw_new - lw_old) + math.log(N_fwd) - math.log(N_rev)

    if log_alpha >= 0 or rng.random() < math.exp(log_alpha):
        return nv, nm, nk, True
    return vals, mults, k, False


# ============================================================
# Validation against exact w_tilde
# ============================================================

def validate(L_max=14, n_mc=300000):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from core import partitions as ep
    except ImportError:
        print("core.py not found; skipping validation.")
        return

    print("Validation: MC vs exact w_tilde")
    print()

    for L in range(4, L_max + 1, 2):
        # Exact weighted distribution
        from partition_numpy import partition_to_arrays
        exact = {}
        Z = 0.0
        for p_t in ep(L):
            vs, ms = partition_to_arrays(list(p_t))
            kk = int(ms.sum())
            w  = math.exp(log_w_np(ms, kk))
            exact[p_t] = w
            Z += w
        for p in exact:
            exact[p] /= Z
        km_ex = sum(len(p) * exact[p] for p in exact)
        k2_ex = sum(len(p)**2 * exact[p] for p in exact)
        ks_ex = math.sqrt(max(k2_ex - km_ex**2, 0))

        rng  = random.Random(7)
        vals, mults = init_arrays(L)
        k = int(mults.sum())
        for _ in range(50000):
            vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)

        k_list = []
        visits = Counter()
        for _ in range(n_mc):
            vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
            k_list.append(k)
            visits[arrays_to_tuple(vals, mults)] += 1

        km_mc = sum(k_list) / n_mc
        ks_mc = math.sqrt(sum((kk - km_mc)**2 for kk in k_list) / n_mc)
        tv    = sum(abs(exact.get(p, 0) - visits.get(p, 0) / n_mc)
                    for p in set(exact) | set(visits)) / 2
        err   = abs(km_mc - km_ex) / km_ex * 100

        print(f"  L={L:3d}  <k>_ex={km_ex:.4f}  <k>_MC={km_mc:.4f}  "
              f"err={err:.2f}%  TV={tv:.4f}  "
              f"sig_ex={ks_ex:.4f}  sig_MC={ks_mc:.4f}")
    print()


# ============================================================
# Run-length distribution
# ============================================================

def run_length_distribution(L, n_steps=300000, seed=42):
    PI     = math.pi
    pred_k = PI / math.sqrt(6) * math.sqrt(L)
    pred_r = math.sqrt(6) / PI * math.sqrt(L)

    rng  = random.Random(seed)
    vals, mults = init_arrays(L)
    k = int(mults.sum())
    for _ in range(n_steps // 5):
        vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)

    dist  = Counter()
    k_sum = 0
    for _ in range(n_steps):
        vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
        k_sum += k
        for v, m in zip(vals.tolist(), mults.tolist()):
            dist[v] += m

    total  = sum(dist.values())
    k_mean = k_sum / n_steps
    r_bar  = L / k_mean
    q      = max(0, 1 - 1 / r_bar)

    print(f"  L={L}: <k>={k_mean:.2f} (pred {pred_k:.2f}), "
          f"r_bar={r_bar:.2f} (pred {pred_r:.2f})")
    print(f"  {'j':>4}  {'f_meas':>10}  {'f_pred':>10}  {'ratio':>7}")
    for j in range(1, min(max(dist) + 1, int(6 * r_bar) + 2)):
        fm = dist.get(j, 0) / total
        fp = (1 / r_bar) * q**(j - 1) if q > 0 else 0
        if fm < 5e-5 and fp < 5e-5:
            break
        print(f"  {j:>4}  {fm:>10.5f}  {fp:>10.5f}  "
              f"  {fm / fp if fp > 1e-9 else 0:>7.4f}")
    print()


# ============================================================
# Single-L scaling task (used standalone and by SLURM array)
# ============================================================

def run_scaling_task(L, n_per_L=None, seed=42):
    """
    Run the scaling study for a single L value.
    Returns a result dict suitable for JSON serialisation.
    """
    PI     = math.pi
    n      = n_per_L or max(100000, 500 * L)
    thin   = max(1, L // 20)
    burnin = max(20000, 200 * L)

    t0   = time.time()
    rng  = random.Random(seed)
    vals, mults = init_arrays(L)
    k = int(mults.sum())

    for _ in range(burnin):
        vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)

    k_list = []
    acc    = 0
    for i in range(n * thin):
        vals, mults, k, a = metropolis_step_fast(vals, mults, k, rng)
        acc += a
        if i % thin == 0:
            k_list.append(k)
    t1 = time.time()

    sqL  = math.sqrt(L)
    L14  = L**0.25
    km   = sum(k_list) / len(k_list)
    kvar = sum((kk - km)**2 for kk in k_list) / len(k_list)
    ksig = math.sqrt(max(kvar, 0))
    acf  = autocorrelation_np(k_list)
    ti   = tau_int_np(acf)
    ar   = acc / (n * thin)

    ratio_logL = km / (sqL * math.log(L)) if L > 1 else 0.0

    return {
        'L':          L,
        'km':         km,
        'km_over_sqL': km / sqL,
        'ratio_logL': ratio_logL,
        'ksig':       ksig,
        'ksig_over_L14': ksig / L14,
        'tau':        ti,
        'acc':        ar,
        'elapsed':    t1 - t0,
    }


# ============================================================
# Scaling study (sequential, for standalone run)
# ============================================================

def scaling_study(L_vals, n_per_L=None, seed=42):
    PI = math.pi
    print("Scaling study (weighted measure w_tilde = k! * 2^k / prod m_j!)")
    print("  NOTE: under this weighted measure E[k] ~ c * sqrt(L) * log(L),")
    print("  not sqrt(L). The ratio <k>/sqrt(L) grows as log(L) and does")
    print("  not converge. This is consistent with the canonical/GC")
    print("  inequivalence established in the paper (Section 6.3).")
    print()
    # Correct Erdos-Lehner coefficient: sqrt(6)/(2*pi)
    c_EL = math.sqrt(6) / (2 * PI)
    print(f"  Erdos-Lehner mean scaling: E[k] ~ {c_EL:.4f} * sqrt(L) * log(L)")
    print()
    hdr = (f"{'L':>6}  {'<k>/sqL':>9}  {'<k>/(sqL*lnL)':>14}  "
           f"{'sig/L^.25':>10}  {'tau':>7}  {'acc':>6}  {'t(s)':>5}")
    print(hdr)
    print("-" * len(hdr))

    results = []
    for L in L_vals:
        r = run_scaling_task(L, n_per_L=n_per_L, seed=seed)
        print(f"{r['L']:>6}  {r['km_over_sqL']:>9.4f}  "
              f"{r['ratio_logL']:>14.5f}  "
              f"{r['ksig_over_L14']:>10.4f}  "
              f"{r['tau']:>7.1f}  {r['acc']:>6.3f}  {r['elapsed']:>5.1f}")
        results.append(r)
    return results


# ============================================================
# Speed benchmark
# ============================================================

def benchmark(L_vals=None):
    if L_vals is None:
        L_vals = [100, 500, 1000, 5000]
    print("Speed benchmark (steps/sec):")
    for L in L_vals:
        rng = random.Random(42)
        vals, mults = init_arrays(L)
        k = int(mults.sum())
        for _ in range(1000):
            vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
        n  = 3000
        t0 = time.time()
        for _ in range(n):
            vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
        t1 = time.time()
        print(f"  L={L:>6}: {n / (t1 - t0):>9.0f} steps/sec  k~{k}")
    print()


# ============================================================
# Entry point
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="mc_fast.py weighted sampler")
    parser.add_argument('--task',   type=int, default=None,
                        help='SLURM array task index into SCALING_L_VALS')
    parser.add_argument('--outdir', type=str, default='results_mc',
                        help='Directory for JSON output (SLURM mode)')
    parser.add_argument('--seed',   type=int, default=42)
    args = parser.parse_args()

    if args.task is not None:
        # ---- SLURM job-array mode: run one L value, write JSON ----
        L = SCALING_L_VALS[args.task]
        print(f"[mc_fast] SLURM task {args.task}: L={L}", flush=True)
        os.makedirs(args.outdir, exist_ok=True)
        r = run_scaling_task(L, seed=args.seed)
        out_path = os.path.join(args.outdir, f"scaling_{L}.json")
        with open(out_path, 'w') as f:
            json.dump(r, f, indent=2)
        print(f"[mc_fast] Written {out_path}", flush=True)
        return

    # ---- Standalone sequential mode ----
    print("\n" + "=" * 62)
    print("Reptation partition MC sampler (numpy-accelerated)")
    print("=" * 62 + "\n")

    benchmark()
    validate(L_max=14, n_mc=300000)

    print("Run-length distributions")
    print()
    for L in [50, 100, 200]:
        run_length_distribution(L, n_steps=300000, seed=42)

    print("=" * 62)
    scaling_study(SCALING_L_VALS, seed=args.seed)


if __name__ == '__main__':
    main()
