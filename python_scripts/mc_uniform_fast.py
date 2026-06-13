"""
mc_uniform_fast.py
==================
Numpy-accelerated Metropolis-Hastings sampler over integer partitions
of L under the UNIFORM measure — every partition of L has equal weight.

Replaces mc_uniform.py with the same algorithms but using the numpy
partition representation from partition_numpy.py, giving ~10-50x
speedup over the pure-Python version.

Can be run in two modes:

  1. Standalone:  python3 mc_uniform_fast.py
     Runs validation, speed benchmark, mode/mean convergence study,
     and run-length distributions with BE vs geometric comparison.

  2. SLURM job-array task:
         python3 mc_uniform_fast.py --task <i> --outdir <dir>
     Runs either:
       i in 0..len(MODE_L_VALS)-1 : one L from MODE_L_VALS (mode study)
       i in len(MODE_L_VALS)..    : one L from RLD_L_VALS  (run-length)
     Writes a JSON result file to <dir>/.
     Used by mc_uniform_slurm.sh.

See also: mc_fast.py    (weighted sampler)
          aggregate.py  (combines SLURM job-array output)
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
    all_merges_np,
    all_splits_np,
    all_transfers_np,
    apply_merge_np,
    apply_split_np,
    apply_transfer_np,
    arrays_to_parts_list,
    arrays_to_tuple,
    autocorrelation_np,
    init_arrays,
    n_merges_np,
    n_splits_np,
    n_transfers_np,
    tau_int_np,
)

PI = math.pi

# L values for mode convergence study (indexed by SLURM_ARRAY_TASK_ID)
MODE_L_VALS = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000]
# L values for run-length distribution study
RLD_L_VALS = [100, 200, 500]


# ============================================================
# MH step (uniform measure)
# ============================================================


def metropolis_step_uniform_fast(vals, mults, k, rng, p_transfer=0.5):
    """
    One MH step under the uniform measure.
    Transfer uses exact enumeration (all_transfers_np) — no rejection loop.
    Split samples proportional to floor(a/2) as in mc_uniform.py.
    Merge samples uniformly over unordered pairs.
    Acceptance: min(1, N_fwd / N_rev).
    Returns (new_vals, new_mults, new_k, accepted).
    """
    u = rng.random()
    n_d = len(vals)

    if u < p_transfer:
        # --- Transfer: exact uniform sampling over valid (a,b) pairs ---
        moves = all_transfers_np(vals, mults)
        if not moves:
            return vals, mults, k, False
        a, b = moves[rng.randrange(len(moves))]
        nv, nm, nk = apply_transfer_np(vals, mults, k, a, b)
        N_fwd = n_transfers_np(vals, mults)
        N_rev = n_transfers_np(nv, nm)

    elif u < p_transfer + (1 - p_transfer) / 2:
        # --- Split: sample (a,c) proportional to floor(a/2) ---
        splittable = [(int(v), int(v) // 2) for v in vals if int(v) >= 2]
        if not splittable:
            return vals, mults, k, False
        total = sum(w for _, w in splittable)
        r = rng.randint(1, total)
        cumul = 0
        a = c = None
        for v, w in splittable:
            cumul += w
            if r <= cumul:
                a = v
                c = rng.randint(1, v // 2)
                break
        nv, nm, nk = apply_split_np(vals, mults, k, a, c)
        N_fwd = n_splits_np(vals)
        N_rev = n_merges_np(nv, nm)

    else:
        # --- Merge: sample uniformly over unordered pairs ---
        n_same = int((mults >= 2).sum())
        n_diff = n_d * (n_d - 1) // 2
        n_merge = n_same + n_diff
        if n_merge == 0:
            return vals, mults, k, False

        r = rng.randint(1, n_merge)
        if r <= n_diff:
            # distinct-value pair: sample via index trick
            i = rng.randint(0, n_d - 1)
            j = rng.randint(0, n_d - 2)
            if j >= i:
                j += 1
            c, d = int(vals[i]), int(vals[j])
        else:
            # same-value pair
            rep = vals[mults >= 2].tolist()
            v = int(rep[rng.randrange(len(rep))])
            c = d = v

        nv, nm, nk = apply_merge_np(vals, mults, k, min(c, d), max(c, d))
        N_fwd = n_merges_np(vals, mults)
        N_rev = n_splits_np(nv)

    if N_rev == 0:
        return vals, mults, k, False

    log_alpha = math.log(N_fwd) - math.log(N_rev)
    if log_alpha >= 0 or rng.random() < math.exp(log_alpha):
        return nv, nm, nk, True
    return vals, mults, k, False


# ============================================================
# Validation against exact uniform distribution
# ============================================================


def validate(L_max=20, n_mc=200000):
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    try:
        from core import partitions as ep
    except ImportError:
        print("core.py not found; skipping validation.")
        return

    print("Validation: MC vs exact uniform distribution")
    print()
    pred_k = PI / math.sqrt(6)

    print(
        f"  {'L':>4}  {'p(L)':>6}  {'modal_k_ex':>12}  "
        f"{'modal_k_HR':>12}  {'<k>_ex':>8}  {'<k>_MC':>8}  "
        f"{'err%':>6}  {'TV':>7}"
    )
    print(f"  {'-' * 78}")

    for L in range(4, L_max + 1, 2):
        all_p = list(ep(L))
        pL = len(all_p)
        k_cnt = Counter(len(p) for p in all_p)
        modal_ex = k_cnt.most_common(1)[0][0]
        km_ex = sum(len(p) for p in all_p) / pL
        sqL = math.sqrt(L)

        rng = random.Random(7)
        vals, mults = init_arrays(L)
        k = int(mults.sum())
        for _ in range(50000):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)

        visits = Counter()
        k_samp = []
        for _ in range(n_mc):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
            k_samp.append(k)
            visits[arrays_to_tuple(vals, mults)] += 1

        km_mc = sum(k_samp) / n_mc
        tv = sum(abs(1 / pL - visits.get(p, 0) / n_mc) for p in all_p) / 2
        err = abs(km_mc - km_ex) / km_ex * 100

        print(
            f"  {L:>4}  {pL:>6}  {modal_ex:>12}  "
            f"{pred_k * sqL:>12.2f}  {km_ex:>8.4f}  {km_mc:>8.4f}  "
            f"{err:>6.2f}%  {tv:>7.4f}"
        )

    print()
    print("Key: modal_k_ex is the exact most-probable k under uniform measure.")
    print("     modal_k_HR  is the Hardy-Ramanujan asymptotic prediction.")
    print("     These agree well even at small L.")
    print()


# ============================================================
# Single-L mode/mean convergence task
# ============================================================


def run_mode_task(L, n_per_L=None, seed=42):
    """
    Run the mode/mean convergence study for a single L value.
    Returns a result dict suitable for JSON serialisation.
    """
    n = n_per_L or max(1000000, 5000 * L)
    thin = max(1, L // 20)
    burnin = max(50000, 300 * L)

    t0 = time.time()
    rng = random.Random(seed)
    vals, mults = init_arrays(L)
    k = int(mults.sum())

    for _ in range(burnin):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)

    k_list = []
    k_cnt = Counter()
    for i in range(n * thin):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        if i % thin == 0:
            k_list.append(k)
            k_cnt[k] += 1
    t1 = time.time()

    sqL = math.sqrt(L)
    L14 = L**0.25
    km = sum(k_list) / len(k_list)
    kvar = sum((kk - km) ** 2 for kk in k_list) / len(k_list)
    ksig = math.sqrt(max(kvar, 0))
    acf = autocorrelation_np(k_list)
    ti = tau_int_np(acf)
    mk = k_cnt.most_common(1)[0][0]
    ratio_logL = mk / (sqL * math.log(L)) if L > 1 else 0.0

    return {
        "task": "mode",
        "L": L,
        "modal_k": mk,
        "modal_over_sqL": mk / sqL,
        "ratio_logL": ratio_logL,
        "km": km,
        "km_over_sqL": km / sqL,
        "ksig": ksig,
        "ksig_over_L14": ksig / L14,
        "tau": ti,
        "elapsed": t1 - t0,
    }


# ============================================================
# Mode convergence study (sequential, for standalone run)
# ============================================================


def mode_convergence_study(L_vals, n_per_L=None, seed=42):
    # Correct Erdos-Lehner coefficient: sqrt(6)/(2*pi)
    c_EL = math.sqrt(6) / (2 * PI)

    print("Mode convergence study")
    print(f"  Erdos-Lehner: E[k]_canon ~ {c_EL:.4f} * sqrt(L) * log(L)")
    print(f"  Modal k (HR saddle): modal_k ~ {PI / math.sqrt(6):.4f} * sqrt(L)")
    print(f"  Column '<k>/(sqL*lnL)' should converge to ~{c_EL:.4f}.")
    print()

    hdr = (
        f"{'L':>7}  {'modal_k':>8}  {'modal/sqL':>10}  "
        f"{'k*/(sqL*lnL)':>13}  {'<k>/sqL':>9}  {'sig/L^.25':>10}  "
        f"{'tau':>7}  {'t(s)':>5}"
    )
    print(hdr)
    print("-" * len(hdr))

    results = []
    for L in L_vals:
        r = run_mode_task(L, n_per_L=n_per_L, seed=seed)
        print(
            f"{r['L']:>7}  {r['modal_k']:>8}  "
            f"{r['modal_over_sqL']:>10.4f}  "
            f"{r['ratio_logL']:>13.5f}  "
            f"{r['km_over_sqL']:>9.4f}  "
            f"{r['ksig_over_L14']:>10.4f}  "
            f"{r['tau']:>7.1f}  {r['elapsed']:>5.1f}"
        )
        results.append(r)
    return results


# ============================================================
# BE distribution helper
# ============================================================


def _be_distribution(j_max, alpha):
    """
    Bose-Einstein run-length frequencies f_BE(j) for j=1..j_max,
    normalised over a sufficiently large tail.
    """
    raw = {}
    for j in range(1, j_max + 1):
        ea = math.exp(j * alpha)
        raw[j] = 1.0 / (ea - 1.0) if ea - 1.0 > 1e-300 else 0.0
    norm = sum(
        1.0 / (math.exp(j * alpha) - 1.0)
        for j in range(1, 10 * j_max + 1)
        if math.exp(j * alpha) - 1.0 > 1e-300
    )
    return {j: v / norm for j, v in raw.items()}


def _kl_divergence(p_dict, q_dict, j_vals):
    """KL(p || q) over j_vals, skipping j where p[j]=0."""
    kl = 0.0
    for j in j_vals:
        pj = p_dict.get(j, 0.0)
        qj = q_dict.get(j, 0.0)
        if pj > 0 and qj > 0:
            kl += pj * math.log(pj / qj)
    return kl


# ============================================================
# Single-L run-length distribution task
# ============================================================


def run_rld_task(L, n_steps=600000, seed=42):
    """
    Measure f(j) under the uniform measure and compare to BE and
    geometric predictions.  Returns a result dict for JSON output.
    """
    # Canonical mean from exact DP (Table 1 of paper)
    k_canon = {100: 21.75, 200: 34.17, 500: 61.37}[L]

    # Solve sum_{j>=1} 1/(exp(j*alpha)-1) = k_canon for alpha
    def gc_mean(a):
        s, j = 0.0, 1
        while True:
            t = 1.0 / (math.exp(j * a) - 1.0)
            s += t
            if t < 1e-12:
                break
            j += 1
        return s

    lo, hi = 1e-6, 10.0
    while gc_mean(lo) < k_canon:
        lo /= 2
    for _ in range(100):
        mid = (lo + hi) / 2
        if gc_mean(mid) > k_canon:
            lo = mid
        else:
            hi = mid
        if hi - lo < 1e-10:
            break
    alpha = (lo + hi) / 2
    pred_k_gc = PI / math.sqrt(6.0) * math.sqrt(L)
    pred_r_gc = math.sqrt(6.0) / PI * math.sqrt(L)

    rng = random.Random(seed)
    vals, mults = init_arrays(L)
    k = int(mults.sum())
    for _ in range(n_steps // 5):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)

    dist = Counter()
    k_sum = 0
    for _ in range(n_steps):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        k_sum += k
        for v, m in zip(vals.tolist(), mults.tolist()):
            dist[v] += m

    total = sum(dist.values())
    k_mean = k_sum / n_steps
    total = sum(dist.values())
    k_mean = k_sum / n_steps
    r_bar_mc = L / k_mean  # MC-measured mean run length
    r_bar = L / k_canon  # canonical mean run length (matches BE normalisation)
    q_geom = max(0.0, 1.0 - 1.0 / r_bar)

    j_max = max(dist.keys())
    j_range = range(1, j_max + 1)

    f_mc = {j: dist.get(j, 0) / total for j in j_range}
    f_geom = {j: (1.0 / r_bar) * q_geom ** (j - 1) for j in j_range}
    f_be = _be_distribution(j_max, alpha)

    support = [j for j in j_range if f_mc.get(j, 0) > 0]
    kl_be = _kl_divergence(f_mc, f_be, support)
    kl_geom = _kl_divergence(f_mc, f_geom, support)
    ratio = kl_geom / kl_be if kl_be > 1e-15 else float("inf")

    return {
        "task": "rld",
        "L": L,
        "k_mean": k_mean,
        "k_canon": k_canon,
        "r_bar_mc": r_bar_mc,
        "r_bar": r_bar,  # canonical, used for both BE and geom
        "alpha": alpha,
        "pred_k_gc": pred_k_gc,
        "pred_r_gc": pred_r_gc,
        "kl_be": kl_be,
        "kl_geom": kl_geom,
        "kl_ratio": ratio,
        "f_mc": {str(j): v for j, v in f_mc.items()},
        "f_be": {str(j): v for j, v in f_be.items()},
        "f_geom": {str(j): v for j, v in f_geom.items()},
    }


def print_rld_result(r):
    L = r["L"]
    print(
        f"  L={L}: <k>_MC={r['k_mean']:.2f}  <k>_canon={r['k_canon']:.2f}"
        f"  r_bar_MC={r['r_bar_mc']:.2f}  r_bar_canon={r['r_bar']:.2f}"
        f"  (GC pred: <k>={r['pred_k_gc']:.2f})"
    )
    print(
        f"  alpha (canonical) = {r['alpha']:.5f}"
        f"  [GC alpha = pi/sqrt(6L) = {math.pi / math.sqrt(6.0 * L):.5f}]"
    )
    print(
        f"  D_KL(MC||BE)={r['kl_be']:.4f}  "
        f"D_KL(MC||geom)={r['kl_geom']:.4f}  "
        f"ratio={r['kl_ratio']:.2f}"
    )
    print()
    print(
        f"  {'j':>4}  {'f_MC':>10}  {'f_BE':>10}  "
        f"{'ratio_BE':>9}  {'f_geom':>10}  {'ratio_geom':>11}"
    )
    f_mc = {int(j): v for j, v in r["f_mc"].items()}
    f_be = {int(j): v for j, v in r["f_be"].items()}
    f_geom = {int(j): v for j, v in r["f_geom"].items()}
    cutoff = int(6 * r["r_bar"]) + 2
    for j in range(1, cutoff):
        fm = f_mc.get(j, 0.0)
        fbe = f_be.get(j, 0.0)
        fg = f_geom.get(j, 0.0)
        if fm < 5e-5 and fbe < 5e-5:
            break
        r_be = fm / fbe if fbe > 1e-9 else 0.0
        r_geo = fm / fg if fg > 1e-9 else 0.0
        print(
            f"  {j:>4}  {fm:>10.5f}  {fbe:>10.5f}  "
            f"{r_be:>9.4f}  {fg:>10.5f}  {r_geo:>11.4f}"
        )
    print()


# ============================================================
# Speed benchmark
# ============================================================


def benchmark(L_vals=None):
    if L_vals is None:
        L_vals = [100, 500, 2000, 5000, 10000]
    print("Speed benchmark (steps/sec):")
    for L in L_vals:
        rng = random.Random(42)
        vals, mults = init_arrays(L)
        k = int(mults.sum())
        for _ in range(2000):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        n = 5000
        t0 = time.time()
        for _ in range(n):
            vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        t1 = time.time()
        print(f"  L={L:>6}: {n / (t1 - t0):>9.0f} steps/sec  k~{k}")
    print()


# ============================================================
# Entry point
# ============================================================


def main():
    parser = argparse.ArgumentParser(description="mc_uniform_fast.py uniform sampler")
    parser.add_argument("--task", type=int, default=None, help="SLURM array task index")
    parser.add_argument(
        "--outdir",
        type=str,
        default="results_mc_uniform",
        help="Directory for JSON output (SLURM mode)",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--n-steps",
        type=int,
        default=600000,
        help="Steps for run-length distribution tasks",
    )
    args = parser.parse_args()

    if args.task is not None:
        # ---- SLURM job-array mode ----
        os.makedirs(args.outdir, exist_ok=True)
        n_mode = len(MODE_L_VALS)

        if args.task < n_mode:
            L = MODE_L_VALS[args.task]
            print(f"[mc_uniform_fast] task {args.task}: mode study L={L}", flush=True)
            r = run_mode_task(L, seed=args.seed)
            out = os.path.join(args.outdir, f"mode_{L}.json")
        else:
            idx = args.task - n_mode
            L = RLD_L_VALS[idx]
            print(f"[mc_uniform_fast] task {args.task}: RLD L={L}", flush=True)
            r = run_rld_task(L, n_steps=args.n_steps, seed=args.seed)
            out = os.path.join(args.outdir, f"rld_{L}.json")

        with open(out, "w") as f:
            json.dump(r, f, indent=2)
        print(f"[mc_uniform_fast] Written {out}", flush=True)
        return

    # ---- Standalone sequential mode ----
    print()
    print("=" * 65)
    print("Uniform partition sampler: Hardy-Ramanujan verification")
    print("(numpy-accelerated)")
    print("=" * 65)
    print()

    validate(L_max=20, n_mc=200000)
    benchmark()

    print("=" * 65)
    mode_convergence_study(MODE_L_VALS, seed=args.seed)

    print()
    print("=" * 65)
    print("Run-length distributions: BE vs geometric (paper Section 8.3)")
    print()
    for L in RLD_L_VALS:
        r = run_rld_task(L, n_steps=args.n_steps, seed=args.seed)
        print_rld_result(r)


if __name__ == "__main__":
    main()
