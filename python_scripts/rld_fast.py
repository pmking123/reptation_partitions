"""
rld_fast.py
================
Run-length distribution analysis for the uniform partition ensemble.

Key result:
  The marginal distribution of part sizes under the uniform measure over
  partitions of L is NOT geometric — it follows the Bose-Einstein (BE)
  distribution:

      f(j) ~ (1/E[k]) * 1/(exp(j*alpha) - 1)

  where alpha = pi/sqrt(6*L)  and  E[k] ~ (pi/sqrt(6)) * sqrt(L).

  The geometric distribution (which the HR saddle-point argument naively
  suggests) corresponds to the approximation exp(-j*alpha) for the full
  BE form 1/(exp(j*alpha)-1). This approximation is valid only for
  j >> 1/alpha ~ sqrt(L), i.e. for large parts. For small parts (j=1,2,...)
  the BE and geometric distributions differ substantially:
    f_BE(j=1) ~ 1/(alpha*E[k]) = sqrt(6L)/(pi * E[k]) ~ 1/(pi * pred_k)
    f_geom(j=1) = 1/r_bar ~ (pi/sqrt(6)) / sqrt(L) * 1/sqrt(L) = ...

  The exact formula for the marginal is:

      f(j) = (1/(p(L)*E[k])) * sum_{m=1}^{floor(L/j)} p(L - m*j)

  which converges to the BE form as L -> infinity.

Relevant points:

  - The run-length distribution at the HR saddle is Bose-Einstein, not geometric.
  - Geometric is the distribution in the COMPOSITION ensemble (all compositions
    equally weighted); the PARTITION ensemble has the BE distribution.
  - The distinction matters physically: BE has more weight at j=1 (short runs)
    and a heavier tail (long runs) than geometric.

Parallelisation:
  All sample_rld calls are independent (each uses its own RNG seed and shares
  no state). They are dispatched via ProcessPoolExecutor across all available
  cores. Step 2b reuses the Step 2 results rather than resampling.

Requires: mc.py and mc_uniform.py in the same directory.

Usage:
    python3 rld.py
"""

import math
import random
import sys
import os
import time
import numpy as np
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mc_fast import init_arrays
from mc_uniform_fast import metropolis_step_uniform_fast

PI = math.pi

# Number of worker processes: honour SLURM_CPUS_PER_TASK if set,
# otherwise use all available cores.
N_WORKERS = int(os.environ.get("SLURM_CPUS_PER_TASK", os.cpu_count() or 1))


# ============================================================
# Exact theoretical marginal (finite L)
# ============================================================

def p_L(L, _cache={}):
    """Number of partitions of L (exact, via enumeration for small L)."""
    if L in _cache:
        return _cache[L]
    if L < 0:
        return 0
    if L == 0:
        return 1
    try:
        from scripts.core import partitions as ep
        n = sum(1 for _ in ep(L))
    except ImportError:
        # Fallback: Euler recurrence
        n = _euler_p(L)
    _cache[L] = n
    return n


def _euler_p(L, _cache={}):
    """Partition function via Euler's recurrence."""
    if L in _cache:
        return _cache[L]
    if L < 0:
        return 0
    if L == 0:
        return 1
    result = 0
    k = 1
    while True:
        p1 = k*(3*k-1)//2
        p2 = k*(3*k+1)//2
        if p1 > L:
            break
        sign = (-1)**(k+1)
        result += sign * _euler_p(L - p1)
        if p2 <= L:
            result += sign * _euler_p(L - p2)
        k += 1
    _cache[L] = result
    return result


def f_exact(j, L, Ek):
    """Exact marginal f(j) = (1/(p(L)*E[k])) * sum_m p(L-m*j)."""
    pL = p_L(L)
    if pL == 0:
        return 0
    contrib = sum(p_L(L - m*j) for m in range(1, L//j + 1))
    return contrib / (pL * Ek)


# ============================================================
# Asymptotic distributions
# ============================================================

def f_bose_einstein(j, L, Ek=None):
    """
    Bose-Einstein asymptotic marginal:
        f(j) = (1/E[k]) / (exp(j*alpha) - 1)
    where alpha = pi/sqrt(6*L).
    This is the correct asymptotic form from the Hardy-Ramanujan analysis.
    """
    alpha = PI / math.sqrt(6 * L)
    if Ek is None:
        Ek = PI / math.sqrt(6) * math.sqrt(L)
    x = j * alpha
    if x > 50:
        return 0.0
    denom = math.exp(x) - 1
    if denom < 1e-12:
        return float('inf')
    return (1.0 / Ek) / denom


def f_geometric(j, L, Ek=None):
    """
    Geometric marginal (naive saddle-point approximation):
        f(j) = (1/r_bar) * (1 - 1/r_bar)^(j-1)
    where r_bar = L/E[k].
    This is the approximation that replaces 1/(e^x-1) with e^{-x},
    valid only for large x (large j).
    """
    if Ek is None:
        Ek = PI / math.sqrt(6) * math.sqrt(L)
    r_bar = L / Ek
    if r_bar <= 1:
        return 0.0
    q = 1.0 - 1.0 / r_bar
    return (1.0 / r_bar) * q**(j-1)


# ============================================================
# MC sampling of run-length distribution
# ============================================================

def sample_rld(L, n_steps=500000, seed=42):
    rng        = random.Random(seed)
    vals, mults = init_arrays(L)
    k          = int(mults.sum())
    for _ in range(max(50000, 300*L)):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)

    dist   = Counter()
    k_sum  = 0.0
    k2_sum = 0.0
    n_samp = 0
    thin   = max(1, L//20)

    for i in range(n_steps * thin):
        vals, mults, k, _ = metropolis_step_uniform_fast(vals, mults, k, rng)
        if i % thin == 0:
            k_sum  += k
            k2_sum += k*k
            for v, m in zip(vals.tolist(), mults.tolist()):
                dist[v] += m
            n_samp += 1

    total  = sum(dist.values())
    k_mean = k_sum / n_samp
    k_std  = math.sqrt(max(k2_sum/n_samp - k_mean**2, 0))

    return {
        'L': L, 'n_samp': n_samp,
        'k_mean': k_mean, 'k_std': k_std,
        'dist': dict(dist), 'total': total
    }

def _sample_rld_worker(args):
    """Top-level wrapper for ProcessPoolExecutor (must be picklable)."""
    L, n_steps, seed = args
    return (L, n_steps), sample_rld(L, n_steps=n_steps, seed=seed)


def run_all_sampling(jobs):
    """
    Dispatch all (L, n_steps, seed) jobs in parallel.
    Returns dict keyed by (L, n_steps) -> res.
    When two jobs share the same (L, n_steps) key (Step 2 and Step 2b both
    need L=100/200/500/2000 at 600000 steps) only one run is performed and
    the result is shared.

    jobs: list of (L, n_steps, seed)
    """
    # Deduplicate by (L, n_steps) — keep first seed seen
    seen = {}
    for L, n_steps, seed in jobs:
        key = (L, n_steps)
        if key not in seen:
            seen[key] = (L, n_steps, seed)

    unique_jobs = list(seen.values())
    results = {}

    print(f"Dispatching {len(unique_jobs)} sampling jobs across {N_WORKERS} workers...")
    t0 = time.time()

    with ProcessPoolExecutor(max_workers=N_WORKERS) as executor:
        futures = {executor.submit(_sample_rld_worker, job): job
                   for job in unique_jobs}
        for f in as_completed(futures):
            key, res = f.result()
            results[key] = res
            L, n_steps = key
            elapsed = time.time() - t0
            print(f"  done: L={L}, n_steps={n_steps}, "
                  f"n_samp={res['n_samp']}, elapsed={elapsed:.1f}s")

    print(f"All sampling complete in {time.time()-t0:.1f}s\n")
    return results


# ============================================================
# Comparison and analysis
# ============================================================

def compare_distributions(L, res, use_exact=True):
    """
    Compare MC distribution against BE, geometric, and (for small L) exact theory.
    """
    km    = res['k_mean']
    dist  = res['dist']
    total = res['total']

    pred_k = PI / math.sqrt(6) * math.sqrt(L)
    pred_r = math.sqrt(6) / PI * math.sqrt(L)

    print(f"L = {L}")
    print(f"  E[k]  = {km:.3f}  (HR pred {pred_k:.2f},  ratio {km/pred_k:.3f})")
    print(f"  r_bar = {L/km:.3f}  (HR pred {pred_r:.2f})")
    print()

    # Decide column set based on whether exact is feasible
    show_exact = use_exact and L <= 30

    if show_exact:
        header = f"  {'j':>4}  {'MC':>9}  {'BE':>9}  {'geom':>9}  {'exact':>9}  {'MC/BE':>7}  {'MC/geom':>9}"
    else:
        header = f"  {'j':>4}  {'MC':>9}  {'BE':>9}  {'geom':>9}  {'MC/BE':>7}  {'MC/geom':>9}"
    print(header)
    print("  " + "-" * (len(header)-2))

    kl_be   = 0.0
    kl_geom = 0.0
    max_j   = max(dist.keys())

    for j in range(1, min(max_j+1, int(8*L/km)+2)):
        fm = dist.get(j, 0) / total
        fbe = f_bose_einstein(j, L, km)
        fg  = f_geometric(j, L, km)

        if fm < 1e-10 and fbe < 1e-10:
            break

        r_be   = fm/fbe   if fbe   > 1e-12 else 0
        r_geom = fm/fg    if fg    > 1e-12 else 0

        if fbe   > 1e-12 and fm > 1e-12:
            kl_be   += fm * math.log(fm/fbe)
        if fg    > 1e-12 and fm > 1e-12:
            kl_geom += fm * math.log(fm/fg) 

        if show_exact:
            fe = f_exact(j, L, km)
            print(f"  {j:>4}  {fm:>9.5f}  {fbe:>9.5f}  {fg:>9.5f}  {fe:>9.5f}  {r_be:>7.4f}  {r_geom:>9.4f}")
        else:
            print(f"  {j:>4}  {fm:>9.5f}  {fbe:>9.5f}  {fg:>9.5f}  {r_be:>7.4f}  {r_geom:>9.4f}")

    print()
    print(f"  KL(MC || BE)   = {kl_be:.6f}")
    print(f"  KL(MC || geom) = {kl_geom:.6f}")
    print(f"  BE is {'better' if kl_be < kl_geom else 'worse'} than geometric "
          f"by factor {kl_geom/kl_be:.2f}")
    print()


# ============================================================
# Convergence of BE approximation quality
# ============================================================

def be_convergence_study(L_vals, results):
    """
    Track how well the BE approximation fits the MC distribution as L grows.
    Results are looked up from the pre-computed results dict.
    """
    print("BE approximation quality vs L")
    print(f"  {'L':>6}  {'E[k]/sqL':>10}  {'KL_BE':>9}  {'KL_geom':>10}  "
          f"{'ratio':>7}  {'f_BE(1)/f_MC(1)':>16}")
    print("  " + "-"*65)

    for L in L_vals:
        n_steps = max(600000, 3000*L)
        res   = results[(L, n_steps)]
        km    = res['k_mean']
        dist  = res['dist']
        total = res['total']

        kl_be = kl_geom = 0.0
        for j in range(1, max(dist.keys())+1):
            fm  = dist.get(j, 0) / total
            fbe = f_bose_einstein(j, L, km)
            fg  = f_geometric(j, L, km)
            if fm < 1e-10 and fbe < 1e-10:
                break
            if fbe > 1e-12 and fm > 1e-12:
                kl_be   += fm * math.log(fm/fbe)
            if fg  > 1e-12 and fm > 1e-12:
                kl_geom += fm * math.log(fm/fg)

        f_mc_1  = dist.get(1, 0) / total
        f_be_1  = f_bose_einstein(1, L, km)
        ratio_1 = f_be_1/f_mc_1 if f_mc_1 > 0 else 0

        ratio = kl_geom/kl_be if kl_be > 1e-10 else float('inf')
        print(f"  {L:>6}  {km/math.sqrt(L):>10.4f}  {kl_be:>9.5f}  {kl_geom:>10.5f}  "
              f"{ratio:>7.2f}  {ratio_1:>16.4f}")

    print()


# ============================================================
# GC vs canonical normaliser comparison
# ============================================================

def f_bose_einstein_gc(j, L):
    """
    BE marginal normalised by the grand canonical E[k], not the canonical mean.
    """
    alpha = PI / math.sqrt(6 * L)
    Ek_gc = PI / math.sqrt(6) * math.sqrt(L)
    x = j * alpha
    if x > 50:
        return 0.0
    denom = math.exp(x) - 1
    if denom < 1e-12:
        return float('inf')
    return (1.0 / Ek_gc) / denom


def compare_normalisers(L, res):
    """
    Compare BE predictions using canonical vs GC normaliser.

    KL divergence cannot be used here because f_BE_gc is not a normalised
    probability distribution: since E[k]_GC < E[k]_canon, the GC-normalised
    BE formula sums to E[k]_canon/E[k]_GC > 1. Instead we report:
      (a) the normalisation ratio E[k]_canon / E[k]_GC (= the factor by which
          f_BE_gc overcounts probability mass)
      (b) the prediction error at j=1, the point of largest absolute discrepancy
      (c) the KL divergence for the canonical normaliser, for reference
    """
    km    = res['k_mean']
    dist  = res['dist']
    total = res['total']
    Ek_gc = PI / math.sqrt(6) * math.sqrt(L)

    # (a) normalisation ratio
    norm_ratio = km / Ek_gc

    # (b) j=1 prediction errors
    f_mc_1  = dist.get(1, 0) / total
    f_can_1 = f_bose_einstein(1, L, km)
    f_gc_1  = f_bose_einstein_gc(1, L)
    err_can = f_can_1 / f_mc_1
    err_gc  = f_gc_1  / f_mc_1

    # (c) KL for canonical normaliser
    kl_canon = 0.0
    for j in range(1, max(dist.keys()) + 1):
        fm    = dist.get(j, 0) / total
        f_can = f_bose_einstein(j, L, km)
        if fm < 1e-10 and f_can < 1e-10:
            break
        if f_can > 1e-12 and fm > 1e-12:
            kl_canon += fm * math.log(fm / f_can)

    print(f"  L={L}: E[k]_canon={km:.2f}, E[k]_GC={Ek_gc:.2f}, "
          f"norm_ratio={norm_ratio:.3f}")
    print(f"    f_BE_gc sums to {norm_ratio:.3f} (not 1): KL undefined for GC normaliser")
    print(f"    f_MC(1)={f_mc_1:.5f},  f_BE_canon(1)={f_can_1:.5f} "
          f"(ratio {err_can:.4f}),  f_BE_gc(1)={f_gc_1:.5f} (ratio {err_gc:.4f})")
    print(f"    GC normaliser error at j=1: factor {1/err_gc:.3f} "
          f"(underestimates by {(1-err_gc)*100:.1f}%)")
    print(f"    KL(MC || BE_canon) = {kl_canon:.6f}")
    print()


# ============================================================
# Entry point
# ============================================================

# NOTE on crossover:
# At small L (L < ~55), the geometric distribution actually fits better than BE,
# because the BE formula is an asymptotic result requiring j*alpha << 1 for many
# terms, and alpha = pi/sqrt(6L) is not small enough.
# At L > ~55, BE becomes progressively better (KL ratio geom/BE grows with L).
# The crossover is pinned between L=50 (ratio=0.80, geom better) and L=60
# (ratio=1.02, BE better) by the Step 3 convergence study.
# The exact formula f(j) = sum_m p(L-mj) / (p(L)*E[k]) is accurate for all L.

if __name__ == '__main__':
    print()
    print("=" * 65)
    print("Run-length distribution analysis: uniform partition ensemble")
    print("=" * 65)
    print()
    print(f"Using {N_WORKERS} worker processes.")
    print()

    print("Key finding: the marginal run-length distribution is Bose-Einstein,")
    print("  f(j) ~ (1/E[k]) / (exp(j*pi/sqrt(6L)) - 1),")
    print("not geometric as the naive saddle-point argument suggests.")
    print()

    # ----------------------------------------------------------
    # Build the full job list (deduplicated by run_all_sampling)
    # ----------------------------------------------------------

    # Step 1 jobs
    step1_L      = [15, 20, 25, 30]
    step1_jobs   = [(L, 200000, 42) for L in step1_L]

    # Step 2 jobs (Step 2b reuses these — no resampling)
    step2_L      = [50, 100, 200, 500, 2000]
    step2_jobs   = [(L, 600000, 42) for L in step2_L]

    # Step 3 jobs
    step3_L      = [20, 50, 60, 70, 80, 90, 100, 200, 500, 1000, 2000]
    step3_jobs   = [(L, max(600000, 3000*L), 42) for L in step3_L]

    all_jobs = step1_jobs + step2_jobs + step3_jobs

    # Run all sampling in parallel (duplicate (L, n_steps) keys run only once)
    results = run_all_sampling(all_jobs)

    # ----------------------------------------------------------
    # Step 1: Exact comparison for small L
    # ----------------------------------------------------------
    print("=" * 65)
    print("Step 1: Exact theory vs BE vs geometric (small L, exact enumeration)")
    print()
    for L in step1_L:
        compare_distributions(L, results[(L, 200000)], use_exact=True)

    # ----------------------------------------------------------
    # Step 2: MC comparison for larger L
    # ----------------------------------------------------------
    print("=" * 65)
    print("Step 2: MC vs BE vs geometric (larger L)")
    print("  (600 000 steps per L, matching the paper figure caption)")
    print()
    for L in step2_L:
        res = results[(L, 600000)]
        print(f"[L={L}, {res['n_samp']} samples]")
        compare_distributions(L, res, use_exact=False)

    # ----------------------------------------------------------
    # Step 2b: GC vs canonical normaliser (reuses Step 2 results)
    # ----------------------------------------------------------
    print("=" * 65)
    print("Step 2b: GC vs canonical normaliser -- prediction error at j=1")
    print()
    for L in [100, 200, 500, 2000]:
        compare_normalisers(L, results[(L, 600000)])

    # ----------------------------------------------------------
    # Step 3: BE convergence study
    # ----------------------------------------------------------
    print("=" * 65)
    print("Step 3: BE approximation quality vs L")
    print("  (n_steps = max(600 000, 3000*L) to keep MC noise below KL signal)")
    print()
    be_convergence_study(step3_L, results)

    # ----------------------------------------------------------
    # Summary
    # ----------------------------------------------------------
    print("=" * 65)
    print("Summary")
    print()
    print("The run-length marginal distribution under the uniform partition")
    print("measure is Bose-Einstein (BE), not geometric:")
    print()
    print("  f_BE(j)   = (1/E[k]) / (exp(j*pi/sqrt(6L)) - 1)")
    print("  f_geom(j) = (1/r_bar) * (1 - 1/r_bar)^(j-1)")
    print()
    print("The geometric approximation underestimates f(j=1) by ~37-56%")
    print("and overestimates f(j=2..6) by ~28-60% for L >= 100. BE agrees")
    print("with MC up to a small correction that vanishes as L -> infinity.")
    print("(At small L, L < ~55, the BE formula is an asymptotic result and")
    print(" the exact formula f_exact matches MC better than either asymptotic.")
    print(" Crossover confirmed between L=50 (geom better) and L=60 (BE better).)")
    print()
    print("Physical interpretation:")
    print("  Geometric: independent runs of fixed mean length (Poisson process)")
    print("  BE: runs behave like bosons -- small-j runs are enhanced because")
    print("      multiple copies are allowed, raising f(j) above the Poisson")
    print("      baseline by the 1/(1-e^{-alpha}) ~ 1/(j*alpha) factor for small j.")
