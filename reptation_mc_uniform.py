"""
reptation_mc_uniform.py
=======================
Metropolis-Hastings sampler over integer partitions of L under
the UNIFORM measure — every partition of L has equal weight 1/p(L).

Purpose
-------
Test the Hardy-Ramanujan prediction for the most probable partition:
    modal k  ~  (pi/sqrt(6)) * sqrt(L)  ~  1.2825 * sqrt(L)
    mean  k  converges much more slowly (not to 1.28) — see notes below.

Key finding from exact enumeration:
    L=10: exact modal k = 4,  HR predicts 4.06  (error < 1.5%)
    L=15: exact modal k = 5,  HR predicts 4.97  (error < 1%)
    L=20: exact modal k = 6,  HR predicts 5.74  (error < 5%)

The mean k/sqrt(L) converges to ~1.65 at L=20 and grows further —
the HR saddle point governs the MODE, not the MEAN.

Algorithm
---------
Metropolis-Hastings with three move types, all O(n_distinct) per step:

  Transfer:   a -> a+1  and  b -> b-1       (k, L preserved)
  Split:      part a -> (c, a-c)            (k+1, L preserved)
  Merge:      parts c,d -> c+d             (k-1, L preserved)

Acceptance: min(1, N_fwd / N_rev) where N_fwd and N_rev are computed
using O(1) count formulas (no move enumeration needed).

Requires: reptation_mc.py in the same directory.

Usage
-----
    python3 reptation_mc_uniform.py
"""

import math
import random
import time
from collections import Counter
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from reptation_mc import (
    Partition, init_partition,
    apply_transfer, apply_split, apply_merge,
    all_transfers, all_splits, all_merges,
)


# ============================================================
# O(n_distinct) move count formulas
# ============================================================

def n_transfers_count(p):
    """Number of valid transfer moves in O(n_distinct)."""
    vals   = set(p.counts.keys())
    n_v    = len(vals)
    dec    = [v for v in vals if v >= 2]
    n_d    = len(dec)
    n_adj  = sum(1 for v in vals if v + 1 in vals)
    n_ei   = sum(1 for v in dec if p.counts[v] == 1)
    return n_v * n_d - n_adj - n_ei


def n_splits_count(p):
    """Number of valid split moves: sum floor(a/2) over distinct parts."""
    return sum(a // 2 for a in p.counts)


def n_merges_count(p):
    """Number of valid merge moves (distinct unordered pairs)."""
    vals  = list(p.counts.keys())
    n_d   = len(vals)
    n_rep = sum(1 for m in p.counts.values() if m >= 2)
    return n_d * (n_d - 1) // 2 + n_rep


# ============================================================
# Fast uniform-measure Metropolis step
# ============================================================

def metropolis_step_uniform(p, rng, p_transfer=0.5):
    """
    One MH step under the uniform measure over partitions of L.
    O(n_distinct) per step using direct sampling and count formulas.
    Acceptance: min(1, N_fwd / N_rev).
    """
    u    = rng.random()
    vals = list(p.counts.keys())
    n_d  = len(vals)

    if u < p_transfer:
        # Transfer: sample (a,b) with a->a+1, b->b-1
        dec = [v for v in vals if v >= 2]
        if not dec:
            return p, False
        for _ in range(30):
            a = rng.choice(vals)
            b = rng.choice(dec)
            if a + 1 == b:
                continue   # null move
            if a == b and p.counts[a] < 2:
                continue   # need two copies for a==b
            break
        else:
            return p, False
        new_p = apply_transfer(p, a, b)
        N_fwd = n_transfers_count(p)
        N_rev = n_transfers_count(new_p)

    elif u < p_transfer + (1 - p_transfer) / 2:
        # Split: sample (a,c) proportional to floor(a/2)
        splittable = [(v, v // 2) for v in vals if v >= 2]
        if not splittable:
            return p, False
        total = sum(w for _, w in splittable)
        r = rng.randint(1, total)
        cumul = 0
        for v, w in splittable:
            cumul += w
            if r <= cumul:
                a = v
                c = rng.randint(1, v // 2)
                break
        new_p = apply_split(p, a, c)
        N_fwd = n_splits_count(p)
        N_rev = n_merges_count(new_p)

    else:
        # Merge: sample unordered pair {c,d}
        n_same = sum(1 for m in p.counts.values() if m >= 2)
        n_diff = n_d * (n_d - 1) // 2
        n_merge = n_same + n_diff
        if n_merge == 0:
            return p, False
        r = rng.randint(1, n_merge)
        if r <= n_diff:
            i = rng.randint(0, n_d - 1)
            j = rng.randint(0, n_d - 2)
            if j >= i:
                j += 1
            c, d = vals[i], vals[j]
        else:
            rep = [v for v, m in p.counts.items() if m >= 2]
            v   = rng.choice(rep)
            c = d = v
        new_p = apply_merge(p, min(c, d), max(c, d))
        N_fwd = n_merges_count(p)
        N_rev = n_splits_count(new_p)

    if N_rev == 0:
        return p, False
    log_alpha = math.log(N_fwd) - math.log(N_rev)
    if log_alpha >= 0 or rng.random() < math.exp(log_alpha):
        return new_p, True
    return p, False


# ============================================================
# Autocorrelation
# ============================================================

def autocorrelation(series, max_lag=500):
    n    = len(series)
    ml   = min(max_lag, n // 4)
    mean = sum(series) / n
    var  = sum((x - mean)**2 for x in series) / n
    if var < 1e-12:
        return [1.0] + [0.0] * ml
    return [
        sum((series[i] - mean) * (series[i + t] - mean)
            for i in range(n - t)) / (var * (n - t))
        for t in range(ml + 1)
    ]


def tau_int(acf):
    tau = 0.5
    for t in range(1, len(acf)):
        tau += acf[t]
        if t >= 5 * tau:
            break
    return tau


# ============================================================
# Validation against exact uniform distribution
# ============================================================

def validate(L_max=20, n_mc=300000):
    try:
        from reptation_core import partitions as ep
    except ImportError:
        print("reptation_core.py not found; skipping.")
        return

    print("Validation: MC vs exact uniform distribution")
    print()

    PI     = math.pi
    pred_k = PI / math.sqrt(6)

    print(f"  {'L':>4}  {'p(L)':>6}  {'modal_k_ex':>12}  "
          f"{'modal_k_HR':>12}  {'<k>_ex':>8}  {'<k>_MC':>8}  "
          f"{'err%':>6}  {'TV':>7}")
    print(f"  {'-'*78}")

    for L in range(4, L_max + 1, 2):
        all_p  = list(ep(L))
        pL     = len(all_p)
        k_cnt  = Counter(len(p) for p in all_p)
        modal_ex = k_cnt.most_common(1)[0][0]
        km_ex  = sum(len(p) for p in all_p) / pL
        sqL    = math.sqrt(L)

        rng  = random.Random(7)
        part = init_partition(L, rng)
        for _ in range(50000):
            part, _ = metropolis_step_uniform(part, rng)

        visits  = Counter()
        k_samp  = []
        for _ in range(n_mc):
            part, _ = metropolis_step_uniform(part, rng)
            k_samp.append(part.k)
            visits[part.to_tuple()] += 1

        km_mc   = sum(k_samp) / n_mc
        tv      = sum(abs(1/pL - visits.get(p, 0)/n_mc) for p in all_p) / 2
        err     = abs(km_mc - km_ex) / km_ex * 100

        print(f"  {L:>4}  {pL:>6}  {modal_ex:>12}  "
              f"{pred_k*sqL:>12.2f}  {km_ex:>8.4f}  {km_mc:>8.4f}  "
              f"{err:>6.2f}%  {tv:>7.4f}")

    print()
    print("Key: modal_k_ex is the exact most-probable k under uniform measure.")
    print("     modal_k_HR is the Hardy-Ramanujan asymptotic prediction.")
    print("     These agree well even at small L.")
    print()


# ============================================================
# Mode convergence study (main test of HR prediction)
# ============================================================

def mode_convergence_study(L_vals, n_per_L=None, seed=42):
    """
    Track the empirical mode of k versus HR prediction.
    This is the central test: modal_k / sqrt(L) -> pi/sqrt(6)?
    """
    PI     = math.pi
    pred_k = PI / math.sqrt(6)

    print("Mode convergence study")
    print(f"  HR prediction: modal k / sqrt(L) -> {pred_k:.4f}")
    print()

    hdr = (f"{'L':>7}  {'modal_k':>8}  {'modal/sqL':>10}  "
           f"{'err%':>7}  {'<k>/sqL':>9}  {'sig/L^.25':>10}  "
           f"{'tau':>7}  {'t(s)':>5}")
    print(hdr)
    print("-" * len(hdr))

    results = []
    for L in L_vals:
        n      = n_per_L or max(1000000, 5000 * L)
        thin   = max(1, L // 20)
        burnin = max(50000, 300 * L)

        t0   = time.time()
        rng  = random.Random(seed)
        part = init_partition(L, rng)
        for _ in range(burnin):
            part, _ = metropolis_step_uniform(part, rng)

        k_list  = []
        k_cnt   = Counter()
        for i in range(n * thin):
            part, _ = metropolis_step_uniform(part, rng)
            if i % thin == 0:
                k_list.append(part.k)
                k_cnt[part.k] += 1
        t1 = time.time()

        sqL    = math.sqrt(L)
        L14    = L ** 0.25
        km     = sum(k_list) / len(k_list)
        kvar   = sum((k - km)**2 for k in k_list) / len(k_list)
        ksig   = math.sqrt(max(kvar, 0))
        acf_   = autocorrelation(k_list)
        ti     = tau_int(acf_)
        mk     = k_cnt.most_common(1)[0][0]
        me     = (mk / sqL - pred_k) / pred_k * 100

        print(f"{L:>7}  {mk:>8}  {mk/sqL:>10.4f}  "
              f"{me:>+7.2f}%  {km/sqL:>9.4f}  {ksig/L14:>10.4f}  "
              f"{ti:>7.1f}  {t1-t0:>5.1f}")

        results.append({
            'L': L, 'modal_k': mk, 'km': km, 'ksig': ksig,
            'tau': ti, 'k_cnt': k_cnt
        })

    return results


# ============================================================
# Run-length distribution
# ============================================================

def run_length_distribution(L, n_steps=500000, seed=42):
    PI     = math.pi
    pred_k = PI / math.sqrt(6) * math.sqrt(L)
    pred_r = math.sqrt(6) / PI * math.sqrt(L)

    rng  = random.Random(seed)
    part = init_partition(L, rng)
    for _ in range(n_steps // 5):
        part, _ = metropolis_step_uniform(part, rng)

    dist  = Counter()
    k_sum = 0
    for _ in range(n_steps):
        part, _ = metropolis_step_uniform(part, rng)
        k_sum  += part.k
        for v, m in part.counts.items():
            dist[v] += m

    total  = sum(dist.values())
    k_mean = k_sum / n_steps
    r_bar  = L / k_mean
    q      = max(0, 1 - 1 / r_bar)

    print(f"  L={L}: <k>={k_mean:.2f} (HR pred {pred_k:.2f}), "
          f"r_bar={r_bar:.2f} (HR pred {pred_r:.2f})")
    print(f"  {'j':>4}  {'f_meas':>10}  {'geom_pred':>10}  {'ratio':>7}")
    for j in range(1, min(max(dist) + 1, int(6 * r_bar) + 2)):
        fm = dist.get(j, 0) / total
        fp = (1 / r_bar) * q ** (j - 1) if q > 0 else 0
        if fm < 5e-5 and fp < 5e-5:
            break
        print(f"  {j:>4}  {fm:>10.5f}  {fp:>10.5f}  "
              f"{fm/fp if fp > 1e-9 else 0:>7.4f}")
    print()


# ============================================================
# Speed benchmark
# ============================================================

def benchmark(L_vals=None):
    if L_vals is None:
        L_vals = [100, 500, 1000, 5000, 10000]

    print("Speed benchmark (steps/sec):")
    for L in L_vals:
        p   = init_partition(L, random.Random(1))
        rng = random.Random(42)
        for _ in range(2000):
            p, _ = metropolis_step_uniform(p, rng)
        n  = 5000
        t0 = time.time()
        for _ in range(n):
            p, _ = metropolis_step_uniform(p, rng)
        t1 = time.time()
        print(f"  L={L:>6}: {n/(t1-t0):>9.0f} steps/sec  k~{p.k}")
    print()


# ============================================================
# Entry point
# ============================================================


# ============================================================
# Exact modal k study (no MC needed for small L)
# ============================================================

def exact_modal_study(L_max=50):
    try:
        from reptation_core import partitions as ep
    except ImportError:
        print("reptation_core.py not found.")
        return

    PI     = math.pi
    pred_k = PI / math.sqrt(6)

    print("Exact modal k vs Hardy-Ramanujan prediction")
    print(f"  HR asymptotic: modal_k / sqrt(L) -> {pred_k:.4f}")
    print()
    print(f"  {'L':>5}  {'p(L)':>7}  {'modal_k':>8}  {'modal/sqL':>11}  {'HR_pred':>9}  {'err%':>7}  {'freq':>7}")
    print("  " + "-"*60)

    for L in range(4, L_max + 1):
        try:
            parts = list(ep(L))
        except Exception:
            break
        k_cnt   = Counter(len(p) for p in parts)
        pL      = len(parts)
        mk      = k_cnt.most_common(1)[0][0]
        mk_freq = k_cnt[mk] / pL
        sqL     = math.sqrt(L)
        hr      = pred_k * sqL
        err     = (mk / sqL - pred_k) / pred_k * 100

        print(f"  {L:>5}  {pL:>7}  {mk:>8}  "
              f"{mk/sqL:>11.4f}  {hr:>9.2f}  {err:>+7.2f}%  {mk_freq:>7.4f}")

    print()


if __name__ == '__main__':
    print()
    print("=" * 65)
    print("Uniform partition sampler: Hardy-Ramanujan verification")
    print("=" * 65)
    print()

    # Step 1: Validate
    validate(L_max=20, n_mc=200000)

    # Step 2: Speed benchmark
    benchmark([100, 500, 2000, 5000, 10000])

    # Step 3: Mode convergence (main result)
    print("=" * 65)
    mode_convergence_study(
        [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000],
        seed=42
    )

    # Step 4: Run-length distributions
    print()
    print("=" * 65)
    print("Run-length distributions (checking geometric shape)")
    print()
    for L in [100, 500, 2000]:
        run_length_distribution(L, n_steps=500000, seed=42)
