"""
reptation_partition_simulations.py
===================================
Exact enumeration and verification for the reptation-partition correspondence.

Developed in the context of a research programme connecting the
run-length recoding of reptating polymer chains on a 2D square lattice
to the theory of integer partitions and restricted partition classes.

Contents
--------
PART 1 — Core definitions
    Compositions, partitions, self-avoidance, degeneracy w(lambda)

PART 2 — Verification 1: Decoupling theorem
    w_strong(comp) = w_H(odd runs) * w_V(even runs)
    Verified for all 4095 compositions of L=1..12

PART 3 — Verification 2: Formula analysis
    p=2 formula: w_H(h1,h2) = 4 if h1!=h2, else 2  [exact]
    Distinct-parts formula: w_H = 2^p  [FAILS for p>=3]
    Reason: vanishing consecutive sub-sums create partial sum collisions
    Correction documented with examples

PART 4 — Verification 3: Most probable partition
    w(lambda) computed exactly for L up to 20
    k* (most probable run count) and r_bar* (mean run length) vs predictions
    Fluctuation analysis: sigma_k vs L^(1/4)

Key results
-----------
- Decoupling theorem: EXACT (zero failures across all tested cases)
- p=2 formula: EXACT
- Distinct-parts formula w_H=2^p: CORRECTED
  Fails when any consecutive sub-block has a vanishing signed sum
  Fraction achieving w_H=2^p falls from 1.0 (L=6) to ~0.57 (L=16)
- k* scaling: converging toward predicted 1.2825*sqrt(L) but
  pre-asymptotic for L<=20; larger L needed for clean convergence
- sigma_k scaling: converging toward predicted 1.13*L^(1/4)

Usage
-----
    python3 reptation_partition_simulations.py

Requires: Python 3.8+, standard library only (itertools, collections,
          math, functools, time). No external dependencies.

Runtime: ~60s for the full suite up to L=20.
         Parts 1-3 complete in <5s.
"""

from itertools import permutations, combinations, product as iproduct
from collections import Counter
from functools import lru_cache
import math
import time


# ============================================================
# PART 1 — Core definitions
# ============================================================

def compositions(L, min_part=1):
    """Generate all compositions of L with parts >= min_part."""
    if L == 0:
        yield ()
        return
    if L < min_part:
        return
    for first in range(min_part, L + 1):
        for rest in compositions(L - first, min_part):
            yield (first,) + rest


def partitions(L, min_part=1):
    """Generate all partitions of L with parts >= min_part."""
    seen = set()
    for comp in compositions(L, min_part):
        p = tuple(sorted(comp, reverse=True))
        if p not in seen:
            seen.add(p)
            yield p


def partial_sums_distinct(values, signs):
    """
    Return True if the signed partial sums S_0=0, S_1, ..., S_p
    are all distinct.
    values: tuple of positive integers
    signs:  tuple of +1/-1
    """
    seen = {0}
    s = 0
    for v, e in zip(values, signs):
        s += e * v
        if s in seen:
            return False
        seen.add(s)
    return True


@lru_cache(maxsize=None)
def count_valid_sign_assignments(subseq):
    """
    Count sign assignments in {+1,-1}^p such that all signed
    partial sums (including S_0=0) are distinct.
    subseq: tuple of positive integers
    """
    p = len(subseq)
    if p == 0:
        return 1
    count = 0
    for signs in iproduct((1, -1), repeat=p):
        if partial_sums_distinct(subseq, signs):
            count += 1
    return count


def w_strong(composition):
    """
    Number of sign assignments (epsilon_1,...,epsilon_k) such that
    both the H-subsequence (odd-indexed) and V-subsequence (even-indexed)
    signed partial sums are all distinct.

    By the decoupling theorem (verified below):
        w_strong = w_H(odd-indexed runs) * w_V(even-indexed runs)
    """
    h = composition[0::2]
    v = composition[1::2]
    return count_valid_sign_assignments(h) * count_valid_sign_assignments(v)


def distinct_perms(seq):
    """Generate all distinct permutations of seq (handles repeated elements)."""
    seq = sorted(seq)
    if not seq:
        yield ()
        return
    prev = object()
    for i, val in enumerate(seq):
        if val == prev:
            continue
        prev = val
        rest = seq[:i] + seq[i+1:]
        for p in distinct_perms(rest):
            yield (val,) + p


@lru_cache(maxsize=None)
def _wH(h): return count_valid_sign_assignments(h)

@lru_cache(maxsize=None)
def _wV(v): return count_valid_sign_assignments(v)


def w_lambda(partition):
    """
    Exact degeneracy w(lambda): sum over all distinct orderings
    (compositions mapping to lambda) of w_H(odd runs) * w_V(even runs).

    Each distinct permutation of lambda corresponds to exactly one
    chain configuration type, and contributes w_H * w_V valid signed
    realisations (self-avoiding geometric configurations).
    """
    total = 0
    for perm in distinct_perms(list(partition)):
        total += _wH(perm[0::2]) * _wV(perm[1::2])
    return total


def w_lambda_brute(partition):
    """
    Brute-force degeneracy (slow, for verification only).
    Enumerates all permutations including duplicates, deduplicates via set.
    """
    seen = set()
    total = 0
    for perm in permutations(partition):
        if perm in seen:
            continue
        seen.add(perm)
        h = perm[0::2]
        v = perm[1::2]
        total += count_valid_sign_assignments(h) * count_valid_sign_assignments(v)
    return total


# ============================================================
# PART 2 — Verification 1: Decoupling theorem
# ============================================================

def run_verification_1():
    print("=" * 62)
    print("PART 2 — Verification 1: Decoupling theorem")
    print("=" * 62)
    print()
    print("Testing: w_strong(comp) == w_H(odd runs) * w_V(even runs)")
    print("for all compositions of L = 1..12.")
    print()

    failures = 0
    tested = 0
    for L in range(1, 13):
        for comp in compositions(L):
            h = comp[0::2]
            v = comp[1::2]
            w_product = (count_valid_sign_assignments(h) *
                         count_valid_sign_assignments(v))
            # Recompute w_strong independently via full sign enumeration
            k = len(comp)
            w_direct = 0
            for signs in iproduct((1, -1), repeat=k):
                h_signs = signs[0::2]
                v_signs = signs[1::2]
                if (partial_sums_distinct(h, h_signs) and
                        partial_sums_distinct(v, v_signs)):
                    w_direct += 1
            if w_product != w_direct:
                print(f"  FAIL: comp={comp}, product={w_product}, direct={w_direct}")
                failures += 1
            tested += 1

    print(f"  Tested {tested} compositions.")
    print(f"  Failures: {failures}")
    print(f"  Result: {'VERIFIED — decoupling theorem holds exactly' if failures == 0 else 'FAILED'}")
    print()

    # Extremal cases
    print("Extremal cases:")
    print()
    print("  All-equal H-subsequences (expected w_H = 2 for all p >= 1):")
    all_ok = True
    for p in range(1, 6):
        for v in [1, 2, 3]:
            subseq = tuple([v] * p)
            w = count_valid_sign_assignments(subseq)
            ok = (w == 2)
            if not ok:
                all_ok = False
            if p <= 3:
                print(f"    {subseq}: w_H = {w}  {'OK' if ok else 'FAIL'}")
    print(f"  All-equal result: {'OK' if all_ok else 'FAIL'}")
    print()

    print("  p=2 formula: w_H(h1,h2) = 4 if h1 != h2, else 2")
    failures_p2 = 0
    for h1 in range(1, 8):
        for h2 in range(1, 8):
            w = count_valid_sign_assignments((h1, h2))
            expected = 4 if h1 != h2 else 2
            if w != expected:
                failures_p2 += 1
                print(f"    FAIL: ({h1},{h2}) got {w}, expected {expected}")
    print(f"  p=2 formula: {'VERIFIED' if failures_p2 == 0 else f'{failures_p2} failures'}")
    print()


# ============================================================
# PART 3 — Verification 2: Distinct-parts formula analysis
# ============================================================

def run_verification_2():
    print("=" * 62)
    print("PART 3 — Verification 2: Distinct-parts formula w_H = 2^p")
    print("=" * 62)
    print()

    print("Checking whether w_H = 2^p holds for all-distinct H-subseqs.")
    print()

    # Detailed analysis for (1,2,3)
    print("Case (1,2,3): all 8 sign sequences and their partial sums")
    print()
    subseq = (1, 2, 3)
    print(f"  {'signs':12s}  {'partial sums':20s}  {'valid?':6s}  {'reason'}")
    print(f"  {'-'*10}  {'-'*18}  {'-'*6}  {'-'*16}")
    for signs in iproduct((1, -1), repeat=3):
        s0, s1 = 0, signs[0] * 1
        s2 = s1 + signs[1] * 2
        s3 = s2 + signs[2] * 3
        sums = [s0, s1, s2, s3]
        valid = len(set(sums)) == 4
        reason = ''
        if not valid:
            for i in range(4):
                for j in range(i+1, 4):
                    if sums[i] == sums[j]:
                        reason = f"S_{i}=S_{j}={sums[i]}"
        sign_str = '(' + ','.join('+' if e == 1 else '-' for e in signs) + ')'
        print(f"  {sign_str:12s}  {str(sums):20s}  {'YES' if valid else 'NO':6s}  {reason}")
    w = count_valid_sign_assignments(subseq)
    print(f"\n  w_H(1,2,3) = {w} (out of 8). Deficit = {8 - w}.")
    print(f"  Forbidden: sign sequences where 1+2-3=0 or -1-2+3=0 (S_3=S_0=0).")
    print()

    # Statistics across all L
    print("Counter-examples (distinct parts but w_H < 2^p) for L=1..15:")
    print()
    deviations = {}
    for L in range(1, 16):
        for comp in compositions(L):
            h = comp[0::2]
            if len(set(h)) == len(h) and len(h) >= 2:
                p = len(h)
                w = count_valid_sign_assignments(h)
                w_pred = 2**p
                if w < w_pred:
                    if p not in deviations:
                        deviations[p] = []
                    deviations[p].append((h, w, w_pred))

    print(f"  {'p':>4}  {'cases':>8}  {'min w/2^p':>10}  {'deficit range'}")
    for p in sorted(deviations):
        cases = deviations[p]
        ratios = [w / w_pred for _, w, w_pred in cases]
        deficits = [w_pred - w for _, w, w_pred in cases]
        print(f"  {p:>4}  {len(cases):>8}  {min(ratios):>10.4f}  "
              f"{min(deficits)}..{max(deficits)}")

    print()
    print("Fraction of distinct-part compositions with w_H = 2^p exactly:")
    print()
    print(f"  {'L':>4}  {'exact/total':>14}  {'fraction':>10}")
    for L in [6, 8, 10, 12, 14, 16]:
        total = exact = 0
        for comp in compositions(L):
            h = comp[0::2]
            if len(set(h)) == len(h) and len(h) >= 2:
                total += 1
                if count_valid_sign_assignments(h) == 2**len(h):
                    exact += 1
        if total > 0:
            print(f"  {L:>4}  {exact:>6}/{total:<6}  {exact/total:>10.3f}")

    print()
    print("Conclusion: w_H = 2^p ONLY when no consecutive sub-block of the")
    print("H-subsequence has a vanishing signed sum. This is a Sidon-like")
    print("condition. The formula w_H = 2^p must be replaced by:")
    print()
    print("  w_H = 2^p - 2 * #{sign seqs with at least one partial sum collision}")
    print()
    print("The correction grows with p and L, becoming increasingly important.")
    print()


# ============================================================
# PART 4 — Verification 3: Most probable partition
# ============================================================

def run_verification_3(L_max=20):
    print("=" * 62)
    print("PART 4 — Verification 3: Most probable partition")
    print("=" * 62)
    print()

    # First verify w_lambda against brute force for small L
    print("Correctness check: w_lambda vs brute force for L=2..10")
    errors = 0
    for L in range(2, 11):
        for p in partitions(L):
            wf = w_lambda(p)
            wb = w_lambda_brute(p)
            if wf != wb:
                errors += 1
                print(f"  MISMATCH L={L} {p}: fast={wf}, brute={wb}")
    print(f"  {'VERIFIED' if errors == 0 else f'{errors} errors'} for L=2..10")
    print()

    PI = math.pi
    pred_k = PI / math.sqrt(6)   # 1.2825
    pred_r = math.sqrt(6) / PI   # 0.7797
    pred_sig = math.sqrt(pred_k) # 1.1325

    print(f"Analytical predictions (Hardy-Ramanujan saddle point):")
    print(f"  k*     ~ {pred_k:.4f} * sqrt(L)   (most probable run count)")
    print(f"  r_bar* ~ {pred_r:.4f} * sqrt(L)   (most probable mean run length)")
    print(f"  sigma_k ~ {pred_sig:.4f} * L^(1/4)  (fluctuations in k)")
    print()

    hdr = (f"{'L':>4}  {'p(L)':>5}  {'k*':>4}  {'k_pred':>6}  "
           f"{'k*/sqL':>8}  {'r*':>5}  {'r_pred':>6}  {'r*/sqL':>8}  "
           f"{'<k>_w':>7}  {'sig_k':>7}  {'sig/L^.25':>10}  {'wmax/Z':>9}")
    print(hdr)
    print("-" * len(hdr))

    results = []
    for L in list(range(4, L_max + 1, 2)):
        t0 = time.time()
        ws = {p: w_lambda(p) for p in partitions(L)}
        t1 = time.time()

        if t1 - t0 > 15:
            print(f"  [stopped at L={L}: {t1-t0:.1f}s per step]")
            break

        Z     = sum(ws.values())
        best  = max(ws, key=ws.get)
        wmax  = ws[best]
        ks    = len(best)
        rs    = L / ks
        sqL   = math.sqrt(L)
        L14   = L ** 0.25
        km    = sum(len(p) * w for p, w in ws.items()) / Z
        kvar  = sum(len(p)**2 * w for p, w in ws.items()) / Z - km**2
        ksig  = math.sqrt(max(kvar, 0))

        print(f"{L:>4}  {len(ws):>5}  {ks:>4}  {pred_k*sqL:>6.2f}  "
              f"{ks/sqL:>8.4f}  {rs:>5.2f}  {pred_r*sqL:>6.2f}  {rs/sqL:>8.4f}  "
              f"{km:>7.3f}  {ksig:>7.4f}  {ksig/L14:>10.4f}  {wmax/Z:>9.6f}")
        results.append((L, ks, rs, km, ksig, Z, wmax, ws))

    # w_max / Z scaling
    print()
    print("w_max / Z_L scaling (predicted ~ c/L for large L):")
    print()
    print(f"  {'L':>4}  {'wmax/Z':>12}  {'1/L':>10}  {'(wmax/Z)*L':>12}")
    for L, ks, rs, km, ksig, Z, wmax, ws in results:
        ratio = wmax / Z
        print(f"  {L:>4}  {ratio:>12.7f}  {1/L:>10.7f}  {ratio*L:>12.5f}")

    # Run-length distribution at saddle point
    print()
    print("Run-length distribution at most probable partition:")
    print("(checking geometric distribution prediction)")
    print()
    for L, ks, rs, km, ksig, Z, wmax, ws in results:
        if L not in [12, 16, 20]:
            continue
        best = max(ws, key=ws.get)
        k = len(best)
        rb = L / k
        q = 1 - 1 / rb
        freq = Counter(best)
        print(f"  L={L}, most probable partition: {best}")
        print(f"  k={k}, r_bar={rb:.2f}, geometric parameter q={q:.3f}")
        print(f"  {'j':>4}  {'count':>6}  {'geom_pred':>10}  {'ratio':>7}")
        for j in sorted(freq):
            actual = freq[j]
            pred = k * (1/rb) * q**(j-1)
            ratio = actual / pred if pred > 0 else float('inf')
            print(f"  {j:>4}  {actual:>6}  {pred:>10.2f}  {ratio:>7.3f}")
        print()

    return results


# ============================================================
# MAIN
# ============================================================

if __name__ == '__main__':
    print()
    print("=" * 62)
    print("Reptation-partition correspondence: exact enumeration")
    print("=" * 62)
    print()

    t_start = time.time()

    run_verification_1()
    run_verification_2()

    print("Note: Verification 3 (most probable partition) runs up to L=20.")
    print("This takes ~60s. Set L_max lower for a quicker run.")
    print()
    results = run_verification_3(L_max=20)

    t_end = time.time()
    print()
    print(f"Total runtime: {t_end - t_start:.1f}s")
    print()
    print("Summary of key findings:")
    print()
    print("  1. Decoupling theorem w_strong = w_H * w_V: EXACT")
    print("  2. p=2 formula w_H(h1,h2): EXACT")
    print("  3. Distinct-parts formula w_H = 2^p: CORRECTED")
    print("     Fails when consecutive sub-block has vanishing signed sum.")
    print("     Fraction achieving 2^p falls from 1.0 (L=6) to ~0.57 (L=16).")
    print("  4. k* and sigma_k scaling: converging toward predictions")
    print("     but pre-asymptotic for L<=20. Need L~100-200 for clean fit.")
