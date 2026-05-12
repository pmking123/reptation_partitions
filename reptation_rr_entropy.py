"""
reptation_rr_entropy.py
=======================
Numerical verification of the Rogers-Ramanujan entropy coefficient
and the M5 distinct-part growth rate.

Two asymptotic results are verified:

RESULT 1: Rogers-Ramanujan (M6, gap-2 partitions)
  The number of partitions of L whose parts differ by >= 2 (gap-2
  partitions) satisfies:
      log q_RR(L) ~ alpha_RR * sqrt(L),   alpha_RR = 2*pi/sqrt(15)
  This is derived in the paper via the product side of the first
  Rogers-Ramanujan identity and the Hardy-Ramanujan circle method.
  Here we verify it by direct exact enumeration up to L=2000.

  The ratio alpha_RR / alpha_0 = sqrt(2/5) ~ 0.6325, corresponding
  to a 37% reduction in conformational entropy relative to the
  unrestricted (M0) case.

RESULT 2: M5 distinct-part partitions
  The number of partitions of L into distinct parts satisfies:
      log d(L) ~ alpha_5 * sqrt(L),   alpha_5 = pi/sqrt(12) = alpha_0/sqrt(2)
  This is the Euler identity result for strict partitions.  The growth
  rate is smaller than alpha_0 by a factor of sqrt(2).

Both are computed by dynamic programming with O(L^2) time and space.

Output
------
  - Convergence table: log q_RR(L)/sqrt(L) -> alpha_RR
  - Convergence table: log d(L)/sqrt(L)    -> alpha_5
  - Ratio table:       alpha_RR / alpha_0 and alpha_5 / alpha_0
  - CSV files for plotting

Requires: no external dependencies beyond standard library.

Usage
-----
    python3 reptation_rr_entropy.py
"""

import math
import csv
import os

PI   = math.pi
# Asymptotic coefficients
ALPHA_0  = PI * math.sqrt(2.0 / 3.0)          # ~ 2.5651  unrestricted
ALPHA_RR = 2.0 * PI / math.sqrt(15.0)          # ~ 1.6223  gap-2 (RR)
ALPHA_5  = PI / math.sqrt(3.0)           # ~ 1.8138  distinct parts = alpha_0/sqrt(2)
# Ratios
RATIO_RR = ALPHA_RR / ALPHA_0                  # = sqrt(2/5) ~ 0.6325
RATIO_5  = ALPHA_5  / ALPHA_0                  # = 1/sqrt(2) ~ 0.7071


# ============================================================
# DP: unrestricted partitions p(L)
# ============================================================

def compute_p(L_max):
    """
    p(L) = number of unrestricted partitions of L.
    Uses the standard recurrence via the pentagonal number theorem,
    or equivalently the triangle dp[L][k] summed over k.
    Here we use the simpler one-dimensional recurrence:
        p(L) = sum_{k>=1} (-1)^{k+1} [p(L - k(3k-1)/2) + p(L - k(3k+1)/2)]
    which follows from Euler's pentagonal theorem.
    """
    p = [0] * (L_max + 1)
    p[0] = 1
    for n in range(1, L_max + 1):
        s = 0
        k = 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > n:
                break
            sign = (-1) ** (k + 1)
            s += sign * p[n - g1]
            if g2 <= n:
                s += sign * p[n - g2]
            k += 1
        p[n] = s
    return p


# ============================================================
# DP: gap-2 partitions q_RR(L)
# Rogers-Ramanujan class: all parts differ by >= 2
# ============================================================

def compute_q_rr(L_max):
    """
    q_RR(L) = number of partitions of L with the property that
    consecutive parts (in non-increasing order) differ by >= 2.
    Equivalently: partitions into parts that are congruent to
    +-1 mod 5 (by the first Rogers-Ramanujan identity), but we
    compute directly via DP.

    Recurrence: let f(n, m) = number of gap-2 partitions of n
    whose largest part is <= m.  Then:
        f(n, m) = f(n, m-1) + f(n-m, m-2)
    with f(0, m) = 1 for all m >= 0, f(n, m) = 0 for n < 0 or m < 1.

    We compute q_RR(L) = f(L, L) for each L.
    """
    # dp[n][m] = number of gap-2 partitions of n with largest part <= m
    # We need dp[L][L] for each L up to L_max.
    # To save memory, note that f(n, m) = f(n, m-1) + f(n-m, m-2).
    # We build a 2D table of size (L_max+1) x (L_max+3).
    L = L_max
    # f[n][m]: n up to L, m up to L+2 (to handle m-2 without underflow)
    f = [[0] * (L + 3) for _ in range(L + 1)]
    for m in range(L + 3):
        f[0][m] = 1
    for n in range(1, L + 1):
        for m in range(1, L + 2):
            f[n][m] = f[n][m - 1]
            if n >= m and m >= 2:
                f[n][m] += f[n - m][m - 2]
            # if m == 1: the only gap-2 partition using only 1s is a single 1
            # (since two 1s would have gap 0). So f(1,1)=1, f(n,1)=0 for n>1.
            elif n >= m and m == 1:
                f[n][m] += (1 if n == 1 else 0)
    q = [f[n][n] for n in range(L + 1)]
    return q


# ============================================================
# DP: distinct-part partitions d(L)
# ============================================================

def compute_d(L_max):
    """
    d(L) = number of partitions of L into distinct parts.
    Equivalently: strict partitions (all parts different).
    Uses the recurrence:
        d(L) = sum over subsets S of {1,...,L} with sum(S) = L of 1
    implemented via DP:
        dp[j] = number of distinct-part partitions of j
    building up part by part.
    """
    dp = [0] * (L_max + 1)
    dp[0] = 1
    for part in range(1, L_max + 1):
        # Traverse in reverse to ensure each part is used at most once
        for j in range(L_max, part - 1, -1):
            dp[j] += dp[j - part]
    return dp


# ============================================================
# Convergence analysis
# ============================================================

def convergence_table(vals, L_vals, name, alpha_theory, alpha_0):
    """
    Print and return convergence data: log(vals[L])/sqrt(L) -> alpha_theory.
    """
    rows = []
    print(f"\n{'='*65}")
    print(f"Convergence of log {name}(L) / sqrt(L)  ->  {alpha_theory:.6f}")
    print(f"{'='*65}")
    print(f"  {'L':>6}  {'log count':>12}  {'ratio':>10}  "
          f"{'theory':>10}  {'% error':>9}  {'ratio/alpha_0':>14}")
    print("  " + "-" * 65)

    for L in L_vals:
        if L < 2 or vals[L] <= 0:
            continue
        sqL    = math.sqrt(L)
        log_v  = math.log(vals[L])
        ratio  = log_v / sqL
        err_pc = 100.0 * (ratio - alpha_theory) / alpha_theory
        r_a0   = ratio / alpha_0
        print(f"  {L:>6}  {log_v:>12.4f}  {ratio:>10.6f}  "
              f"{alpha_theory:>10.6f}  {err_pc:>+9.4f}  {r_a0:>14.6f}")
        rows.append({'L': L, 'log_count': log_v, 'ratio': ratio,
                     'theory': alpha_theory, 'pct_error': err_pc,
                     'ratio_over_alpha0': r_a0})
    return rows


# ============================================================
# Main
# ============================================================

def main():
    L_max = 2000

    print()
    print("=" * 65)
    print("Rogers-Ramanujan and M5 entropy coefficients")
    print("Numerical verification via exact DP")
    print("=" * 65)
    print()
    print(f"  alpha_0  = pi*sqrt(2/3)   = {ALPHA_0:.6f}  (unrestricted, M0-M3)")
    print(f"  alpha_RR = 2*pi/sqrt(15)  = {ALPHA_RR:.6f}  (gap-2, M6)")
    print(f"  alpha_5  = pi/sqrt(3)     = {ALPHA_5:.6f}  (distinct parts, M5 = alpha_0/sqrt(2))")
    print()
    print(f"  alpha_RR / alpha_0 = sqrt(2/5) = {RATIO_RR:.6f}"
          f"  ({100*(1-RATIO_RR):.1f}% reduction)")
    print(f"  alpha_5  / alpha_0 = 1/sqrt(2) = {RATIO_5:.6f}"
          f"  ({100*(1-RATIO_5):.1f}% entropy reduction vs unrestricted)")

    # --- Compute counts ---
    print(f"\nComputing DP tables up to L = {L_max} ...")
    p   = compute_p(L_max)
    q   = compute_q_rr(L_max)
    d   = compute_d(L_max)
    print("Done.\n")

    # Spot-check against known values
    # p(10)=42, p(20)=627, p(50)=204226
    # d(10)=10, d(20)=64, d(50)=3658
    # q_RR: first few values are 1,1,1,1,2,2,3,3,4,5,6,7,9,10,12,14,17,19,23,...
    print("Spot checks:")
    print(f"  p(10)={p[10]} (expected 42),  "
          f"p(20)={p[20]} (expected 627),  "
          f"p(50)={p[50]} (expected 204226)")
    print(f"  d(10)={d[10]} (expected 10),  "
          f"d(20)={d[20]} (expected 64),   "
          f"d(50)={d[50]} (expected 3658)")
    # q_RR(10): partitions of 10 with gap>=2
    # They are: 10, 9+1, 8+2, 7+3, 6+4, 6+3+1, 5+4+1 -> 7? Let's just print
    print(f"  q_RR(10)={q[10]},  q_RR(20)={q[20]},  q_RR(50)={q[50]}")

    # --- L values for convergence tables ---
    L_vals = (
        [5, 10, 20, 30, 50, 75, 100, 150, 200, 300, 400,
         500, 750, 1000, 1250, 1500, 1750, 2000]
    )
    L_vals = [L for L in L_vals if L <= L_max]

    # --- Convergence tables ---
    rows_rr = convergence_table(q, L_vals, "q_RR", ALPHA_RR, ALPHA_0)
    rows_d  = convergence_table(d, L_vals, "d",    ALPHA_5,  ALPHA_0)  # ALPHA_5=pi/sqrt(3)

    # Also show unrestricted p(L) as reference
    rows_p  = convergence_table(p, L_vals, "p",    ALPHA_0,  ALPHA_0)

    # --- Combined ratio table ---
    print(f"\n{'='*65}")
    print("Combined: entropy coefficients vs L")
    print(f"{'='*65}")
    print(f"  {'L':>6}  {'alpha_p':>9}  {'alpha_RR':>10}  "
          f"{'alpha_d':>9}  {'RR/p':>8}  {'d/p':>8}")
    print("  " + "-" * 55)

    p_map  = {r['L']: r['ratio'] for r in rows_p}
    rr_map = {r['L']: r['ratio'] for r in rows_rr}
    d_map  = {r['L']: r['ratio'] for r in rows_d}

    for L in L_vals:
        if L not in p_map or L not in rr_map or L not in d_map:
            continue
        ap  = p_map[L]
        arr = rr_map[L]
        ad  = d_map[L]
        print(f"  {L:>6}  {ap:>9.6f}  {arr:>10.6f}  "
              f"{ad:>9.6f}  {arr/ap:>8.6f}  {ad/ap:>8.6f}")
    print()
    print(f"  Theoretical limits:   {ALPHA_0:.6f}   "
          f"{ALPHA_RR:.6f}   {ALPHA_5:.6f}   "
          f"{RATIO_RR:.6f}   {RATIO_5:.6f}")

    # --- CSV output ---
    out_dir = os.path.dirname(os.path.abspath(__file__))

    # RR convergence
    rr_csv = os.path.join(out_dir, "rr_entropy_convergence.csv")
    with open(rr_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "L", "log_q_RR", "alpha_eff", "alpha_RR_theory",
            "pct_error", "ratio_over_alpha0"])
        w.writeheader()
        for r in rows_rr:
            w.writerow({
                "L":               r["L"],
                "log_q_RR":        f"{r['log_count']:.6f}",
                "alpha_eff":       f"{r['ratio']:.6f}",
                "alpha_RR_theory": f"{ALPHA_RR:.6f}",
                "pct_error":       f"{r['pct_error']:.4f}",
                "ratio_over_alpha0": f"{r['ratio_over_alpha0']:.6f}",
            })
    print(f"RR convergence data written to: {rr_csv}")

    # Distinct-part convergence
    d_csv = os.path.join(out_dir, "distinct_part_convergence.csv")
    with open(d_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "L", "log_d", "alpha_eff", "alpha_5_theory",
            "pct_error", "ratio_over_alpha0"])
        w.writeheader()
        for r in rows_d:
            w.writerow({
                "L":             r["L"],
                "log_d":         f"{r['log_count']:.6f}",
                "alpha_eff":     f"{r['ratio']:.6f}",
                "alpha_5_theory": f"{ALPHA_5:.6f}",
                "pct_error":     f"{r['pct_error']:.4f}",
                "ratio_over_alpha0": f"{r['ratio_over_alpha0']:.6f}",
            })
    print(f"Distinct-part convergence data written to: {d_csv}")

    # Combined table CSV (all three models)
    combined_csv = os.path.join(out_dir, "entropy_coefficients.csv")
    with open(combined_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=[
            "L", "log_p", "log_q_RR", "log_d",
            "alpha_p", "alpha_RR", "alpha_d",
            "alpha_RR_over_alpha_p", "alpha_d_over_alpha_p"])
        w.writeheader()
        for L in L_vals:
            if L not in p_map or L not in rr_map or L not in d_map:
                continue
            w.writerow({
                "L":                     L,
                "log_p":                 f"{math.log(p[L]):.6f}" if p[L]>0 else "",
                "log_q_RR":              f"{math.log(q[L]):.6f}" if q[L]>0 else "",
                "log_d":                 f"{math.log(d[L]):.6f}" if d[L]>0 else "",
                "alpha_p":               f"{p_map[L]:.6f}",
                "alpha_RR":              f"{rr_map[L]:.6f}",
                "alpha_d":               f"{d_map[L]:.6f}",
                "alpha_RR_over_alpha_p": f"{rr_map[L]/p_map[L]:.6f}",
                "alpha_d_over_alpha_p":  f"{d_map[L]/p_map[L]:.6f}",
            })
    print(f"Combined entropy coefficients written to: {combined_csv}")

    # --- Summary ---
    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    print()
    # Use L=2000 for final estimates
    L_ref = 2000
    if p[L_ref] > 0 and q[L_ref] > 0 and d[L_ref] > 0:
        sqL = math.sqrt(L_ref)
        ap  = math.log(p[L_ref]) / sqL
        arr = math.log(q[L_ref]) / sqL
        ad  = math.log(d[L_ref]) / sqL
        print(f"  At L={L_ref}:")
        print(f"    alpha_p  (unrestricted)  = {ap:.5f}  "
              f"(theory {ALPHA_0:.5f}, err {100*(ap-ALPHA_0)/ALPHA_0:+.2f}%)")
        print(f"    alpha_RR (gap-2)         = {arr:.5f}  "
              f"(theory {ALPHA_RR:.5f}, err {100*(arr-ALPHA_RR)/ALPHA_RR:+.2f}%)")
        print(f"    alpha_d  (distinct)      = {ad:.5f}  "
              f"(theory {ALPHA_5:.5f}, err {100*(ad-ALPHA_5)/ALPHA_5:+.2f}%)")
        print()
        print(f"    Observed ratio alpha_RR/alpha_p = {arr/ap:.5f}  "
              f"(theory sqrt(2/5) = {RATIO_RR:.5f})")
        print(f"    Observed ratio alpha_d/alpha_p  = {ad/ap:.5f}  "
              f"(theory 1/sqrt(2) = {RATIO_5:.5f})")
        print()
        print(f"    Entropy reduction (RR):      {100*(1-arr/ap):.1f}%  "
              f"(theory 37.0%)")
        print(f"    Entropy reduction (distinct): {100*(1-ad/ap):.1f}%  "
              f"(theory 29.3%)")


if __name__ == "__main__":
    main()
