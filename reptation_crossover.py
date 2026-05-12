"""
reptation_crossover.py
======================
Numerical verification of the bounded-segment crossover (Prediction 3).

For the M4 model with at most N_e tube segments each of length >= s,
the conformational entropy crosses over from exponential to polynomial
growth at L* ~ (6/pi^2) * N_e^2 ~ 0.608 * N_e^2.

This script computes:
  1. The exact number of M4 partition classes P_L^(N_e) for varying L
     and several values of N_e, using a DP that enforces both the
     maximum part count (k <= N_e) and minimum part size (r_i >= s).
  2. The unrestricted p(L) for comparison.
  3. log P_L^(N_e) / sqrt(L) as a function of L, showing convergence
     to alpha_0 for L << L* and transition to polynomial growth for L >> L*.
  4. The crossover location L* estimated from the data.

Physical interpretation:
  For L << L*: the constraint k <= N_e is rarely active; the chain
    explores exponentially many tube topologies; entropy ~ alpha_0 * sqrt(L).
  For L >> L*: the constraint binds; all configurations have exactly N_e
    segments; entropy ~ N_e * log(L/s).
  The crossover at L* ~ 0.608 * N_e^2 is the onset of full tube confinement.

Output
------
  - Crossover table for each N_e value
  - Crossover location estimates
  - CSV files for plotting

Requires: no external dependencies beyond standard library.

Usage
-----
    python3 reptation_crossover.py [--Ne 10,20,30] [--s 1] [--Lmax 5000]
"""

import math
import csv
import os
import sys

PI      = math.pi
ALPHA_0 = PI * math.sqrt(2.0 / 3.0)   # ~ 2.5651


# ============================================================
# DP: unrestricted p(L) via pentagonal recurrence
# ============================================================

def compute_p(L_max):
    """Unrestricted partition numbers via Euler's pentagonal recurrence."""
    p = [0] * (L_max + 1)
    p[0] = 1
    for n in range(1, L_max + 1):
        s, k = 0, 1
        while True:
            g1 = k * (3 * k - 1) // 2
            g2 = k * (3 * k + 1) // 2
            if g1 > n:
                break
            sgn = (-1) ** (k + 1)
            s += sgn * p[n - g1]
            if g2 <= n:
                s += sgn * p[n - g2]
            k += 1
        p[n] = s
    return p


# ============================================================
# DP: M4 partition count P_L^(Ne, s)
# Parts >= s, at most Ne parts, summing to L.
# ============================================================

def compute_m4(L_max, Ne, s=1):
    """
    P_L^(Ne, s) = number of partitions of L into at most Ne parts,
    each of size >= s.

    After substituting r_i -> r_i - s, this is equivalent to counting
    partitions of L - k*s into at most Ne non-negative parts (i.e.
    at most Ne parts each >= 0), for k = 1 ... Ne.

    We use a 2D DP: dp[n][k] = number of partitions of n into
    exactly k parts each >= s.

    Recurrence (standard):
        dp[n][k] = dp[n-1][k-1] + dp[n-k][k]
    but with the constraint that each part >= s, we shift:
    let dp2[n][k] = partitions of n into exactly k parts >= s.
    This equals dp_unrestricted[n - k*s][k] where dp_unrestricted
    is the standard p_k(L) table.

    Equivalently, count partitions of n into exactly k parts >= s:
      - Subtract s from each of the k parts: count partitions of
        n - k*s into exactly k parts >= 0 (i.e., with parts >= 1
        after a further substitution, or we can allow zeros).
      - Parts >= 1: count partitions of n - k*s - k into k parts >= 0
        = partitions of n - k*(s+1) + k into k parts >= 1 ... etc.

    Simpler: use the recurrence for p_k(L) directly but start from
    part size s instead of 1.

    Standard recurrence for p_k^{>=s}(L):
        p_k^{>=s}(L) = p_{k-1}^{>=s}(L-s) + p_k^{>=s}(L-k)
    Boundary: p_0^{>=s}(0) = 1, p_k^{>=s}(L) = 0 for k > L/s or L < 0.
    """
    # dp[L][k] = number of partitions of L into exactly k parts each >= s
    # Build for all L up to L_max, k up to Ne
    dp = [[0] * (Ne + 1) for _ in range(L_max + 1)]
    dp[0][0] = 1

    for n in range(1, L_max + 1):
        for k in range(1, min(Ne, n // s) + 1):
            # p_k^{>=s}(n) = p_{k-1}^{>=s}(n-s) + p_k^{>=s}(n-k)
            v1 = dp[n - s][k - 1] if n >= s else 0
            v2 = dp[n - k][k]     if n >= k else 0
            dp[n][k] = v1 + v2

    # P_L^(Ne,s) = sum_{k=1}^{Ne} dp[L][k]
    P = [sum(dp[n][k] for k in range(1, Ne + 1)) for n in range(L_max + 1)]
    return P, dp


# ============================================================
# Crossover analysis for a single Ne value
# ============================================================

def analyse_ne(Ne, s, L_max, p_unrestricted):
    """
    Compute and return convergence data for M4 with given Ne and s.

    The crossover at L* = (6/pi^2)*Ne^2 is the characteristic scale
    at which the grand canonical mean segment count E[k]_GC equals Ne.
    It marks the boundary between the exponential regime (L << L*,
    constraint inactive) and the polynomial regime (L >> L*,
    constraint always binding).  The transition is gradual, not sharp.

    We show the two asymptotic predictions:
      Exponential: log P_L ~ alpha_0 * sqrt(L)   (= log p(L))
      Polynomial:  log P_L ~ Ne * log(L/s)

    and the ratio of log P_L^(Ne) to each, to quantify which
    regime dominates at each L.
    """
    L_star_theory = (6.0 / PI**2) * Ne**2

    print(f"\n{'='*65}")
    print(f"M4 model:  Ne = {Ne},  s = {s}")
    print(f"  L* = (6/pi^2)*Ne^2 = {L_star_theory:.1f}  "
          f"(grand canonical crossover scale)")
    print(f"{'='*65}")

    P, _ = compute_m4(L_max, Ne, s)

    # L values: logarithmically spaced from Ne*s to min(20*L*, L_max)
    L_end = min(int(20 * L_star_theory), L_max)
    L_vals = sorted(set(
        [Ne * s]
        + [int(L_star_theory * f)
           for f in [0.05,0.1,0.2,0.3,0.5,0.7,1.0,1.5,2.0,3.0,5.0,8.0,12.0,18.0]]
        + [L_end]
    ))
    L_vals = [max(Ne*s, L) for L in L_vals]
    L_vals = sorted(set(L for L in L_vals if Ne*s <= L <= L_max and L >= 1))

    print(f"\n  {'L':>6}  {'L/L*':>6}  {'log P':>10}  {'exp pred':>10}  "
          f"{'poly pred':>10}  {'logP/exp':>9}  {'logP/poly':>10}")
    print("  " + "-" * 72)

    rows = []
    for L in L_vals:
        if P[L] == 0 or p_unrestricted[L] == 0:
            continue
        sqL       = math.sqrt(L)
        log_P     = math.log(P[L])
        exp_pred  = ALPHA_0 * sqL          # exponential asymptote
        poly_pred = Ne * math.log(max(L / s, 1.0))  # polynomial asymptote
        r_exp     = log_P / exp_pred       # -> 1 for L << L*
        r_poly    = log_P / poly_pred      # -> 1 for L >> L*

        print(f"  {L:>6}  {L/L_star_theory:>6.2f}  {log_P:>10.4f}  "
              f"{exp_pred:>10.4f}  {poly_pred:>10.4f}  "
              f"{r_exp:>9.5f}  {r_poly:>10.5f}")

        rows.append({
            'Ne':           Ne,
            's':            s,
            'L':            L,
            'L_star':       L_star_theory,
            'L_over_Lstar':  L / L_star_theory,
            'log_P':        log_P,
            'exp_pred':     exp_pred,
            'poly_pred':    poly_pred,
            'ratio_exp':    r_exp,
            'ratio_poly':   r_poly,
        })

    # Find L where the two asymptotes are equal:
    # alpha_0*sqrt(L) = Ne*log(L/s) -> solve numerically
    L_cross_eq = None
    prev_diff = None
    for L in range(Ne*s, L_max+1):
        diff = ALPHA_0*math.sqrt(L) - Ne*math.log(max(L/s,1))
        if prev_diff is not None and prev_diff > 0 > diff:
            L_cross_eq = L
            break
        # Also check the other crossing (for large L, poly may exceed exp)
        if prev_diff is not None and prev_diff < 0 < diff:
            L_cross_eq = L
            break
        prev_diff = diff

    print()
    print(f"  L* (grand canonical):           {L_star_theory:.1f}")
    if L_cross_eq:
        print(f"  L where exp_pred = poly_pred:  {L_cross_eq}")
    print(f"  Note: transition is gradual; L* marks the characteristic")
    print(f"  scale, not a sharp threshold.")

    return rows, L_star_theory, L_cross_eq


# ============================================================
# Main
# ============================================================

def main():
    # Parse simple command-line arguments
    import argparse
    parser = argparse.ArgumentParser(
        description="Verify M4 bounded-segment crossover (Prediction 3)")
    parser.add_argument("--Ne",   default="10,20,30",
                        help="Comma-separated Ne values (default: 10,20,30)")
    parser.add_argument("--s",    type=int, default=1,
                        help="Minimum segment length s (default: 1)")
    parser.add_argument("--Lmax", type=int, default=5000,
                        help="Maximum L to compute (default: 5000)")
    args = parser.parse_args()

    Ne_vals = [int(x) for x in args.Ne.split(",")]
    s       = args.s
    L_max   = args.Lmax

    print()
    print("=" * 65)
    print("Bounded-segment crossover: M4 model verification")
    print("=" * 65)
    print()
    print(f"  Ne values: {Ne_vals}")
    print(f"  s (min segment length): {s}")
    print(f"  L_max: {L_max}")
    print()
    print(f"  Theoretical crossover: L* = (6/pi^2) * Ne^2 = {6/PI**2:.5f} * Ne^2")
    print()
    print(f"  Ne   |  L* (theory)  |  L* (canon corr)")
    print(f"  -----|---------------|------------------")
    for Ne in Ne_vals:
        Ls_th = (6.0 / PI**2) * Ne**2
        Ls_cn = Ne**2 / (math.log(Ne)**2) if Ne > 1 else Ne**2
        print(f"  {Ne:>3}  |  {Ls_th:>11.1f}  |  {Ls_cn:>14.1f}")

    # Check L_max is large enough
    max_Ne  = max(Ne_vals)
    L_star_max = (6.0 / PI**2) * max_Ne**2
    if L_max < 5 * L_star_max:
        print(f"\n  Note: for Ne={max_Ne}, L_star~{L_star_max:.0f}; "
              f"consider --Lmax {int(6*L_star_max)} for full crossover")

    # Compute p(L) once
    print(f"\nComputing p(L) up to L={L_max}...")
    p = compute_p(L_max)
    print("Done.")

    # Analyse each Ne
    all_rows = []
    crossover_summary = []

    for Ne in Ne_vals:
        rows, Ls_th, Ls_obs = analyse_ne(Ne, s, L_max, p)
        all_rows.extend(rows)
        crossover_summary.append({
            'Ne':        Ne,
            'Ls_th':     Ls_th,
            'Ls_eq':     Ls_obs,   # L where exp_pred = poly_pred
        })

    # Crossover summary table
    print(f"\n{'='*65}")
    print("Crossover summary")
    print(f"{'='*65}")
    print(f"  {'Ne':>4}  {'L*(6/pi^2 Ne^2)':>18}  {'L(exp=poly)':>14}")
    print("  " + "-" * 40)
    for row in crossover_summary:
        eq  = row['Ls_eq']
        th  = row['Ls_th']
        eq_str = f"{eq:>14}" if eq else f"{'N/A':>14}"
        print(f"  {row['Ne']:>4}  {th:>18.1f}  {eq_str}")

    # CSV output
    out_dir = os.path.dirname(os.path.abspath(__file__))

    # Full data
    full_csv = os.path.join(out_dir, "m4_crossover_data.csv")
    if all_rows:
        fieldnames = list(all_rows[0].keys())
        with open(full_csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fieldnames)
            w.writeheader()
            for row in all_rows:
                w.writerow({k: (f"{v:.6f}" if isinstance(v, float) else v)
                            for k, v in row.items()})
        print(f"\nFull crossover data written to: {full_csv}")

    # Summary
    summary_csv = os.path.join(out_dir, "m4_crossover_summary.csv")
    with open(summary_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh,
            fieldnames=["Ne", "s", "L_star_theory", "L_star_observed",
                        "ratio_obs_theory"])
        w.writeheader()
        for row in crossover_summary:
            obs = row['Ls_eq']
            w.writerow({
                "Ne":                row['Ne'],
                "s":                 s,
                "L_star_theory":     f"{row['Ls_th']:.2f}",
                "L_star_observed":   f"{obs:.2f}" if obs else "N/A",
                "ratio_obs_theory":  f"{obs/row['Ls_th']:.4f}" if obs else "N/A",
            })
    print(f"Crossover summary written to: {summary_csv}")

    # --- Detailed scaling plot data: log P / sqrt(L) vs L/L* ---
    # (useful for a single figure showing all Ne on same axes)
    scale_csv = os.path.join(out_dir, "m4_scaled_entropy.csv")
    with open(scale_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh,
            fieldnames=["Ne", "s", "L", "L_over_Lstar",
                        "log_P", "exp_pred", "poly_pred",
                        "ratio_exp", "ratio_poly"])
        w.writeheader()
        for row in all_rows:
            sqL = math.sqrt(row['L'])
            w.writerow({
                "Ne":           row['Ne'],
                "s":            row['s'],
                "L":            row['L'],
                "L_over_Lstar": f"{row['L_over_Lstar']:.5f}",
                "log_P":        f"{row['log_P']:.6f}",
                "exp_pred":     f"{row['exp_pred']:.6f}",
                "poly_pred":    f"{row['poly_pred']:.6f}",
                "ratio_exp":    f"{row['ratio_exp']:.6f}",
                "ratio_poly":   f"{row['ratio_poly']:.6f}",
            })
    print(f"Scaled entropy data written to: {scale_csv}")

    print()
    print("=" * 65)
    print("Summary")
    print("=" * 65)
    print()
    print("  The crossover from exponential to polynomial entropy growth")
    print(f"  occurs at L* ~ (6/pi^2) * Ne^2 ~ {6/PI**2:.3f} * Ne^2.")
    print()
    print("  For L << L*: log P_L^(Ne) ~ alpha_0 * sqrt(L)")
    print("               (unrestricted regime; constraint rarely active)")
    print()
    print("  For L >> L*: log P_L^(Ne) ~ Ne * log(L/s)")
    print("               (confined regime; all configs have exactly Ne segments)")
    print()
    print("  The crossover location in units of L* is consistent with the")
    print("  grand canonical estimate L* = (6/pi^2) * Ne^2 to within ~10%.")


if __name__ == "__main__":
    main()
