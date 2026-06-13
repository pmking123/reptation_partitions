"""
aggregate.py
============
Collect JSON result files from SLURM job-array tasks and print the
same formatted tables that the standalone (sequential) scripts produce.

Usage
-----
After mc_slurm.sh completes:
    python3 aggregate.py --mode weighted --indir results_mc

After mc_uniform_slurm.sh completes:
    python3 aggregate.py --mode uniform  --indir results_mc_uniform

Options
-------
--mode    weighted | uniform
--indir   directory containing the JSON result files (default: auto)
--paper   if given, print a side-by-side comparison with the paper's
          quoted KL divergence table (uniform mode only)
"""

import argparse
import json
import math
import os
import sys


PI = math.pi


# ============================================================
# Weighted-measure scaling study table
# ============================================================

def aggregate_weighted(indir):
    from mc_fast import SCALING_L_VALS

    print()
    print("=" * 62)
    print("Reptation partition MC sampler — scaling study results")
    print("(weighted measure  w_tilde = k! * 2^k / prod m_j!)")
    print("=" * 62)
    print()

    c_EL = math.sqrt(6) / (2 * PI)
    print(f"  Erdos-Lehner mean scaling: E[k] ~ {c_EL:.4f} * sqrt(L) * log(L)")
    print()

    hdr = (f"{'L':>6}  {'<k>/sqL':>9}  {'<k>/(sqL*lnL)':>14}  "
           f"{'sig/L^.25':>10}  {'tau':>7}  {'acc':>6}  {'t(s)':>5}")
    print(hdr)
    print("-" * len(hdr))

    missing = []
    for L in SCALING_L_VALS:
        path = os.path.join(indir, f"scaling_{L}.json")
        if not os.path.exists(path):
            missing.append(L)
            print(f"  {'(missing)':>6}  L={L}")
            continue
        with open(path) as f:
            r = json.load(f)
        print(f"{r['L']:>6}  {r['km_over_sqL']:>9.4f}  "
              f"{r['ratio_logL']:>14.5f}  "
              f"{r['ksig_over_L14']:>10.4f}  "
              f"{r['tau']:>7.1f}  {r['acc']:>6.3f}  {r['elapsed']:>5.1f}")

    if missing:
        print()
        print(f"  WARNING: missing results for L = {missing}")
    print()


# ============================================================
# Uniform-measure mode convergence table
# ============================================================

def aggregate_mode(indir):
    from mc_uniform_fast import MODE_L_VALS

    c_EL = math.sqrt(6) / (2 * PI)

    print()
    print("=" * 65)
    print("Mode/mean convergence study — uniform measure")
    print("=" * 65)
    print()
    print(f"  Erdos-Lehner: E[k]_canon ~ {c_EL:.4f} * sqrt(L) * log(L)")
    print(f"  Modal k (HR saddle): modal_k ~ {PI/math.sqrt(6):.4f} * sqrt(L)")
    print()

    hdr = (f"{'L':>7}  {'modal_k':>8}  {'modal/sqL':>10}  "
           f"{'k*/(sqL*lnL)':>13}  {'<k>/sqL':>9}  {'sig/L^.25':>10}  "
           f"{'tau':>7}  {'t(s)':>5}")
    print(hdr)
    print("-" * len(hdr))

    missing = []
    for L in MODE_L_VALS:
        path = os.path.join(indir, f"mode_{L}.json")
        if not os.path.exists(path):
            missing.append(L)
            print(f"  {'(missing)':>7}  L={L}")
            continue
        with open(path) as f:
            r = json.load(f)
        print(f"{r['L']:>7}  {r['modal_k']:>8}  "
              f"{r['modal_over_sqL']:>10.4f}  "
              f"{r['ratio_logL']:>13.5f}  "
              f"{r['km_over_sqL']:>9.4f}  "
              f"{r['ksig_over_L14']:>10.4f}  "
              f"{r['tau']:>7.1f}  {r['elapsed']:>5.1f}")

    if missing:
        print()
        print(f"  WARNING: missing results for L = {missing}")
    print()


# ============================================================
# Uniform-measure run-length distribution tables + KL summary
# ============================================================

# Paper's quoted KL divergences (Section 8.3 table)
PAPER_KL = {
    100: (0.040, 0.072),
    200: (0.023, 0.099),
    500: (0.026, 0.145),
}

def aggregate_rld(indir, compare_paper=False):
    from mc_uniform_fast import RLD_L_VALS

    print()
    print("=" * 65)
    print("Run-length distributions: BE vs geometric")
    print("(uniform measure, paper Section 8.3)")
    print("=" * 65)
    print()

    kl_rows = []
    missing = []

    for L in RLD_L_VALS:
        path = os.path.join(indir, f"rld_{L}.json")
        if not os.path.exists(path):
            missing.append(L)
            continue
        with open(path) as f:
            r = json.load(f)

        print(f"  L={L}: <k>_MC={r['k_mean']:.2f}  r_bar={r['r_bar']:.2f}"
              f"  (GC pred: <k>={r['pred_k_gc']:.2f}, "
              f"r_bar={r['pred_r_gc']:.2f})")
        print(f"  alpha = pi/sqrt(6L) = {r['alpha']:.5f}")
        print(f"  D_KL(MC||BE)={r['kl_be']:.4f}  "
              f"D_KL(MC||geom)={r['kl_geom']:.4f}  "
              f"ratio={r['kl_ratio']:.2f}")
        print()
        print(f"  {'j':>4}  {'f_MC':>10}  {'f_BE':>10}  "
              f"{'ratio_BE':>9}  {'f_geom':>10}  {'ratio_geom':>11}")

        f_mc   = {int(j): v for j, v in r['f_mc'].items()}
        f_be   = {int(j): v for j, v in r['f_be'].items()}
        f_geom = {int(j): v for j, v in r['f_geom'].items()}
        cutoff = int(6 * r['r_bar']) + 2

        for j in range(1, cutoff):
            fm  = f_mc.get(j, 0.0)
            fbe = f_be.get(j, 0.0)
            fg  = f_geom.get(j, 0.0)
            if fm < 5e-5 and fbe < 5e-5:
                break
            r_be  = fm / fbe if fbe  > 1e-9 else 0.0
            r_geo = fm / fg  if fg   > 1e-9 else 0.0
            print(f"  {j:>4}  {fm:>10.5f}  {fbe:>10.5f}  "
                  f"{r_be:>9.4f}  {fg:>10.5f}  {r_geo:>11.4f}")
        print()
        kl_rows.append((L, r['kl_be'], r['kl_geom'], r['kl_ratio']))

    if missing:
        print(f"  WARNING: missing RLD results for L = {missing}")
        print()

    # KL summary table (mirrors the paper's table)
    if kl_rows:
        print("KL divergence summary")
        print()
        if compare_paper:
            print(f"  {'L':>5}  {'KL_BE (MC)':>12}  {'KL_BE (paper)':>14}  "
                  f"{'KL_geom (MC)':>13}  {'KL_geom (paper)':>16}  "
                  f"{'ratio (MC)':>11}  {'ratio (paper)':>14}")
            print("  " + "-" * 95)
            for L, kl_be, kl_geom, ratio in kl_rows:
                p_be, p_geo = PAPER_KL.get(L, (None, None))
                p_ratio = p_geo / p_be if (p_be and p_geo) else None
                be_str  = f"{p_be:.3f}" if p_be   else "  n/a"
                geo_str = f"{p_geo:.3f}" if p_geo  else "  n/a"
                rat_str = f"{p_ratio:.1f}" if p_ratio else " n/a"
                print(f"  {L:>5}  {kl_be:>12.4f}  {be_str:>14}  "
                      f"{kl_geom:>13.4f}  {geo_str:>16}  "
                      f"{ratio:>11.2f}  {rat_str:>14}")
        else:
            print(f"  {'L':>5}  {'D_KL(MC||BE)':>14}  "
                  f"{'D_KL(MC||geom)':>16}  {'ratio':>7}")
            print("  " + "-" * 48)
            for L, kl_be, kl_geom, ratio in kl_rows:
                print(f"  {L:>5}  {kl_be:>14.4f}  {kl_geom:>16.4f}  "
                      f"{ratio:>7.2f}")
        print()
        print("  Expected: ratio grows with L (BE improves relative to")
        print("  geometric as chain lengthens — paper Section 8.3).")
        print()


# ============================================================
# Cross-check: weighted vs uniform <k> at matching L values
# ============================================================

def cross_check(weighted_indir, uniform_indir):
    from mc_fast        import SCALING_L_VALS as W_VALS
    from mc_uniform_fast import MODE_L_VALS   as U_VALS

    shared = sorted(set(W_VALS) & set(U_VALS))
    if not shared:
        return

    print()
    print("=" * 65)
    print("Cross-check: weighted vs uniform <k> at shared L values")
    print("(weighted <k> should be >> uniform <k>, per Section 6.3)")
    print("=" * 65)
    print()
    print(f"  {'L':>6}  {'<k>_weighted':>14}  {'<k>_uniform':>13}  {'ratio':>7}")
    print("  " + "-" * 46)

    for L in shared:
        wp = os.path.join(weighted_indir, f"scaling_{L}.json")
        up = os.path.join(uniform_indir,  f"mode_{L}.json")
        if not os.path.exists(wp) or not os.path.exists(up):
            continue
        with open(wp) as f: wr = json.load(f)
        with open(up) as f: ur = json.load(f)
        km_w = wr['km']
        km_u = ur['km']
        print(f"  {L:>6}  {km_w:>14.2f}  {km_u:>13.2f}  "
              f"{km_w / km_u:>7.2f}")
    print()


# ============================================================
# Entry point
# ============================================================

def main():
    parser = argparse.ArgumentParser(description="Aggregate SLURM JSON output")
    parser.add_argument('--mode',    choices=['weighted', 'uniform', 'both'],
                        default='both')
    parser.add_argument('--indir',   type=str, default=None,
                        help='Result directory (overrides defaults)')
    parser.add_argument('--paper',   action='store_true',
                        help='Compare KL table against paper values')
    parser.add_argument('--weighted-indir', type=str,
                        default='results_mc',
                        help='Weighted results dir (--mode both)')
    parser.add_argument('--uniform-indir',  type=str,
                        default='results_mc_uniform',
                        help='Uniform results dir (--mode both)')
    args = parser.parse_args()

    if args.mode == 'weighted':
        indir = args.indir or args.weighted_indir
        aggregate_weighted(indir)

    elif args.mode == 'uniform':
        indir = args.indir or args.uniform_indir
        aggregate_mode(indir)
        aggregate_rld(indir, compare_paper=args.paper)

    else:  # both
        w_indir = args.indir or args.weighted_indir
        u_indir = args.indir or args.uniform_indir
        aggregate_weighted(w_indir)
        aggregate_mode(u_indir)
        aggregate_rld(u_indir, compare_paper=args.paper)
        cross_check(w_indir, u_indir)


if __name__ == '__main__':
    main()
