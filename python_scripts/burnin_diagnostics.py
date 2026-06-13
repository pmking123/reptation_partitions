"""
burnin_diagnostics.py
=====================
Four burn-in sufficiency checks for mc_fast.py at a given L.

Usage
-----
    python3 burnin_diagnostics.py              # defaults: L=500, burnin=100000
    python3 burnin_diagnostics.py --L 500 --burnin 100000 --chains 4

Output
------
  burnin_trace_L{L}.png        Trace plot of k for first `burnin` steps (check 1)
  burnin_runmean_L{L}.png      Running mean of k over burn-in (check 4)
  burnin_convergence_L{L}.txt  Gelman-Rubin R-hat and starting-point summary (checks 2, 3)
"""

import argparse
import math
import random
import sys
import os
import time

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Allow running from any directory that contains partition_numpy.py / mc_fast.py
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from partition_numpy import init_arrays, arrays_to_tuple
from mc_fast import metropolis_step_fast


# ============================================================
# Chain initialisation options
# ============================================================

def init_low(L):
    """Default: near the GC mode, k ~ 1.28 sqrt(L)."""
    return init_arrays(L)


def init_high(L):
    """All parts equal to 1: k = L (maximally many parts)."""
    from partition_numpy import partition_to_arrays
    import numpy as np
    vals  = np.array([1], dtype=np.int32)
    mults = np.array([L], dtype=np.int32)
    return vals, mults


def init_mid(L, n_warmup=50000, seed=99):
    """
    Start from the mode: run a long chain from the low start,
    use the final state as the 'mid' initialisation.
    This approximates k*_{w_tilde}.
    """
    rng = random.Random(seed)
    vals, mults = init_low(L)
    k = int(mults.sum())
    for _ in range(n_warmup):
        vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
    return vals, mults


# ============================================================
# Check 1 & 4: trace and running mean over the burn-in period
# ============================================================

def trace_and_runmean(L, burnin, seed=42, thin=1):
    """
    Run one chain for `burnin` steps from the low start,
    recording k at every `thin` steps.
    Returns array of k values.
    """
    rng = random.Random(seed)
    vals, mults = init_low(L)
    k = int(mults.sum())
    k_trace = []
    for i in range(burnin):
        vals, mults, k, _ = metropolis_step_fast(vals, mults, k, rng)
        if i % thin == 0:
            k_trace.append(k)
    return np.array(k_trace)


def plot_trace(k_trace, L, burnin, thin, outdir):
    steps = np.arange(len(k_trace)) * thin
    fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=False)

    # Panel 1: full trace
    axes[0].plot(steps, k_trace, lw=0.4, alpha=0.7, color='steelblue')
    axes[0].axhline(k_trace[len(k_trace)//2:].mean(), color='red',
                    lw=1.2, linestyle='--', label=f'2nd-half mean = {k_trace[len(k_trace)//2:].mean():.1f}')
    axes[0].set_xlabel('Step')
    axes[0].set_ylabel('k')
    axes[0].set_title(f'Trace of k during burn-in  (L={L}, {burnin} steps)')
    axes[0].legend(fontsize=9)

    # Panel 2: running mean
    runmean = np.cumsum(k_trace) / (np.arange(len(k_trace)) + 1)
    axes[1].plot(steps, runmean, lw=1.0, color='darkorange')
    axes[1].axhline(k_trace[len(k_trace)//2:].mean(), color='red',
                    lw=1.2, linestyle='--', label='2nd-half mean')
    # mark 1% convergence band
    final = k_trace[len(k_trace)//2:].mean()
    axes[1].axhspan(final * 0.99, final * 1.01, alpha=0.15, color='red',
                    label='±1% band')
    axes[1].set_xlabel('Step')
    axes[1].set_ylabel('Running mean of k')
    axes[1].set_title('Running mean convergence')
    axes[1].legend(fontsize=9)

    plt.tight_layout()
    path = os.path.join(outdir, f'burnin_trace_L{L}.png')
    plt.savefig(path, dpi=120)
    plt.close()
    print(f'  Saved: {path}')

    # Find when running mean first enters ±1% of final value and stays there
    in_band = np.abs(runmean - final) / final < 0.01
    # Find last exit from band (i.e., first stable entry)
    crossings = np.where(~in_band)[0]
    if len(crossings) == 0:
        convergence_step = 0
    else:
        convergence_step = (crossings[-1] + 1) * thin
    return convergence_step, final


# ============================================================
# Check 2: Gelman-Rubin R-hat
# ============================================================

def gelman_rubin(chains):
    """
    Compute Gelman-Rubin R-hat from a list of 1D numpy arrays (chains).
    All chains must have the same length n.
    Returns R-hat scalar.
    """
    m = len(chains)
    n = len(chains[0])
    chain_means = np.array([c.mean() for c in chains])
    chain_vars  = np.array([c.var(ddof=1) for c in chains])

    B = n * chain_means.var(ddof=1)          # between-chain variance * n
    W = chain_vars.mean()                    # within-chain variance

    Vhat = (n - 1) / n * W + (m + 1) / (m * n) * B
    Rhat = math.sqrt(Vhat / W) if W > 0 else float('inf')
    return Rhat, W, B / n


# ============================================================
# Check 3: starting-point sensitivity
# ============================================================

def run_chain(L, n_steps, init_fn, seed, thin=1):
    rng = random.Random(seed)
    vals, mults = init_fn(L)
    k = int(mults.sum())
    k_vals = []
    acc = 0
    for i in range(n_steps):
        vals, mults, k, a = metropolis_step_fast(vals, mults, k, rng)
        acc += a
        if i % thin == 0:
            k_vals.append(k)
    return np.array(k_vals), acc / n_steps


# ============================================================
# Main
# ============================================================

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--L',       type=int, default=500)
    parser.add_argument('--burnin',  type=int, default=None,
                        help='Burn-in length (default: max(20000, 200*L))')
    parser.add_argument('--chains',  type=int, default=4,
                        help='Number of chains for Gelman-Rubin (default 4)')
    parser.add_argument('--outdir',  type=str, default='.')
    args = parser.parse_args()

    L      = args.L
    burnin = args.burnin if args.burnin is not None else max(20000, 200 * L)
    n_chains = args.chains
    thin   = max(1, L // 20)   # same thinning as run_scaling_task
    os.makedirs(args.outdir, exist_ok=True)

    report_lines = []
    report_lines.append(f'Burn-in diagnostics  L={L}  burnin={burnin}  thin={thin}')
    report_lines.append('=' * 60)

    # ----------------------------------------------------------
    # Checks 1 & 4: trace + running mean (single chain, low start)
    # ----------------------------------------------------------
    print(f'\n[1/3] Trace and running mean  (L={L}, burnin={burnin}) ...')
    t0 = time.time()
    k_trace = trace_and_runmean(L, burnin, seed=42, thin=thin)
    conv_step, k_final = plot_trace(k_trace, L, burnin, thin, args.outdir)
    elapsed = time.time() - t0

    report_lines.append(f'\n--- Check 1 & 4: trace and running mean ---')
    report_lines.append(f'  Initial k (step 0):        {k_trace[0]}')
    report_lines.append(f'  2nd-half mean k:           {k_final:.2f}')
    report_lines.append(f'  Running mean enters ±1%'
                        f' band at step:  {conv_step}')
    report_lines.append(f'  Fraction of burnin used:   '
                        f'{conv_step/burnin:.3f}')
    report_lines.append(f'  Elapsed: {elapsed:.1f}s')
    if conv_step < 0.5 * burnin:
        verdict = 'PASS  (running mean converges in first half of burn-in)'
    else:
        verdict = 'WARN  (running mean not converged until second half; consider longer burn-in)'
    report_lines.append(f'  Verdict: {verdict}')

    # ----------------------------------------------------------
    # Check 2: Gelman-Rubin R-hat across n_chains independent chains
    # ----------------------------------------------------------
    print(f'\n[2/3] Gelman-Rubin R-hat  ({n_chains} chains, '
          f'each {2*burnin} steps, using 2nd half) ...')
    t0 = time.time()
    chains_second_half = []
    for i in range(n_chains):
        k_chain, ar = run_chain(L, 2 * burnin, init_low, seed=42 + i, thin=thin)
        chains_second_half.append(k_chain[len(k_chain)//2:])
        print(f'  Chain {i+1}/{n_chains}: mean k = {k_chain[len(k_chain)//2:].mean():.2f}, '
              f'acc = {ar:.3f}')
    elapsed = time.time() - t0

    Rhat, W, B = gelman_rubin(chains_second_half)
    report_lines.append(f'\n--- Check 2: Gelman-Rubin R-hat ({n_chains} chains) ---')
    for i, c in enumerate(chains_second_half):
        report_lines.append(f'  Chain {i+1}: mean k = {c.mean():.2f}  std = {c.std():.2f}')
    report_lines.append(f'  Within-chain variance W:   {W:.2f}')
    report_lines.append(f'  Between-chain variance B:  {B:.2f}')
    report_lines.append(f'  R-hat:                     {Rhat:.5f}')
    report_lines.append(f'  Elapsed: {elapsed:.1f}s')
    if Rhat < 1.01:
        verdict = 'PASS  (R-hat < 1.01)'
    elif Rhat < 1.05:
        verdict = 'WARN  (1.01 <= R-hat < 1.05; marginal convergence)'
    else:
        verdict = 'FAIL  (R-hat >= 1.05; chains have not mixed)'
    report_lines.append(f'  Verdict: {verdict}')

    # ----------------------------------------------------------
    # Check 3: starting-point sensitivity
    # ----------------------------------------------------------
    print(f'\n[3/3] Starting-point sensitivity  (low / mid / high) ...')
    t0 = time.time()

    print(f'  Running low-start chain ...')
    k_low,  ar_low  = run_chain(L, burnin, init_low,  seed=42, thin=thin)
    print(f'  Running high-start chain ...')
    k_high, ar_high = run_chain(L, burnin, init_high, seed=42, thin=thin)
    print(f'  Running mid-start chain (requires 50k warmup) ...')
    def init_mid_L(L_):
        return init_mid(L_, n_warmup=50000, seed=99)
    k_mid,  ar_mid  = run_chain(L, burnin, init_mid_L, seed=42, thin=thin)
    elapsed = time.time() - t0

    n2 = len(k_low) // 2
    means = {
        'low  (k0 ~ 1.28 sqrt(L))': k_low[n2:].mean(),
        'high (k = L, all-ones)   ': k_high[n2:].mean(),
        'mid  (approx k*_wtilde)  ': k_mid[n2:].mean(),
    }
    report_lines.append(f'\n--- Check 3: starting-point sensitivity ---')
    report_lines.append(f'  2nd-half means after {burnin} steps:')
    for label, m in means.items():
        report_lines.append(f'    {label}: {m:.2f}')
    vals_list = list(means.values())
    spread = max(vals_list) - min(vals_list)
    rel_spread = spread / np.mean(vals_list)
    report_lines.append(f'  Max spread:    {spread:.2f}  ({rel_spread*100:.2f}%)')
    report_lines.append(f'  Elapsed: {elapsed:.1f}s')
    if rel_spread < 0.02:
        verdict = 'PASS  (all starts agree within 2%)'
    elif rel_spread < 0.05:
        verdict = 'WARN  (spread 2–5%; consider longer burn-in for low start)'
    else:
        verdict = 'FAIL  (>5% spread; low-start chain not converged)'
    report_lines.append(f'  Verdict: {verdict}')

    # Plot starting-point comparison
    fig, ax = plt.subplots(figsize=(10, 4))
    steps = np.arange(len(k_low)) * thin
    ax.plot(steps, k_low,  lw=0.5, alpha=0.7, label='low start',  color='steelblue')
    ax.plot(steps, k_high, lw=0.5, alpha=0.7, label='high start', color='tomato')
    ax.plot(steps, k_mid,  lw=0.5, alpha=0.7, label='mid start',  color='seagreen')
    ax.axvline(burnin // 2, color='k', lw=1, linestyle=':', label='half-burnin')
    ax.set_xlabel('Step')
    ax.set_ylabel('k')
    ax.set_title(f'Starting-point sensitivity  (L={L})')
    ax.legend(fontsize=9)
    plt.tight_layout()
    path2 = os.path.join(args.outdir, f'burnin_starts_L{L}.png')
    plt.savefig(path2, dpi=120)
    plt.close()
    print(f'  Saved: {path2}')

    # ----------------------------------------------------------
    # Write report
    # ----------------------------------------------------------
    report_lines.append('\n' + '=' * 60)
    # Overall verdict
    all_pass = ('PASS' in verdict)
    # (we check the last verdict; a fuller check would track all three)
    report_lines.append('Overall: re-read individual verdicts above.')

    report_text = '\n'.join(report_lines)
    print('\n' + report_text)

    out_txt = os.path.join(args.outdir, f'burnin_convergence_L{L}.txt')
    with open(out_txt, 'w') as f:
        f.write(report_text + '\n')
    print(f'\n  Report saved: {out_txt}')


if __name__ == '__main__':
    main()
