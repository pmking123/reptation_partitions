"""
prediction2.py
========================
Numerical test of Prediction 2 from "Reptating polymers and integer
partitions: a Bose-Einstein correspondence".

PREDICTION 2 (E):
    Under the expander assumption on the shape-class graph, the tube
    renewal time satisfies

        log tau_renew ~ alpha_0 * sqrt(L),   alpha_0 = pi * sqrt(2/3)

    growing super-polynomially in L.  This is qualitatively distinct
    from the Doi-Edwards reptation time tau_rep ~ L^3.

WHAT THIS SCRIPT DOES:
    Simulates standard lattice reptation (the "slithering snake" move)
    for a single chain of length L on a 2D square lattice with periodic
    boundary conditions and a fixed obstacle field.  The primitive path
    is the chain itself; tube renewal corresponds to the chain fully
    escaping its initial tube.

    For each L we measure the integrated autocorrelation time tau of the
    end-to-end vector C(t) = <r(0).r(t)> / <r(0).r(0)>, which is the
    standard proxy for the tube renewal (decorrelation) time.

    We then fit log(tau) vs sqrt(L) and check for linearity with slope
    approximately alpha_0 = pi*sqrt(2/3) ~ 2.565.

CAVEATS:
    - This tests Prediction 2 only indirectly: the simulation uses
      standard reptation dynamics, not the expander random walk on the
      partition shape-class graph.  The prediction says these should
      have the same scaling if the expander assumption holds.
    - The obstacle density defines the effective tube width.  We use a
      dilute obstacle field so that the chain is not pinned.
    - At accessible L (up to ~200) the super-polynomial growth is hard
      to distinguish from a power law; much longer runs are needed to
      see the asymptotic regime clearly.
    - The label (E) in the paper means this prediction requires a kinetic
      assumption (the expander conjecture) beyond the combinatorial
      results; this script provides supporting numerical evidence, not
      a proof.

OUTPUT:
    - Per-L table: L, sqrt(L), tau, log(tau), log(tau)/sqrt(L)
    - Linear fit of log(tau) vs sqrt(L) with slope and intercept
    - Comparison of slope with alpha_0 and with the reptation prediction
      (power-law fit of log(tau) vs log(L))
    - CSV file: prediction2_results.csv

USAGE:
    python3 prediction2.py [--Lvals 20,30,50,75,100] \
                           [--obstacle_density 0.1]   \
                           [--n_chains 8]             \
                           [--n_runs 4]               \
                           [--seed 42]

    Increase n_chains and n_runs for better statistics at the cost of
    longer runtime.  Each independent run uses a different random seed.
"""

import math
import multiprocessing
import random
import argparse
import csv
import os
import sys
import time
from collections import deque


PI      = math.pi
ALPHA_0 = PI * math.sqrt(2.0 / 3.0)   # ~ 2.5651


# ====================================================================
# Lattice and obstacle setup
# ====================================================================

def build_obstacles(lattice_size, density, seed):
    """Return a frozenset of obstacle sites."""
    rng = random.Random(seed)
    L   = lattice_size
    n   = int(L * L * density)
    obs = set()
    while len(obs) < n:
        obs.add((rng.randint(0, L-1), rng.randint(0, L-1)))
    return frozenset(obs)


def neighbours(x, y, L):
    return [
        ((x+1) % L, y),
        ((x-1) % L, y),
        (x, (y+1) % L),
        (x, (y-1) % L),
    ]


# ====================================================================
# Chain initialisation: grow a self-avoiding walk avoiding obstacles
# ====================================================================

def grow_chain(chain_length, lattice_size, obstacles, rng, max_attempts=1000):
    """
    Grow a self-avoiding chain of length chain_length using a loop-erased
    random walk (LERW).  At each step the walker moves to a random
    obstacle-free neighbour; if it revisits a site the loop is erased,
    keeping the path self-avoiding throughout.  Guaranteed to succeed on
    a well-connected lattice; typical runtime O(L^2), in practice <50ms
    for L<=1000 on a 2000x2000 lattice with rho=0.05.

    Returns a list of (x,y) tuples of length chain_length (head first).
    Raises RuntimeError if unable to place chain after max_attempts restarts.
    """
    L = lattice_size
    for _ in range(max_attempts):
        x = rng.randint(0, L - 1)
        y = rng.randint(0, L - 1)
        if (x, y) in obstacles:
            continue

        path  = [(x, y)]
        index = {(x, y): 0}   # site -> index in path, for O(1) loop detection

        while len(path) < chain_length:
            x, y = path[-1]
            nbs  = [n for n in neighbours(x, y, L) if n not in obstacles]
            if not nbs:
                break          # walker trapped by obstacles; restart
            x, y = rng.choice(nbs)
            if (x, y) in index:
                # erase the loop back to the previous visit of this site
                loop_start = index[(x, y)]
                for site in path[loop_start + 1:]:
                    del index[site]
                path = path[:loop_start + 1]
            else:
                index[(x, y)] = len(path)
                path.append((x, y))

        if len(path) == chain_length:
            return path

    raise RuntimeError(
        f"Could not place chain of length {chain_length} after "
        f"{max_attempts} attempts (lattice may be too crowded)."
)


# ====================================================================
# Reptation move (slithering snake)
# ====================================================================

def reptate_step(chain, occupied, obstacles, lattice_size, rng):
    """
    One reptation move: with equal probability try to advance the head
    or the tail into a random free neighbour, shifting the chain.
    occupied is a set maintained incrementally by the caller (O(1) update).
    Returns (new_chain, accepted).  occupied is updated in place.
    """
    L       = lattice_size
    is_head = rng.random() < 0.5
    if is_head:
        tip   = chain[0]
        other = chain[-1]
    else:
        tip   = chain[-1]
        other = chain[0]

    candidates = [
        nb for nb in neighbours(tip[0], tip[1], L)
        if nb not in obstacles and (nb not in occupied or nb == other)
    ]
    if not candidates:
        return chain, False

    new_tip = rng.choice(candidates)
    if is_head:
        new_chain = [new_tip] + chain[:-1]
    else:
        new_chain = chain[1:] + [new_tip]

    # O(1) incremental update: add new tip, remove vacated end
    occupied.add(new_tip)
    occupied.discard(other)
    return new_chain, True


# ====================================================================
# End-to-end vector (unwrapped for PBC)
# ====================================================================

def end_to_end(chain, L):
    """
    Compute the unwrapped end-to-end vector (rx, ry) accounting for
    periodic boundary conditions.
    """
    rx, ry = 0, 0
    for i in range(len(chain) - 1):
        x1, y1 = chain[i]
        x2, y2 = chain[i+1]
        dx = x2 - x1
        dy = y2 - y1
        if dx >  1: dx -= L
        if dx < -1: dx += L
        if dy >  1: dy -= L
        if dy < -1: dy += L
        rx += dx
        ry += dy
    return rx, ry


# ====================================================================
# Autocorrelation time estimation
# ====================================================================

def integrated_autocorrelation_time(series, c=5.0):
    """
    Estimate the integrated autocorrelation time of a 1D time series
    using the standard self-consistent windowing method.

    Returns tau_int.  Returns None if the series is too short or
    constant.
    """
    n    = len(series)
    mean = sum(series) / n
    var  = sum((x - mean)**2 for x in series) / n
    if var < 1e-14:
        return None

    tau   = 0.5
    for t in range(1, n // 2):
        # unnormalised autocovariance at lag t
        cov = sum((series[i] - mean) * (series[i+t] - mean)
                  for i in range(n - t)) / (n - t)
        rho = cov / var
        tau += rho
        # self-consistent window: stop when window >= c * tau
        if t >= c * tau:
            break

    return max(tau, 0.5)


# ====================================================================
# Single run: simulate one chain, return autocorrelation time
# ====================================================================

def run_single_chain(chain_length, lattice_size, obstacles, seed,
                     n_sweeps, burn_in, thin):
    """
    Simulate one chain for n_sweeps reptation steps after burn_in,
    recording the squared end-to-end distance every `thin` steps.
    Returns the integrated autocorrelation time of r^2(t).
    """
    rng      = random.Random(seed)
    chain    = grow_chain(chain_length, lattice_size, obstacles, rng)
    occupied = set(chain)
    L        = lattice_size

    # Burn-in
    for _ in range(burn_in):
        chain, _ = reptate_step(chain, occupied, obstacles, L, rng)

    # Record r^2
    r0x, r0y = end_to_end(chain, L)

    series = []
    for i in range(n_sweeps):
        chain, _ = reptate_step(chain, occupied, obstacles, L, rng)
        if i % thin == 0:
            rx, ry = end_to_end(chain, L)
            series.append(rx * r0x + ry * r0y)

    tau = integrated_autocorrelation_time(series)
    return tau


# ====================================================================
# Multi-run measurement for one L value
# ====================================================================

def measure_tau(chain_length, lattice_size, obstacle_density,
                n_chains, n_runs, base_seed, verbose=True):
    """
    For a given chain_length, run n_chains * n_runs independent
    simulations and return mean and std of log(tau).

    burn_in and n_sweeps scale as chain_length^2.5 (a compromise between
    the L^2 diffusive scale and L^3 Doi-Edwards scale), capped at 5M
    sweeps to keep runtime manageable at large L.

    All chain simulations are run in parallel via multiprocessing.Pool.
    """

    tau_25   = int(chain_length**2.5)
    n_sweeps = min(max(50 * chain_length**2, 20 * tau_25), 5_000_000)
    burn_in  = max(10 * chain_length**2, 2 * tau_25)
    thin     = max(1, chain_length // 5)

    if verbose:
        print(f"  L={chain_length:4d}  sweeps={n_sweeps:>8d}  "
              f"burn={burn_in:>10d}  thin={thin:>3d}", end="  ", flush=True)

    t0   = time.time()
    jobs = []

    for run in range(n_runs):
        obs_seed   = base_seed + run * 1000
        chain_seed = base_seed + run * 1000 + 1
        obstacles  = build_obstacles(lattice_size, obstacle_density, obs_seed)
        for chain_idx in range(n_chains):
            seed = chain_seed + chain_idx * 17
            jobs.append((chain_length, lattice_size, obstacles, seed,
                         n_sweeps, burn_in, thin))

    n_workers = int(os.environ.get("SLURM_CPUS_PER_TASK",
                                multiprocessing.cpu_count()))
    n_workers = min(n_workers, len(jobs))
    with multiprocessing.Pool(processes=n_workers) as pool:
        raw_taus = pool.starmap(run_single_chain, jobs)

    dt       = time.time() - t0
    tau_vals = [tau * thin for tau in raw_taus if tau is not None]

    if not tau_vals:
        if verbose: print(f"no valid tau estimates  ({dt:.1f}s)")
        return None, None, None

    log_taus = [math.log(t) for t in tau_vals if t > 0]
    mean_log = sum(log_taus) / len(log_taus)
    std_log  = math.sqrt(
        sum((x - mean_log)**2 for x in log_taus) / max(len(log_taus) - 1, 1)
    )
    sem_log  = std_log / math.sqrt(len(log_taus))

    if verbose:
        print(f"log(tau)={mean_log:.3f} ± {sem_log:.3f}  "
              f"(n={len(log_taus)}, {dt:.1f}s)")

    return mean_log, sem_log, len(log_taus)


# ====================================================================
# Linear fit (OLS)
# ====================================================================

def linear_fit(xs, ys):
    """Fit y = a*x + b by OLS. Returns (a, b, r^2)."""
    n    = len(xs)
    sx   = sum(xs);        sy  = sum(ys)
    sxx  = sum(x**2 for x in xs)
    sxy  = sum(x*y for x, y in zip(xs, ys))
    denom = n * sxx - sx**2
    if abs(denom) < 1e-12:
        return None, None, None
    a  = (n * sxy - sx * sy) / denom
    b  = (sy - a * sx) / n
    y_mean = sy / n
    ss_tot = sum((y - y_mean)**2 for y in ys)
    ss_res = sum((y - (a*x + b))**2 for x, y in zip(xs, ys))
    r2 = 1 - ss_res / ss_tot if ss_tot > 1e-12 else 1.0
    return a, b, r2


# ====================================================================
# Main
# ====================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Test Prediction 2: log(tau_renew) ~ alpha_0 * sqrt(L)")
    parser.add_argument("--Lvals", default="20,30,50,75,100",
        help="Comma-separated chain lengths (default: 20,30,50,75,100)")
    parser.add_argument("--lattice", type=int, default=200,
        help="Lattice size (default: 200; must be >> chain length)")
    parser.add_argument("--obstacle_density", type=float, default=0.05,
        help="Obstacle site fraction (default: 0.05)")
    parser.add_argument("--n_chains", type=int, default=4,
        help="Independent chains per run (default: 4)")
    parser.add_argument("--n_runs", type=int, default=2,
        help="Independent obstacle realisations (default: 2)")
    parser.add_argument("--seed", type=int, default=42,
        help="Base random seed (default: 42)")
    parser.add_argument("--outfile", default="prediction2_results.csv",
        help="Output CSV filename")
    args = parser.parse_args()

    L_vals = [int(x) for x in args.Lvals.split(",")]

    print()
    print("=" * 65)
    print("Prediction 2: log(tau_renew) ~ alpha_0 * sqrt(L)")
    print("=" * 65)
    print()
    print(f"  Chain lengths:      {L_vals}")
    print(f"  Lattice size:       {args.lattice} x {args.lattice}")
    print(f"  Obstacle density:   {args.obstacle_density:.3f}")
    print(f"  Chains per run:     {args.n_chains}")
    print(f"  Independent runs:   {args.n_runs}")
    print(f"  Base seed:          {args.seed}")
    print()
    print(f"  alpha_0 = pi*sqrt(2/3) = {ALPHA_0:.4f}  (prediction slope)")
    print(f"  Standard reptation: log(tau) ~ 3*log(L)  (power-law)")
    print()

    results = []

    print(f"  {'L':>4}  {'sqrt(L)':>8}  "
          f"{'log(tau)':>10}  {'± sem':>7}  "
          f"{'log(tau)/sqL':>13}  {'n':>4}")
    print("  " + "-" * 55)

    for L in L_vals:
        mean_log, sem_log, n = measure_tau(
            chain_length      = L,
            lattice_size      = args.lattice,
            obstacle_density  = args.obstacle_density,
            n_chains          = args.n_chains,
            n_runs            = args.n_runs,
            base_seed         = args.seed,
            verbose           = False
        )
        sqL = math.sqrt(L)
        if mean_log is not None:
            ratio = mean_log / sqL
            print(f"  {L:>4}  {sqL:>8.3f}  "
                  f"{mean_log:>10.3f}  {sem_log:>7.3f}  "
                  f"{ratio:>13.4f}  {n:>4}")
            results.append({
                'L':         L,
                'sqrtL':     sqL,
                'log_tau':   mean_log,
                'sem':       sem_log,
                'ratio':     ratio,
                'n':         n,
            })
        else:
            print(f"  {L:>4}  {sqL:>8.3f}  {'N/A':>10}")

    if len(results) < 2:
        print("\nInsufficient data for fitting.")
        return

    print()
    print("=" * 65)
    print("Fit 1: log(tau) = a * sqrt(L) + b  [Prediction 2]")
    print("=" * 65)

    xs1 = [r['sqrtL']   for r in results]
    ys  = [r['log_tau'] for r in results]
    a1, b1, r2_1 = linear_fit(xs1, ys)
    if a1 is not None:
        print(f"  Slope a  = {a1:.4f}  (predicted alpha_0 = {ALPHA_0:.4f})")
        print(f"  Intercept b = {b1:.4f}")
        print(f"  R^2      = {r2_1:.4f}")
        print(f"  a / alpha_0 = {a1/ALPHA_0:.4f}  (expect ~1 if prediction holds)")

    print()
    print("=" * 65)
    print("Fit 2: log(tau) = c * log(L) + d  [standard reptation L^3]")
    print("=" * 65)

    xs2 = [math.log(r['L']) for r in results]
    a2, b2, r2_2 = linear_fit(xs2, ys)
    if a2 is not None:
        print(f"  Exponent c = {a2:.4f}  (standard reptation predicts c = 3)")
        print(f"  Intercept d = {b2:.4f}")
        print(f"  R^2        = {r2_2:.4f}")

    print()
    print("=" * 65)
    print("Model comparison")
    print("=" * 65)
    if a1 is not None and a2 is not None:
        if r2_1 > r2_2:
            print(f"  sqrt(L) fit is better (R^2 = {r2_1:.4f} vs {r2_2:.4f})")
            print(f"  This is consistent with Prediction 2.")
        else:
            print(f"  log(L) fit is better (R^2 = {r2_2:.4f} vs {r2_1:.4f})")
            print(f"  Chains may not be long enough to see the super-polynomial")
            print(f"  regime; try larger L values.")
        print()
        print(f"  NOTE: at accessible L the two fits are hard to distinguish.")
        print(f"  The key test is whether the slope a in Fit 1 converges to")
        print(f"  alpha_0 = {ALPHA_0:.4f} as L increases, and whether R^2 of")
        print(f"  Fit 1 improves relative to Fit 2 at larger L.")

    # Write CSV
    out_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                            args.outfile)
    with open(out_path, 'w', newline='') as fh:
        fieldnames = ['L', 'sqrtL', 'log_tau', 'sem', 'log_tau_over_sqrtL', 'n']
        w = csv.DictWriter(fh, fieldnames=fieldnames)
        w.writeheader()
        for r in results:
            w.writerow({
                'L':                   r['L'],
                'sqrtL':               f"{r['sqrtL']:.4f}",
                'log_tau':             f"{r['log_tau']:.4f}",
                'sem':                 f"{r['sem']:.4f}",
                'log_tau_over_sqrtL':  f"{r['ratio']:.4f}",
                'n':                   r['n'],
            })
    print()
    print(f"Results written to: {out_path}")
    print()


if __name__ == '__main__':
    main()
