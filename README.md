# Reptation Partition Scripts

This repository contains five standalone Python scripts used to analyze the connection between reptating polymer run-length encodings and integer partition ensembles.

## Requirements

- Python 3.8+
- Standard library only (`math`, `random`, `collections`, `itertools`, `functools`, `time`, etc.)
- Optional: `reptation_core.py` in the same directory (used for exact partition enumeration in some validation paths)

No third-party packages are required.

## Repository Layout

- `scripts/reptation_mc.py`
- `scripts/reptation_mc_uniform.py`
- `scripts/reptation_partition_simulations.py`
- `scripts/reptation_be_saddle.py`
- `scripts/reptation_rld.py`

## Script Overview

### 1) `reptation_mc.py`
Metropolis-Hastings sampler over partitions of fixed `L` under weighted measure

- Weight: `w_tilde(lambda) = (k! / prod_j m_j!) * 2^k`
- Moves: transfer, split, merge
- Includes:
  - detailed-balance-correct proposal counting
  - validation against exact weighted distribution (when `reptation_core.py` is available)
  - scaling diagnostics and run-length distribution reports

Run:

```bash
python scripts/reptation_mc.py
```

### 2) `reptation_mc_uniform.py`
Fast Metropolis-Hastings sampler under the **uniform** measure over partitions of `L`

- Target: each partition has equal probability
- Acceptance uses forward/reverse move-count ratios
- Includes:
  - exact-vs-MC validation (optional `reptation_core.py`)
  - mode convergence study for `k`
  - speed benchmark
  - run-length diagnostics

Run:

```bash
python scripts/reptation_mc_uniform.py
```

### 3) `reptation_partition_simulations.py`
Exact enumeration and verification suite for core combinatorial claims

Main sections:

- Verification of decoupling theorem: `w_strong = w_H * w_V`
- Analysis of the `p=2` formula
- Counterexample-driven correction to distinct-parts heuristic `w_H = 2^p`
- Exact degeneracy `w(lambda)` analysis and most-probable partition study up to `L=20`

Run:

```bash
python scripts/reptation_partition_simulations.py
```

Notes:

- Full run includes a heavier final step (around ~60s depending on machine)
- Early verification sections complete much faster

### 4) `reptation_be_saddle.py`
Saddle-point analysis comparing grand-canonical vs canonical behavior

- Computes exact `p_k(L)` table via dynamic programming
- Compares:
  - grand-canonical prediction (`~ sqrt(L)` scale)
  - canonical statistics from exact counts (including `sqrt(L) log(sqrt(L))` behavior)
- Reports entropy comparisons and asymptotic diagnostics

Run:

```bash
python scripts/reptation_be_saddle.py
```

### 5) `reptation_rld.py`
Run-length distribution analysis for uniform partition ensemble

- Compares MC marginals to:
  - Bose-Einstein (BE) asymptotic form
  - geometric approximation
  - exact finite-`L` marginal formula (for small `L`)
- Includes KL-divergence comparisons and convergence-by-`L` study

Run:

```bash
python scripts/reptation_rld.py
```

## Typical Usage

From repo root:

```bash
python scripts/reptation_partition_simulations.py
python scripts/reptation_mc.py
python scripts/reptation_mc_uniform.py
python scripts/reptation_be_saddle.py
python scripts/reptation_rld.py
```

If your environment uses `python3` instead of `python`, replace accordingly.

## Notes on Inter-Script Dependencies

- `reptation_mc_uniform.py` imports utilities from `reptation_mc.py`
- `reptation_rld.py` imports from both `reptation_mc.py` and `reptation_mc_uniform.py`
- Keep all five scripts together in the same `scripts/` directory

## Reproducibility

Most simulations set explicit RNG seeds in-script (commonly `seed=42` or `7`).
For alternative experiments, edit the `__main__` blocks or function call arguments.

## License

Add your preferred license (for example, MIT) before publishing the repository.
