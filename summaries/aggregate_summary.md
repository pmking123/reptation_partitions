# SLURM Output Aggregator

## Summary of `aggregate.py`

---

## (a) Theory

`aggregate.py` contains no physics logic. It is a post-processing utility that
collects JSON result files written by the SLURM job-array tasks of `mc_fast.py`
and `mc_uniform_fast.py` and reproduces the same formatted tables that the
standalone (sequential) versions of those scripts print. Its purpose is to
allow individual $L$-value tasks — which run independently in parallel across
cluster nodes — to be combined into a single coherent output after the job
array completes.

---

## (b) What the script does

The script operates in three modes, selected via `--mode`:

- **`weighted`**: reads `scaling_{L}.json` files produced by `mc_fast.py` SLURM
  tasks and prints the weighted-measure scaling study table
  ($\langle k\rangle/\sqrt{L}$, $\langle k\rangle/(\sqrt{L}\log L)$,
  $\sigma_k/L^{1/4}$, $\tau_\mathrm{int}$, acceptance rate, elapsed time).

- **`uniform`**: reads `mode_{L}.json` and `rld_{L}.json` files produced by
  `mc_uniform_fast.py` SLURM tasks and prints (a) the mode/mean convergence
  table for the uniform measure and (b) the per-$L$ run-length distribution
  tables with BE and geometric comparisons and KL divergences. With
  `--paper`, also prints a side-by-side comparison of MC KL values against
  hardcoded paper reference values (see §(e) below).

- **`both`** (default): runs all three of the above, plus a cross-check table
  comparing weighted-measure and uniform-measure $\langle k\rangle$ at
  shared $L$ values.

All output is printed to stdout. No CSV or other files are written.

---

## (c) How the script does it

### JSON result format

Each JSON file is a result dict written by a single SLURM task. The script
reads each file and extracts the relevant fields:

- `scaling_{L}.json` (from `mc_fast.py`): keys `L`, `km`, `km_over_sqL`,
  `ratio_logL`, `ksig_over_L14`, `tau`, `acc`, `elapsed`.
- `mode_{L}.json` (from `mc_uniform_fast.py`): keys `L`, `modal_k`,
  `modal_over_sqL`, `ratio_logL`, `km`, `km_over_sqL`, `ksig_over_L14`,
  `tau`, `elapsed`.
- `rld_{L}.json` (from `mc_uniform_fast.py`): keys `L`, `k_mean`, `r_bar`,
  `alpha`, `pred_k_gc`, `pred_r_gc`, `kl_be`, `kl_geom`, `kl_ratio`,
  `f_mc`, `f_be`, `f_geom` (the last three as string-keyed dicts of
  $j \to f(j)$ values).

### L-value lists

The $L$ values to aggregate are imported directly from the producing scripts:
`SCALING_L_VALS` from `mc_fast.py` and `MODE_L_VALS`, `RLD_L_VALS` from
`mc_uniform_fast.py`. This ensures the aggregator always uses exactly the
same $L$ list as the job array, with no hardcoded duplication. Missing
result files are flagged with a warning rather than causing an error.

### Directory defaults

Default result directories are `results_mc` (weighted) and
`results_mc_uniform` (uniform), matching the `--outdir` defaults of the
producing scripts. Both can be overridden via `--indir` (for single-mode
use) or `--weighted-indir` / `--uniform-indir` (for `--mode both`).

### Cross-check

`cross_check` identifies $L$ values present in both `SCALING_L_VALS` and
`MODE_L_VALS`, reads the corresponding JSON files from both result
directories, and prints the ratio $\langle k\rangle_\mathrm{weighted} /
\langle k\rangle_\mathrm{uniform}$. Since the weighted measure
$\tilde{w}(\lambda) \propto k!\,2^k/\prod_j m_j!$ strongly favours
high-$k$ partitions, this ratio should substantially exceed 1 at all $L$,
confirming the distinction between the two ensembles discussed in §6.3
of the paper.

---

## (d) Output produced

The script produces no files. All output is printed to stdout. The table
formats mirror those of the producing scripts exactly, allowing the SLURM
and sequential outputs to be compared line-by-line.

---

## (e) One issue: stale `PAPER_KL` constant

The script contains a hardcoded constant:

```python
PAPER_KL = {
    100: (0.040, 0.072),
    200: (0.023, 0.099),
    500: (0.026, 0.145),
}
```

These values are from a version of `mc_uniform_fast.py` that used the GC mean
as the BE normaliser. The correct values in the current paper (§8.3, from
`rld_fast.py` with canonical normaliser) are substantially different:

| $L$ | $D_\mathrm{KL}(\cdot\|f_\mathrm{BE})$ (current) | $D_\mathrm{KL}(\cdot\|f_\mathrm{geom})$ (current) |
|---|---|---|
| 100 | $\approx 0.0018$ | $\approx 0.080$ |
| 200 | $\approx 0.0013$ | $\approx 0.103$ |
| 500 | $\approx 0.0010$ | $\approx 0.147$ |

**The `--paper` flag should not be used** until `PAPER_KL` is updated to
reflect the current paper values. In the meantime, the plain `--mode uniform`
output (without `--paper`) is correct and unaffected by this constant.
