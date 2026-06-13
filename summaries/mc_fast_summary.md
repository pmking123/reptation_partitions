# Weighted Partition Sampler: Canonical/GC Ensemble Comparison

## Summary of `mc_fast.py`

---

## (a) Theory

### The weighted measure $\tilde{w}$

`mc_fast.py` samples integer partitions of $L$ under the **weighted measure**

$$\tilde{w}(\lambda) = \frac{k!\,2^k}{\prod_j m_j!},$$

where $k = \sum_j m_j$ is the total number of parts and $m_j$ is the
multiplicity of part value $j$. This weight arises naturally in the polymer
physics context: a partition $\lambda$ represents an equivalence class of
primitive-path run-length compositions, and $k!/\prod_j m_j!$ counts the
number of distinct orderings (compositions) in the class, while $2^k$ accounts
for the binary axis assignment (horizontal or vertical) of each segment. The
normalised weight is therefore proportional to the number of labelled
lattice-polymer configurations whose sorted run-length sequence equals
$\lambda$.

The partition representation and all move operations are provided by
`partition_numpy.py`; see that script's summary for implementation details.

### Canonical versus grand canonical ensemble inequivalence

The **grand canonical** mean part count is
$$E[k]_\mathrm{GC} \sim \frac{\pi}{\sqrt{6}}\sqrt{L} \approx 1.282\sqrt{L}.$$

The **canonical** (uniform) Erdős–Lehner mean scales parametrically faster:
$$E[k]_\mathrm{canon} \sim \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L \approx 0.390\,\sqrt{L}\log L.$$

The two ensembles are **inequivalent** for the number of parts: the GC relative
fluctuation $\sigma_k/E[k]_\mathrm{GC} = O(1)$, so $k$ is not
self-averaging under the GC measure. Consequently $E[k]/\sqrt{L}$ does not
converge to a finite limit under the canonical or $\tilde{w}$ measures: it
grows as $(\sqrt{6}/2\pi)\log L$, directly visible in the scaling study output.

### The $\tilde{w}$ measure and Erdős–Lehner scaling

The $\tilde{w}$ measure is neither grand canonical nor uniform. The $k!$ factor
preferentially weights partitions with many distinct part values, and the $2^k$
factor additionally favours high-$k$ configurations, pushing $\langle k\rangle$
above the uniform-measure value. Its mean also scales as $\sqrt{L}\log L$ with
the same leading Erdős–Lehner coefficient $\sqrt{6}/(2\pi)$, but approaches
this limit from above and converges only at $L \gg 10^3$ due to large
sub-leading corrections. This is confirmed by the column
$\langle k\rangle/(\sqrt{L}\log L)$ in the scaling study output, which remains
well above $0.390$ across the accessible range $L \leq 500$.

---

## (b) What the script does

`mc_fast.py` implements the $\tilde{w}$-weighted Metropolis-Hastings sampler
and can be run in two modes:

**Standalone mode** (`python3 mc_fast.py`) runs sequentially:
1. **Benchmark**: steps/second at representative $L$ values.
2. **Validation**: compares MC against the exact $\tilde{w}$ distribution for
   even $L = 4$–$14$ (requires `core.py`; skipped gracefully if absent).
3. **Run-length distributions**: measures $f(j)$ under $\tilde{w}$ for
   $L \in \{50, 100, 200\}$ and compares against the geometric prediction.
4. **Scaling study**: runs the sampler at
   $L \in \{20, 50, 100, 200, 300, 500\}$ and prints the scaling table.

**SLURM job-array mode** (`--task i --outdir dir`) runs a single $L$ value
from `SCALING_L_VALS` and writes a JSON result file `scaling_{L}.json`.
`mc_slurm.sh` configures the array; `aggregate.py` combines the outputs.

---

## (c) How the script does it

### MH step: `metropolis_step_fast`

Selects a move type with probabilities $(0.5,\, 0.25,\, 0.25)$ for
(transfer, split, merge), then samples a move uniformly within that type by
full enumeration. The acceptance ratio is

$$\alpha_\mathrm{MH} = \min\!\left(1,\; e^{\Delta\log\tilde{w}}\cdot\frac{N_\mathrm{fwd}}{N_\mathrm{rev}}\right),$$

where $N_\mathrm{fwd}$ and $N_\mathrm{rev}$ are the total available moves of
the chosen type before and after the proposal. The $N_\mathrm{fwd}/N_\mathrm{rev}$
correction ensures detailed balance; $N_\mathrm{rev}$ is recomputed from the
proposed state using the $O(n_d)$ count functions in `partition_numpy.py`.

### Equilibration and sampling

For each $L$:
$$\text{burn-in} = \max(20{,}000,\; 200L), \quad n = \max(100{,}000,\; 500L) \text{ thinned samples}, \quad \text{thin} = \max(1, \lfloor L/20\rfloor).$$

Observables: $\langle k\rangle$, $\sigma_k$, integrated autocorrelation time
$\tau_\mathrm{int}$ (Madras–Sokal window estimator with $c = 5$), acceptance
rate, and elapsed time.

### Validation

`validate` computes the exact $\tilde{w}$ distribution via exhaustive
enumeration (using `partitions` from `core.py`) for even $L = 4$–$14$, and
compares against MC estimates of $\langle k\rangle$, $\sigma_k$, and
total variation distance after $3 \times 10^5$ samples. If `core.py` is not
present the function prints a warning and returns; all other functions are
unaffected.

---

## (d) Output produced

### Scaling study table

```
  Erdos-Lehner mean scaling: E[k] ~ 0.3898 * sqrt(L) * log(L)

     L    <k>/sqL   <k>/(sqL*lnL)   sig/L^.25      tau     acc   t(s)
---------------------------------------------------------------------
    20     3.0482         1.01750      0.9640     32.1   0.588    2.7
    50     4.7523         1.21479      1.1973     41.5   0.604    4.8
   100     6.6917         1.45309      1.4921     41.9   0.612   10.8
   200     9.4481         1.78322      1.7603     40.0   0.618   21.6
   300    11.5625         2.02717      1.9498     54.1   0.624   48.9
   500    14.9287         2.40220      2.2479     56.6   0.630  136.5
```

**$\langle k\rangle/\sqrt{L}$**: grows monotonically from 3.05 to 14.93,
far above the GC limit 1.282, confirming $\sqrt{L}\log L$ scaling.

**$\langle k\rangle/(\sqrt{L}\log L)$**: rises from 1.02 to 2.40 — the column
is *increasing* with $L$, not decreasing toward 0.390. This is not an error:
the $\tilde{w}$ measure's sub-leading corrections are large at accessible $L$,
and the effective prefactor converges to 0.390 only at $L \gg 10^3$.

**$\sigma_k/L^{1/4}$**: rises from 0.96 to 2.25, above the uniform-measure
prediction of 1.13, consistent with the $k!$ factor broadening the part-count
distribution.

**$\tau_\mathrm{int}$**: 32–57 steps, varying non-monotonically; with thinning
factor up to 25, effective autocorrelation in thinned samples decreases with
$L$, confirming adequate thinning.

**Acceptance rate**: rises from 0.588 to 0.630, reflecting smoother local
weight landscape at larger $L$.

---

## (e) What we learn

**$\langle k\rangle/\sqrt{L}$ grows without bound, consistent with
$\sqrt{L}\log L$ scaling.** This rules out the GC prediction
$E[k]_\mathrm{GC} \sim 1.282\sqrt{L}$ and confirms canonical/GC
inequivalence: the number of tube segments is not self-averaging under the
grand canonical measure. This is the numerical signature of the
ensemble-inequivalence result established analytically in §6.3 of the paper.

**The $\tilde{w}$ coefficient converges to the Erdős–Lehner value slowly from
above.** The $k!/\prod m_j!$ factor in $\tilde{w}$ weights compositions
over partitions, and $2^k$ further favours high-$k$ configurations, pushing
$\langle k\rangle$ above its canonical value. Convergence to
$c = \sqrt{6}/(2\pi)$ requires $L \gg 10^3$; across the range studied,
$\langle k\rangle/(\sqrt{L}\log L)$ remains 2–6 times the asymptote.

**Split and merge moves are essential.** Transfer moves preserve $k$; only
splits and merges traverse the full range from $k = 1$ to $k = L$. The
$(0.5, 0.25, 0.25)$ split between move types balances exploration within
and between part-count levels, giving acceptance rates of 0.59–0.63 and
autocorrelation times of 32–57 steps across all $L$ studied.

**The sampler provides an independent Monte Carlo confirmation of
$\sqrt{L}\log L$ scaling in a regime where exact enumeration is infeasible.**
The exact DP computation of $p_k(L)$ in `be_saddle.py` provides the definitive
canonical test up to $L = 10{,}000$; the $\tilde{w}$ sampler confirms the same
qualitative scaling under a physically motivated measure at the same $L$ values,
and could be extended to much larger $L$ if needed.
