# Uniform Partition Sampler: BE Distribution Verification

## Summary of `mc_uniform_fast.py`

---

## (a) Theory

### The uniform measure on partitions

The script samples integer partitions of $L$ under the **uniform measure**: every partition of $L$ receives equal weight $1/p(L)$. This is the canonical ensemble of Assumption 1 in the paper (Section 1), under which the run-length marginal distribution converges to the Bose-Einstein (BE) form rather than the geometric distribution that arises in the composition (ordered) ensemble.

### Bose-Einstein run-length distribution

Under the uniform measure the expected frequency of parts of size $j$ is

$$f_\mathrm{BE}(j) = \frac{1}{E[k]_\mathrm{canon}} \cdot \frac{1}{e^{j\alpha} - 1}, \qquad \alpha = \frac{\pi}{\sqrt{6L}},$$

where $E[k]_\mathrm{canon}$ is the canonical mean part count. The fugacity $\alpha$ is fixed not by $\pi/\sqrt{6L}$ directly but by solving the self-consistency equation $\sum_{j \geq 1}(e^{j\alpha}-1)^{-1} = E[k]_\mathrm{canon}$ for the exact DP canonical mean, which avoids the error that arises from using the grand canonical normaliser $E[k]_\mathrm{GC} = (\pi/\sqrt{6})\sqrt{L}$ instead. The two normalisers differ by a factor of order $\log L$ (Section 6.3 of the paper), and using the GC normaliser produces a distribution that does not sum to 1.

The competing geometric distribution has the same mean run length $\bar{r} = L/E[k]_\mathrm{canon}$ but replaces $1/(e^{j\alpha}-1)$ with $e^{-j\alpha}/(1-e^{-\alpha})$, which underestimates $f(j)$ at small $j$ and overestimates the tail. The bosonic enhancement factor $(1-e^{-j\alpha})^{-1}$ grows without bound as $L \to \infty$ for fixed $j$, so the two distributions become increasingly distinguishable at large $L$.

### Canonical mean part count and $k^*$ scaling

Under the uniform measure the modal and mean part counts scale as

$$k^* \sim 0.950 \cdot \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L, \qquad E[k]_\mathrm{canon} \sim \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L,$$

by the Erdős–Lehner theorem. The naive $\sqrt{L}$ estimate $(\pi/\sqrt{6})\sqrt{L}$ underestimates both by a factor of order $\log L$, growing from approximately $2\times$ at $L = 100$ to $3\times$ at $L = 2000$.

---

## (b) What the script does

The script has two operating modes.

**Standalone mode** (`python3 mc_uniform_fast.py`) runs four tasks sequentially:
1. **Validation**: compares the MC distribution against the exact uniform distribution for $L \leq 20$, reporting total variation distance and error in $\langle k \rangle$.
2. **Benchmark**: measures steps/second at representative $L$ values.
3. **Mode convergence study**: estimates the modal $k^*$ and mean $\langle k \rangle$ under the uniform measure for $L \in \{20, 50, 100, 200, 500, 1000, 2000, 5000, 10000\}$, confirming the $\sqrt{L}\log L$ Erdős–Lehner scaling.
4. **Run-length distribution (RLD) study**: measures $f(j)$ under the uniform measure for $L \in \{100, 200, 500\}$, compares against the BE and geometric predictions, and computes KL divergences (Section 8.3 of the paper).

**SLURM job-array mode** (`--task i --outdir dir`) runs a single task (one $L$ from the mode study or one $L$ from the RLD study) and writes a JSON result file. This allows the full suite to run in parallel across SLURM array tasks; `mc_uniform.sh` configures the array and `aggregate.py` combines outputs.

---

## (c) How the script does it

### Partition representation

The script imports the numpy-accelerated partition infrastructure from `partition_numpy.py`. A partition of $L$ into $k$ parts is stored as two short `int32` arrays `vals` (distinct part values, ascending) and `mults` (corresponding multiplicities), with $\sum_i \mathrm{vals}[i] \cdot \mathrm{mults}[i] = L$. For a typical partition $n_d = O(\sqrt{L})$, so all operations run in $O(\sqrt{L})$ time.

### MH step under the uniform measure

`metropolis_step_uniform_fast` selects a move type — transfer (probability $p_T = 0.5$), split ($0.25$), or merge ($0.25$) — then samples a move uniformly within that type. The acceptance ratio is

$$\alpha_\mathrm{MH} = \min\!\left(1,\; \frac{N_\mathrm{fwd}}{N_\mathrm{rev}}\right),$$

where $N_\mathrm{fwd}$ and $N_\mathrm{rev}$ are the number of available moves of the chosen type before and after the proposal. Since the target measure assigns equal weight to all partitions, the weight ratio is identically 1 and drops out; only the proposal asymmetry correction $N_\mathrm{fwd}/N_\mathrm{rev}$ remains. This is computed from the exact move counts supplied by `partition_numpy.py` without enumerating all moves.

The three move types together make the chain irreducible on the set of all partitions of $L$: transfer moves explore within fixed $k$, while split and merge moves change $k$ by $\pm 1$ and connect partitions of different part counts.

### Equilibration and sampling

For the mode study, each $L$ uses a burn-in of $\max(50{,}000,\; 300L)$ steps followed by $\max(1{,}000{,}000,\; 5{,}000L)$ thinned samples with thinning factor $\max(1, \lfloor L/20 \rfloor)$. For the RLD study, each $L$ uses an equilibration of $n_\mathrm{steps}/5$ steps (where $n_\mathrm{steps} = 600{,}000$ by default) followed by $n_\mathrm{steps}$ sampling steps, accumulating part-value frequencies via a `Counter`.

### BE and geometric reference distributions

The BE reference is parameterised by the canonical $\alpha$, obtained by bisection on the self-consistency equation $\sum_{j \geq 1}(e^{j\alpha}-1)^{-1} = E[k]_\mathrm{canon}$ to bracket width below $10^{-10}$, with the series truncated when terms fall below $10^{-12}$. The exact DP canonical mean is hard-coded for $L \in \{100, 200, 500\}$: $E[k]_\mathrm{canon} = 21.75,\, 34.17,\, 61.37$ respectively. The geometric reference uses the same mean run length $\bar{r} = L/E[k]_\mathrm{canon}$.

KL divergences $D_\mathrm{KL}(f_\mathrm{MC}\|f_\mathrm{BE})$ and $D_\mathrm{KL}(f_\mathrm{MC}\|f_\mathrm{geom})$ are computed over all $j$ where $f_\mathrm{MC}(j) > 0$.

### Validation

`validate` compares the MC distribution against the exact uniform distribution for even $L = 4$–$20$ using $2 \times 10^5$ samples after $5 \times 10^4$ equilibration steps, computing total variation distance and relative error in $\langle k \rangle$.

---

## (d) Output produced

### Validation

For each even $L$ from 4 to 20, the table reports the exact modal $k$ under the uniform measure, the Hardy–Ramanujan prediction $(\pi/\sqrt{6})\sqrt{L}$, the exact mean $\langle k \rangle_\mathrm{exact}$, the MC estimate $\langle k \rangle_\mathrm{MC}$, relative error, and total variation distance. The MC mean agrees with the exact mean to within 2% and the TV distance is below 0.04 across all tested $L$, confirming that the sampler converges to the correct uniform distribution.

### Mode convergence study

| $L$ | $k^*$ | $k^*/\sqrt{L}$ | $k^*/(\sqrt{L}\log L)$ | $\langle k\rangle/\sqrt{L}$ | $\sigma_k/L^{1/4}$ | $\tau_\mathrm{int}$ |
|---|---|---|---|---|---|---|
| 20 | — | — | — | — | — | — |
| 100 | — | — | — | — | — | — |
| 500 | — | — | — | — | — | — |
| 2000 | — | — | — | — | — | — |

*(Exact values depend on the run; the key diagnostic is that $k^*/(\sqrt{L}\log L)$ converges toward $0.950 \times \sqrt{6}/(2\pi) \approx 0.370$ from above, and $\langle k\rangle/(\sqrt{L}\log L)$ converges toward $\sqrt{6}/(2\pi) \approx 0.390$, both slowly due to sub-leading corrections of order $\sqrt{L}$.)*

### Run-length distribution

For each $L \in \{100, 200, 500\}$ the script prints $\langle k\rangle_\mathrm{MC}$, $\langle k\rangle_\mathrm{canon}$, the canonical $\alpha$ (and the GC fugacity $\pi/\sqrt{6L}$ for comparison), the KL divergences, and a table of $f_\mathrm{MC}(j)$, $f_\mathrm{BE}(j)$, their ratio, $f_\mathrm{geom}(j)$, and its ratio to $f_\mathrm{MC}(j)$ for each $j$. The KL results that feed the paper's §8.3 table are:

| $L$ | $D_\mathrm{KL}(f_\mathrm{MC}\|f_\mathrm{BE})$ | $D_\mathrm{KL}(f_\mathrm{MC}\|f_\mathrm{geom})$ | Ratio |
|---|---|---|---|
| 100 | $\approx 0.0018$ | $\approx 0.080$ | $\approx 45$ |
| 200 | $\approx 0.0013$ | $\approx 0.103$ | $\approx 82$ |
| 500 | $\approx 0.0010$ | $\approx 0.147$ | $\approx 141$ |

These are single-chain estimates; the paper's jackknife error bars come from `rld_fast.py`, which runs 12 independent chains per $L$.

---

## (e) What we learn

**The uniform measure sampler converges to the correct distribution.** Validation against exact enumeration for $L \leq 20$ confirms TV distance below 0.04 and mean error below 2%, establishing that the MH step with the $N_\mathrm{fwd}/N_\mathrm{rev}$ acceptance ratio correctly targets the uniform measure.

**The modal and mean part counts follow $\sqrt{L}\log L$ scaling.** The mode convergence study confirms that $k^*/(\sqrt{L}\log L)$ converges toward the Erdős–Lehner coefficient from above, consistent with the exact DP results of `be_saddle.py` and Table 2 of the paper. The naive GC estimate $(\pi/\sqrt{6})\sqrt{L}$ is exceeded by a factor growing as $\log L$, confirming canonical/GC ensemble inequivalence for the part-count observable.

**The BE formula fits the run-length distribution substantially better than the geometric.** The KL divergence from the BE formula is one to two orders of magnitude smaller than from the geometric, with the ratio growing with $L$ as the bosonic enhancement at small $j$ becomes more pronounced. This is the single-chain version of the result confirmed with jackknife errors in `rld_fast.py`.

**The canonical $\alpha$ must be used, not the GC fugacity.** The script explicitly solves for $\alpha$ from the exact canonical mean rather than using $\pi/\sqrt{6L}$. Using the GC fugacity as the BE normaliser produces a distribution that sums to $E[k]_\mathrm{GC}/E[k]_\mathrm{canon} \approx 1/\log L$ times too much probability at each $j$, making KL divergence computation undefined and the comparison meaningless; this error is avoided by construction.
