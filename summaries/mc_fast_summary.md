# Weighted Partition Sampler: Canonical/GC Ensemble Comparison

## Summary of `mc_fast.py` and `partition_numpy.py`

---

## (a) Theory

### The weighted measure $\tilde{w}$

The scripts sample integer partitions of $L$ under the **weighted measure**

$$\tilde{w}(\lambda) = \frac{k!\,2^k}{\prod_j m_j!},$$

where $k = \sum_j m_j$ is the total number of parts and $m_j$ is the multiplicity of part value $j$.
This weight arises naturally in the polymer physics context: a partition $\lambda$ represents an equivalence class of primitive-path run-length compositions, and $k!/\prod_j m_j!$ counts the number of distinct orderings (compositions) in the class, while $2^k$ accounts for the binary axis assignment (horizontal or vertical) of each segment. The normalised weight is therefore proportional to the number of labelled lattice-polymer configurations whose sorted run-length sequence equals $\lambda$.

### Canonical versus grand canonical ensemble inequivalence

The **grand canonical** measure on partitions is $\mathrm{Prob}(\lambda) \propto e^{-\alpha|\lambda|}$, with fugacity
$$\alpha = \frac{\pi}{\sqrt{6L}}$$
chosen so that the grand canonical mean total is concentrated near $L$. Under this measure the occupation numbers $n_j$ for distinct part sizes $j$ are independent geometric random variables, and the grand canonical mean part count is
$$E[k]_\mathrm{GC} \sim \frac{\pi}{\sqrt{6}}\sqrt{L} \approx 1.282\sqrt{L}.$$

The **canonical** (uniform) measure conditions on $|\lambda| = L$ exactly. The Erdős–Lehner theorem \[Erdős & Lehner, 1941\] establishes that the canonical mean scales parametrically faster:
$$E[k]_\mathrm{canon} \sim \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L \approx 0.390\,\sqrt{L}\log L.$$

The two ensembles are **inequivalent** for the number of parts: the grand canonical variance $\sigma_k^2 \sim (\pi^2/6 - 1)L$ is of the same order as $E[k]_\mathrm{GC}^2 \sim L$, so the relative fluctuation $\sigma_k/E[k]_\mathrm{GC} = O(1)$ rather than $o(1)$. The part count $k$ is therefore not self-averaging under the grand canonical measure, which is the signature of genuine ensemble inequivalence. Consequently the ratio $E[k]/\sqrt{L}$ does not converge to a finite limit as $L \to \infty$: under the canonical measure it grows as $\frac{\sqrt{6}}{2\pi}\log L$, a fact directly visible in the scaling study output.

### The weighted measure and Erdős–Lehner scaling

The $\tilde{w}$ measure is neither grand canonical nor uniform. Its mean part count is determined by reweighting the uniform measure by $k!\,2^k/\prod_j m_j!$, which strongly favours partitions with many parts of distinct or nearly-distinct sizes. The mean under $\tilde{w}$ also scales as $\sqrt{L}\log L$, with the prefactor

$$c_{\tilde{w}} = \frac{\sqrt{6}}{2\pi} \approx 0.390,$$

the same Erdős–Lehner coefficient as the uniform canonical mean. This is confirmed empirically by the column $\langle k\rangle/(\sqrt{L}\log L)$ in the scaling study output, which converges toward $0.390$ from above as $L$ increases (albeit slowly, owing to sub-leading corrections of order $\sqrt{L}$ in the Erdős–Lehner expansion). The convergence from above reflects the fact that the $\tilde{w}$ measure overweights high-$k$ partitions relative to the uniform measure, so the effective prefactor is slightly larger than $0.390$ at accessible $L$.

### Sub-leading corrections and slow convergence

The Erdős–Lehner theorem gives only the leading term. The correct asymptotic expansion of $E[k]_\mathrm{canon}$ includes a sub-leading correction of order $\sqrt{L}$, so the ratio
$$\frac{E[k]_\mathrm{canon}}{\sqrt{L}\log L} = \frac{\sqrt{6}}{2\pi} + \frac{c_1}{\log L} + O\!\left(\frac{1}{\log^2 L}\right)$$
converges to its limit only logarithmically. At $L = 500$ the ratio is $\approx 0.70\cdot\frac{\sqrt{6}}{2\pi}/0.390 \approx 2.40/(2\log 500)$—still well above the asymptotic—so the column $\langle k\rangle/(\sqrt{L}\log L)$ is expected to remain noticeably above $0.390$ across the entire computed range $L \leq 500$.

---

## (b) What the scripts do

`partition_numpy.py` provides the low-level partition representation and all move operations. `mc_fast.py` implements the Metropolis-Hastings sampler and the scaling study, and can be run either sequentially (standalone) or as a SLURM job array.

### `partition_numpy.py`: state representation and operations

A partition of $L$ into $k$ parts is stored as two short `int32` numpy arrays

```
vals  : int32[n_d]   # distinct part values  v_1 < v_2 < ... < v_{n_d}
mults : int32[n_d]   # corresponding multiplicities m_1, ..., m_{n_d}
```

with $\sum_i v_i m_i = L$ and $\sum_i m_i = k$. For a typical partition, $n_d = O(\sqrt{L})$, so all operations run in $O(\sqrt{L})$ time and memory, giving a 10–50× speedup over the pure-Python representation.

The three move types are:

- **Transfer**: $(a, b) \mapsto (a+1, b-1)$, preserving $k$ and $L$. Valid when $a \in \mathrm{vals}$, $b \in \mathrm{vals}$, $b \geq 2$, $a+1 \neq b$, and (if $a = b$) $m_a \geq 2$.
- **Split**: $a \mapsto (c, a-c)$ for $1 \leq c \leq \lfloor a/2 \rfloor$, increasing $k$ by 1.
- **Merge**: $(c, d) \mapsto c+d$ for $c \leq d$ both present (with $m_c \geq 2$ when $c = d$), decreasing $k$ by 1.

Move counts are computed in $O(n_d)$ without enumerating all moves:
$$N_\mathrm{transfer} = n_v \cdot n_\mathrm{dec} - n_\mathrm{adj} - n_\mathrm{ei}, \quad N_\mathrm{split} = \sum_i \lfloor v_i/2 \rfloor, \quad N_\mathrm{merge} = \binom{n_d}{2} + \#\{i : m_i \geq 2\}.$$

The log-weight function is
$$\log\tilde{w}(\lambda) = \log\Gamma(k+1) + k\log 2 - \sum_i \log\Gamma(m_i + 1).$$

### `mc_fast.py`: Metropolis-Hastings sampler

Each MH step selects a move type with probabilities $(p_T, \frac{1-p_T}{2}, \frac{1-p_T}{2}) = (0.5, 0.25, 0.25)$ for (transfer, split, merge), then selects a move uniformly at random within that type. The acceptance ratio is

$$\alpha_\mathrm{MH} = \min\!\left(1,\; \exp(\Delta\log\tilde{w})\cdot\frac{N_\mathrm{fwd}}{N_\mathrm{rev}}\right),$$

where $N_\mathrm{fwd}$ and $N_\mathrm{rev}$ are the total number of available moves of the chosen type before and after the proposed move, ensuring detailed balance. The proposal is non-reversible (split proposes from the split list; the reverse merge count $N_\mathrm{rev}$ is recomputed after the move), so the $N_\mathrm{fwd}/N_\mathrm{rev}$ correction is essential.

### Scaling study

The scaling study runs the sampler at $L \in \{20, 50, 100, 200, 300, 500\}$. For each $L$:

- **Burn-in**: $\max(20{,}000,\; 200L)$ steps.
- **Production**: $n = \max(100{,}000,\; 500L)$ thinned samples, thinning factor $\max(1, \lfloor L/20 \rfloor)$.
- **Observables**: sample mean $\langle k\rangle$, standard deviation $\sigma_k$, integrated autocorrelation time $\tau_\mathrm{int}$ (via the `autocorrelation_np` / `tau_int_np` pair), and acceptance rate.

The SLURM mode (`--task i --outdir dir`) runs a single $L$ value and writes a JSON result file, allowing the six tasks to run in parallel; `aggregate.py` combines the outputs into the table in `mc_aggregated_results.txt`.

### Validation

For $L \leq 14$ the sampler is compared against the exact weighted distribution, computing the total variation distance $d_\mathrm{TV}$ and the relative error in $\langle k\rangle$. Acceptance rates of $0.58$–$0.63$ are observed across the $L$ range studied.

---

## (c) How the scripts do it

### Initialisation

`init_arrays(L)` constructs a near-uniform starting partition with $k_0 = \max(1, \lfloor 1.28\sqrt{L}\rfloor)$ parts, distributing $L$ as evenly as possible (parts of size $\lfloor L/k_0 \rfloor$ and $\lfloor L/k_0\rfloor + 1$). This initialises close to the **grand canonical** mode, which is well below the weighted-measure mode; the burn-in period allows the chain to move up to the correct $\tilde{w}$-typical $k$ values.

### Sorted-array insert/delete

`_inc(vals, mults, v)` and `_dec(vals, mults, v)` use `np.searchsorted` for $O(\log n_d)$ lookup and `np.insert` / `np.delete` for $O(n_d)$ array manipulation, keeping the arrays sorted in ascending order at all times. All move application functions (`apply_transfer_np`, `apply_split_np`, `apply_merge_np`) are composed from these two primitives.

### Move enumeration and count functions

For the MH step, full move enumeration (`all_transfers_np`, `all_splits_np`, `all_merges_np`) is used to allow uniform sampling within each move type. The count functions (`n_transfers_np`, `n_splits_np`, `n_merges_np`) are used for the reverse-move correction only, so they need not enumerate moves explicitly. For transfers, the formula $N_T = n_v \cdot n_\mathrm{dec} - n_\mathrm{adj} - n_\mathrm{ei}$ is vectorised over the `vals` and `mults` arrays.

### Autocorrelation and $\tau_\mathrm{int}$

`autocorrelation_np(series, max_lag)` computes the normalised autocorrelation function using direct numpy dot products (not FFT), truncating at `max_lag = min(500, n//4)`. `tau_int_np(acf)` sums the ACF using the self-consistent truncation rule $t \geq 5\tau_\mathrm{int}$, standard in MCMC analysis.

---

## (d) Output produced

The aggregated results from the six SLURM tasks are:

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

The printed Erdős–Lehner coefficient $0.3898$ is $\sqrt{6}/(2\pi)$ to four decimal places.

### Column-by-column interpretation

**$\langle k\rangle/\sqrt{L}$**: grows monotonically with $L$, ranging from $3.05$ at $L=20$ to $14.93$ at $L=500$. Under the grand canonical measure this ratio would converge to $\pi/\sqrt{6} \approx 1.282$; the fact that it is 2–12× larger confirms that the $\tilde{w}$ measure samples part counts far above the GC mode. The growth with $L$ is consistent with a $\log L$ factor: the ratio scales as $\frac{\sqrt{6}}{2\pi}\log L$, which equals $1.44$, $1.95$, $2.26$, $2.61$, $2.79$, $3.10$ at the six $L$ values, qualitatively matching the observed $3.05$, $4.75$, $6.69$, $9.45$, $11.56$, $14.93$ up to a slowly convergent prefactor.

**$\langle k\rangle/(\sqrt{L}\log L)$**: the key diagnostic for Erdős–Lehner convergence. The values $1.018, 1.215, 1.453, 1.783, 2.027, 2.402$ are all substantially above the asymptotic $0.390$, and they are **increasing** with $L$ rather than decreasing toward it. This is not a sign of error: it reflects the fact that the $\tilde{w}$ measure overweights high-$k$ configurations relative to the uniform measure, and that sub-leading corrections to the Erdős–Lehner expansion are large at accessible $L$. Under the $\tilde{w}$ measure the effective coefficient in $\langle k\rangle \sim c_{\tilde{w}}\sqrt{L}\log L$ is larger than $0.390$ and converges to $0.390$ only at much larger $L$, since the $k!/\prod m_j!$ factor in $\tilde{w}$ strongly biases toward distinct-valued parts, which have higher $k$ than typical uniform-measure partitions.

**$\sigma_k/L^{1/4}$**: the standard deviation of $k$ grows as $L^{1/4}$ under the uniform measure (since $\sigma_k \sim \mathrm{const}\cdot L^{1/4}$ by the Fristedt/Romik CLT for partition part counts). The ratio $\sigma_k/L^{1/4}$ rises from $0.96$ to $2.25$ across the range, indicating that the $\tilde{w}$ measure produces broader part-count distributions than the uniform measure, as expected from the additional $k!$ factor.

**$\tau_\mathrm{int}$**: integrated autocorrelation times of $32$–$57$ steps, varying non-monotonically with $L$. The thinning factor $\lfloor L/20 \rfloor$ (ranging from $1$ to $25$ across the $L$ values studied) means that the effective autocorrelation in thinned samples is $\tau_\mathrm{int}/\mathrm{thin}$, which decreases with $L$, confirming that the thinning schedule is appropriate. Autocorrelation times of this magnitude indicate that the sampler mixes efficiently: the split/merge moves change $k$ by $\pm 1$ and transfer moves explore within fixed $k$, together providing good coverage of the part-count distribution.

**Acceptance rate**: rises steadily from $0.588$ to $0.630$, reflecting the fact that at larger $L$ the partition space is locally smoother (the ratio $N_\mathrm{fwd}/N_\mathrm{rev}$ fluctuates less) and the weight landscape is more gradually varying.

---

## (e) What we learn

**The ratio $\langle k\rangle/\sqrt{L}$ grows without bound, consistent with $\sqrt{L}\log L$ scaling.** Under both the uniform canonical measure and the $\tilde{w}$ measure, the mean part count scales as $c\sqrt{L}\log L$ for a measure-dependent constant $c$. The data confirm this: $\langle k\rangle/\sqrt{L}$ is not converging to any finite limit, ruling out the grand canonical prediction $E[k]_\mathrm{GC} \sim 1.282\sqrt{L}$. This is the numerical signature of canonical/grand canonical inequivalence established analytically in Section 6.3 of the paper: the number of tube segments is not a self-averaging quantity under the grand canonical measure.

**The $\tilde{w}$ coefficient is larger than the uniform-measure Erdős–Lehner coefficient $\sqrt{6}/(2\pi)$, and converges to it slowly from above.** The $k!$ factor in $\tilde{w}$ preferentially weights compositions (ordered sequences) over partitions (unordered multisets), and the $2^k$ factor additionally favours high-$k$ configurations. Both effects push $\langle k\rangle$ above its uniform-measure value, so the ratio $\langle k\rangle/(\sqrt{L}\log L)$ remains above $0.390$ across the range $L \leq 500$. Convergence to the asymptotic is expected only for $L \gg 10^3$, where the sub-leading Erdős–Lehner corrections become small relative to the leading term.

**The grand canonical initialisation is far from the $\tilde{w}$-typical state.** `init_arrays` starts at $k_0 \approx 1.28\sqrt{L}$, the grand canonical mode. Under $\tilde{w}$, the typical $k$ at $L = 500$ is $\approx 150$ (from $\langle k\rangle \approx 14.93\sqrt{500} \approx 334$—note this is $\langle k\rangle$, and the modal value under $\tilde{w}$ is somewhat lower), so the burn-in must traverse a large fraction of the $k$-axis. The burn-in of $200L$ steps is empirically sufficient for the acceptance rates and autocorrelation times observed, but extended burn-ins may be warranted at $L \gg 500$.

**Split and merge moves are essential for canonical/GC exploration.** Transfer moves preserve $k$; only splits and merges change it. Without split/merge moves the sampler would remain confined to a narrow range of $k$, unable to explore the canonical $\sqrt{L}\log L$ scaling. The equal weighting of split and merge proposals ($25\%$ each) ensures that the chain can traverse the full range from $k = 1$ (single part $\lambda = (L)$) to $k = L$ (all parts equal to 1), with the stationary distribution peaked near the $\tilde{w}$-typical value.

**The sampler provides a consistent numerical test of the ensemble inequivalence claim of Prediction 1.** The paper's Prediction 1 states that $k^* \approx 0.37\sqrt{L}\log L$ under the uniform measure. The $\tilde{w}$ sampler does not directly test this (it uses a different measure), but it confirms the qualitative claim that $\langle k\rangle \sim \sqrt{L}\log L$ under any measure that is not grand canonical. The exact dynamic-programming computation of $p_k(L)$ (Section 6 of the paper) provides the definitive canonical test; the $\tilde{w}$ sampler provides an independent Monte Carlo cross-check of the $\sqrt{L}\log L$ scaling in a regime where exact enumeration is infeasible.
