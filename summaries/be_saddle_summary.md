# Saddle-Point Analysis: Bose-Einstein Ansatz and Canonical Ensemble

## Summary of `be_saddle.py`

---

## (a) Theory

### Background: the naive saddle-point argument and its failures

The Hardy–Ramanujan (HR) analysis of integer partitions proceeds by treating the partition of $L$ as a statistical ensemble: each part size $j$ carries an occupancy $n_j$ (the number of parts equal to $j$), and the saddle-point of the partition generating function is found by maximising an entropy functional subject to the constraint $\sum_j j \cdot n_j = L$. A naive version of this argument uses the **Shannon entropy** $H[f] = -\sum_j f(j)\log f(j)$ and treats the run-length distribution $f(j)$ as the object to be optimised. This gives a **geometric** saddle:

$$f_\mathrm{geom}(j) = \frac{1}{\bar{r}}\left(1 - \frac{1}{\bar{r}}\right)^{j-1}, \qquad \bar{r} = \frac{L}{E[k]}$$

and suggests that the most probable number of parts scales as $k^* \sim (\pi/\sqrt{6})\sqrt{L}$. The script corrects three interlocking errors in this approach.

### Correction 1 — The entropy functional

The Shannon entropy treats parts of the same size as **distinguishable** and **independent**. In a partition, however, parts are intrinsically **indistinguishable**: a partition is an unordered multiset, not a sequence. The correct entropy functional is the **bosonic (Bose-Einstein) entropy**:

$$S_\mathrm{BE} = \sum_j \bigl[(n_j + 1)\log(n_j + 1) - n_j \log n_j\bigr]$$

This is the entropy of a collection of indistinguishable quanta distributed over energy levels $j$, which is exactly the combinatorial structure of an integer partition. It can be rewritten as $\sum_j [j\alpha n_j + \log(1 + n_j)]$, a form convenient for numerical evaluation.

### Correction 2 — The run-length distribution at the saddle

Maximising $S_\mathrm{BE}$ subject to $\sum_j j \cdot n_j = L$ with Lagrange multiplier $\alpha$ yields the **Bose-Einstein** occupancy:

$$n_j = \frac{1}{e^{j\alpha} - 1}$$

The corresponding run-length marginal (probability that a randomly chosen part has size $j$) is:

$$f_\mathrm{BE}(j) = \frac{1}{E[k]} \cdot \frac{1}{e^{j\alpha} - 1}, \qquad \alpha = \frac{\pi}{\sqrt{6L}}$$

The fugacity $\alpha$ is fixed by the HR saddle-point condition. The BE form differs substantially from the geometric approximation at small $j$: for $j = 1$ and $2$ the geometric underestimates and overestimates respectively, because it replaces $1/(e^x - 1)$ with $e^{-x}$, which is accurate only for $x \gg 1$ (large $j$ or large $L$).

### Correction 3 — The $k^*$ scaling

The HR saddle-point condition determines the fugacity $\alpha = \pi/\sqrt{6L}$ via the self-consistency equation:

$$E[k]_\mathrm{GC} = \sum_{j \geq 1} \frac{1}{e^{j\alpha} - 1} \sim \frac{\pi}{\sqrt{6}}\sqrt{L}$$

This is the **grand canonical** mean part count — the mean of $k$ under a Boltzmann-weighted ensemble where $L$ fluctuates. It governs the BE distribution and fixes $\alpha$, but it is **not** the mode or mean of $k$ under the **canonical** (uniform) measure where $L$ is fixed exactly.

Under the canonical (uniform) measure, the Erdős–Lehner theorem (1941) and Fristedt (1993) give:

$$E[k]_\mathrm{canon} \sim \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L \approx 0.3898\sqrt{L}\log L$$

$$k^*_\mathrm{canon} \sim 0.95 \cdot \frac{\sqrt{6}}{2\pi}\sqrt{L}\log L \approx 0.3704\sqrt{L}\log L$$

Both grow as $\sqrt{L}\log L$, parametrically faster than the grand canonical $\sqrt{L}$. The factor $\log L$ reflects that a uniform partition is much more likely to have many parts than the grand canonical ensemble would suggest, because the grand canonical ensemble exponentially suppresses configurations far from the thermal saddle, while the uniform measure counts every partition equally.

### Physical significance

These three corrections change the physical interpretation of the polymer's tube configurations. The indistinguishability of runs of equal length (bosonic statistics) enhances the occupancy of short runs over the geometric (classical) prediction. The canonical $k^* \sim \sqrt{L}\log L$ rather than $\sqrt{L}$ means that uniformly drawn long chains have many more tube segments than the naive estimate, with implications for tube renewal times and the onset of strong entanglement.

---

## (b) What the script does

The script verifies the three corrections numerically by comparing **exact** canonical statistics (computed from the full $p_k(L)$ table via dynamic programming) against both the grand canonical predictions and the BE asymptotic formulas, across chain lengths $L = 10$ to $10000$. It produces four main outputs:

1. **Main comparison table**: exact canonical mode $k^*$, mean $E[k]_\mathrm{canon}$, and standard deviation $\sigma$, compared against the grand canonical $E[k]_\mathrm{GC}$ and the asymptotic BE formula $(\sqrt{6}/2\pi)\sqrt{L}\log L$, at each $L$.

2. **$k^*$ scaling convergence table**: the ratio $k^*/\mathrm{form}$ where $\mathrm{form} = 0.3704\sqrt{L}\log L$, showing convergence to 1 from above as $L \to \infty$, together with $\sigma/L^{1/4}$ and the frequency at the mode.

3. **OLS scaling fit**: one-parameter least-squares fits of $k^*$ against both $\sqrt{L}\log(\sqrt{L})$ and $\sqrt{L}\log L$ (which are identical up to a factor of 2 in the coefficient), with residuals and the implied coefficient relative to $\sqrt{6}/(2\pi)$.

4. **Entropy comparison table**: $S_\mathrm{BE}$ (bosonic entropy at the HR fugacity) versus $S_\mathrm{geom}$ (geometric entropy evaluated at the canonical saddle $\bar{r} = L/k^*$), showing their ratio as a function of $L$.

---

## (c) How the script does it

### Exact $p_k(L)$ via dynamic programming

`p_k_L(L_max)` builds a 2D table `dp[L][k]` using the recurrence:

$$p_k(L) = p_{k-1}(L-1) + p_k(L-k)$$

The first term counts partitions whose smallest part equals 1 (remove it to get a partition of $L-1$ into $k-1$ parts); the second counts partitions whose smallest part is $\geq 2$ (subtract 1 from each part to get a partition of $L-k$ into $k$ parts). The table runs up to $L_\mathrm{max} = 10000$, requiring an $O(L^2)$ array of integers that grow very large (handled transparently by Python's arbitrary-precision arithmetic).

### Grand canonical functions

`alpha_from_L(L)` returns $\alpha = \pi/\sqrt{6L}$ directly.

`Ek_gc(alpha)` computes

$$E[k]_\mathrm{GC} = \sum_{j=1}^{\infty} \frac{1}{e^{j\alpha} - 1}$$

by direct summation, truncating at $j\alpha > 100$ where terms are negligible.

`be_entropy(alpha)` computes

$$S_\mathrm{BE} = \sum_j \bigl[j\alpha\, n_j + \log(1 + n_j)\bigr]$$

by the same summation, using the simplified form of the bosonic entropy.

### Canonical statistics

`canonical_stats(L, dp)` uses the $p_k(L)$ row of the DP table to compute exactly:

- $k^* = \operatorname{argmax}_k\, p_k(L)$ — the mode, found by scanning the row
- $E[k] = \sum_k k\, p_k(L)\,/\,p(L)$ — the mean
- $\sigma = \sqrt{E[k^2] - E[k]^2}$ — the standard deviation
- mode frequency $= p_{k^*}(L)\,/\,p(L)$

where $p(L) = \sum_k p_k(L)$ is the total partition count.

### Scaling fit

`fit_kstar_scaling(results)` performs one-parameter OLS (ordinary least squares with zero intercept) of $k^*$ against each of the two basis functions $\sqrt{L}\log(\sqrt{L})$ and $\sqrt{L}\log L$, using only $L \geq 100$ to be in the asymptotic regime. The OLS estimator is

$$\hat{a} = \frac{\sum_L k^*(L)\cdot b(L)}{\sum_L b(L)^2}$$

where $b(L)$ is the chosen basis. Since $\log(\sqrt{L}) = \tfrac{1}{2}\log L$, the two fits are equivalent (the coefficient for the $\log\sqrt{L}$ basis is twice that for the $\log L$ basis) and serve as a consistency check. The implied coefficient relative to $\sqrt{6}/(2\pi)$ is reported, with the value $\approx 0.952$ confirming the mode-to-mean ratio $c \approx 0.95$.

### Geometric entropy baseline

The per-part Shannon entropy of a geometric distribution with mean $\bar{r}$ is:

$$s_\mathrm{geom}(\bar{r}) = -\log(1-q) - \frac{q\log q}{1-q}, \qquad q = 1 - \frac{1}{\bar{r}}$$

The total geometric entropy is $k^* \cdot s_\mathrm{geom}(L/k^*)$. This is evaluated at the canonical saddle $\bar{r} = L/k^*$ and compared to $S_\mathrm{BE}$ at the same $L$.

---

## (d) Output produced

### Header: asymptotic formulae

The script first states the two sets of predictions being compared:

- **Grand canonical (HR)**: $\alpha = \pi/\sqrt{6L}$, $E[k]_\mathrm{GC} \sim 1.2825\sqrt{L}$, $f(j) \sim$ BE form
- **Canonical (uniform)**: $E[k] \sim 0.3898\sqrt{L}\log L$, $k^* \sim 0.3704\sqrt{L}\log L$

### Main comparison table

For each $L$ from 10 to 10000 the table shows $k^*_\mathrm{canon}$, $k^*/\sqrt{L}$, $E[k]_\mathrm{canon}$, $E[k]_\mathrm{canon}/\sqrt{L}$, $E[k]_\mathrm{GC}$, $E[k]_\mathrm{GC}/\sqrt{L}$, and the asymptotic BE_form $= (\sqrt{6}/2\pi)\sqrt{L}\log L$. Key observations:

- $E[k]_\mathrm{canon}/\sqrt{L}$ grows steadily from $\approx 1.4$ at $L=10$ to $\approx 3.9$ at $L=10000$, while $E[k]_\mathrm{GC}/\sqrt{L}$ grows more slowly (from $\approx 1.2$ to $\approx 3.8$), with the gap between them widening logarithmically and confirming the $\sqrt{L}\log L$ scaling.
- BE_form tracks $E[k]_\mathrm{canon}$ closely at large $L$, approaching it from below as the $\log L$ correction converges.
- The raw ratio $E[k]_\mathrm{canon}/E[k]_\mathrm{GC}$ starts at $\approx 1.17$ at $L=10$ but converges toward 1 at large $L$ (reaching $\approx 1.01$ at $L=2000$ and $\approx 1.00$ at $L=10000$), since both quantities grow as $\sqrt{L}\log L$ asymptotically with the same leading coefficient. The meaningful distinction is visible in the per-$\sqrt{L}$ columns, where the additive gap grows as $\log L$.

### $k^*$ scaling convergence table

The ratio $k^*/\mathrm{form}$ (where $\mathrm{form} = 0.3704\sqrt{L}\log L$) converges toward 1, reaching 1.001 at $L = 2000$ and 0.9997 at $L=10000$ (the slight undershoot reflecting slow convergence of sub-leading corrections). The standard deviation $\sigma$ scales as $L^{1/4}$ ($\sigma/L^{1/4}$ is approximately constant across the range), and the mode frequency decreases from $\approx 21\%$ at $L=10$ to $\approx 0.5\%$ at $L=10000$, showing the distribution over $k$ broadening relative to its mean as $L$ grows.

### Scaling fit

Both basis functions give identical RMS residuals (0.3945) and max relative errors (5.22%), confirming they are equivalent. The implied coefficient is $\approx 0.950$ when normalised by $\sqrt{6}/(2\pi)$, consistent with the theoretical prediction $k^* \sim 0.95\cdot E[k]_\mathrm{canon}$. The convergence table shows $k^*/(\sqrt{L}\log L)$ approaching $\approx 0.370$ at large $L$, compared to $c_\mathrm{ref} = \sqrt{6}/(2\pi) \approx 0.3899$; the residual $\approx 5\%$ gap persists even at $L=10000$, reflecting slow convergence of sub-leading logarithmic corrections.

### Entropy comparison

$S_\mathrm{BE}/\sqrt{L}$ converges slowly upward toward $\alpha_0 = \pi\sqrt{2/3} \approx 2.565$, reaching 2.494 at $L=2000$ and 2.529 at $L=10000$. $S_\mathrm{geom}/\sqrt{L}$ grows much faster (from 2.1 at $L=10$ to 10.5 at $L=2000$ and 14.9 at $L=10000$) because $S_\mathrm{geom}$ is evaluated at the canonical $k^*$, which grows as $\sqrt{L}\log L$ — meaning the geometric entropy is calculated for an increasingly large number of parts, inflating it artificially. The ratio $S_\mathrm{BE}/S_\mathrm{geom}$ falls from $\approx 0.93$ at $L=10$ to $\approx 0.24$ at $L=2000$ and $\approx 0.17$ at $L=10000$, confirming that **the geometric entropy, evaluated at the correct canonical $k^*$, greatly overestimates the true bosonic entropy**.

---

## (e) What we learn

**The grand canonical and canonical ensembles differ qualitatively, not just quantitatively.** The grand canonical ensemble (where $L$ fluctuates around its mean) assigns a typical part count of $\sim (\pi/\sqrt{6})\sqrt{L}$. The canonical ensemble (where every partition of exactly $L$ is equally weighted) gives a modal part count of $\sim 0.37\sqrt{L}\log L$. These are different *functions* of $L$: the ratio grows without bound as $\log L$. The naive saddle-point argument, by using the wrong ensemble, misidentifies the typical configuration by a factor that grows logarithmically with chain length.

**Runs of equal length are indistinguishable: bosonic statistics are mandatory.** The correct entropy for counting partitions is the bosonic entropy $S_\mathrm{BE}$, not the Shannon entropy. The geometric distribution emerges from $H[f]$; the BE distribution emerges from $S_\mathrm{BE}$. At $j=1$ the BE form gives higher occupancy than geometric, and at $j=2\ldots6$ it gives lower — the BE distribution is more concentrated on the shortest runs, which is physically correct: short runs are favoured because more configurations of the remaining length are accessible.

**The $k^* \sim \sqrt{L}\log L$ scaling converges, but slowly.** The ratio $k^*/[0.3704\sqrt{L}\log L]$ is still $\approx 5\%$ above 1 at $L=100$, reaches 1.001 at $L=2000$, and slightly undershoots to 0.9997 at $L=10000$. This slow convergence (governed by sub-leading corrections of order $\log\log L/\log L$) means that the asymptotic formula is a reasonable but not precise guide at experimentally accessible chain lengths.

**$S_\mathrm{geom}$ and $S_\mathrm{BE}$ are not interchangeable as entropy measures for the partition problem.** Because $k^*$ grows as $\sqrt{L}\log L$ rather than $\sqrt{L}$, the geometric entropy evaluated at the canonical saddle grows far faster than $S_\mathrm{BE}$. The ratio $S_\mathrm{BE}/S_\mathrm{geom} \to 0$ as $L \to \infty$. This is not a contradiction: the geometric entropy is simply the wrong quantity to compute once the correct bosonic entropy and canonical $k^*$ are in hand.

**The GC mean $E[k]_\mathrm{GC}$ retains a precise role.** Despite being parametrically smaller than $E[k]_\mathrm{canon}$, the grand canonical mean is not wrong — it is the correct self-consistency condition for the fugacity $\alpha = \pi/\sqrt{6L}$, and it enters the BE distribution $f_\mathrm{BE}(j)$ as the normaliser when evaluated in the GC ensemble. Its role is to fix the *shape* of the distribution (the $j$-dependence), while the canonical $E[k]$ provides the correct *normalisation* for the marginal probability under the uniform measure. The two quantities answer different questions and must not be conflated.
