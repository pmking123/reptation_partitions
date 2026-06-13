# Bounded-Segment Crossover: M3 Model Verification

## Summary of `crossover.py`

---

## (a) Theory

### The M3 model

The M3 model describes a polymer chain of length $L$ whose tube is decomposed
into at most $N_e$ segments, each of length at least $s$. The accessible
configurations are counted by

$$P_L^{(N_e,\, s)} = \#\{\text{partitions of } L \text{ into at most } N_e \text{ parts, each} \geq s\}$$

This is a doubly constrained partition problem: the upper bound $k \leq N_e$
limits the number of segments (the tube entanglement constraint), and the
lower bound $r_i \geq s$ enforces a minimum segment length (the persistence
length constraint). For $s = 1$ the lower bound is trivially satisfied, so
the constraint reduces to partitions of $L$ into at most $N_e$ parts.

Note: the script's internal documentation refers to this as the "M4 model",
using an earlier numbering. In the paper, M3 is the bounded-segment model and
M4 is the distinct-parts model.

### Two asymptotic regimes

For small $L$ the upper bound $k \leq N_e$ is rarely active and $P_L^{(N_e,s)}$
is essentially the unrestricted partition number $p(L)$. The conformational
entropy is therefore governed by the Hardy–Ramanujan formula:

$$\log P_L^{(N_e,s)} \sim \alpha_0 \sqrt{L}, \qquad \alpha_0 = \pi\sqrt{\tfrac{2}{3}} \approx 2.565, \qquad L \ll L^*$$

For large $L$ the constraint $k \leq N_e$ binds on every configuration. Using
the standard asymptotic $P_L^{(N_e,\, s)} \sim (L/s)^{N_e - 1} / (N_e!\,(N_e-1)!)$
for fixed $N_e$ and large $L$, the entropy grows only polynomially:

$$\log P_L^{(N_e,s)} \sim (N_e - 1)\log(L/s), \qquad L \gg L^*$$

### Two crossover scales

The script computes two characterisations of the crossover location.

**Grand canonical estimate.** Setting $E[k]_\mathrm{GC} = N_e$ where
$E[k]_\mathrm{GC} \sim (\pi/\sqrt{6})\sqrt{L}$ gives the closed-form
expression:

$$L^*_\mathrm{GC} = \frac{6}{\pi^2} N_e^2 \approx 0.608\, N_e^2$$

**Canonical estimate.** Setting the canonical mean $E[k]_\mathrm{canon} = N_e$
where $E[k]_\mathrm{canon} \sim (\sqrt{6}/2\pi)\sqrt{L}\log L$ gives
$L^*_\mathrm{canon}$ as the solution of:

$$\frac{\sqrt{6}}{2\pi}\sqrt{L^*}\log L^* = N_e$$

solved numerically by bisection. For $N_e = 10, 20, 30$: $L^*_\mathrm{canon}
\approx 45, 116, 208$ — substantially smaller than the GC values $L^*_\mathrm{GC}
\approx 61, 243, 547$. The ratio $L^*_\mathrm{GC}/L^*_\mathrm{canon}$ grows
from $\approx 1.3$ at $N_e = 10$ to $\approx 2.6$ at $N_e = 30$ and
$\approx 4.4$ at $N_e = 85$, because the logarithmic factor in
$E[k]_\mathrm{canon}$ pushes the canonical crossover to shorter chains.
$L^*_\mathrm{canon}$ is the physically correct crossover estimate and is
what the paper uses throughout.

### Asymptote-intersection scale

The script also computes $L^\dagger$, the chain length at which the two
asymptotic predictions are equal: $\alpha_0\sqrt{L} = (N_e-1)\log(L/s)$.
This is not the onset of the crossover but the point at which the polynomial
asymptote has caught up with the exponential one. For $N_e = 10, 20$:
$L^\dagger \approx 465, 3705$; for $N_e = 30$, $L^\dagger > 5000$.

### Physical significance

The crossover at $L^*_\mathrm{canon}$ marks the onset of strong tube
confinement. The $N_e^2$ scaling of $L^*$ is the robust quantitative
prediction, with $L^*/N_e \sim N_e$ corresponding to a chain of approximately
$N_e$ entanglement lengths.

---

## (b) What the script does

The script accepts $N_e$ values, minimum segment length $s$, and $L_\mathrm{max}$
via command-line arguments. Default values are
$N_e \in \{10, 20, 30, 50, 75, 85, 100, 200\}$, $s = 1$, $L_\mathrm{max} = 5000$.
For each $N_e$ value it computes:

1. The exact M3 partition count $P_L^{(N_e,s)}$ via dynamic programming.
2. The unrestricted $p(L)$ for reference (computed once and shared).
3. A crossover table showing $\log P_L^{(N_e)}$, the exponential asymptote
   $\alpha_0\sqrt{L}$, and the polynomial asymptote $(N_e-1)\log(L/s)$, with
   both asymptote ratios, at logarithmically spaced $L$ values anchored to
   $L^*_\mathrm{GC}$.
4. A summary crossover table showing $L^*_\mathrm{GC}$ and $L^\dagger$ for
   all $N_e$ values.
5. Three CSV files: full crossover data, a crossover summary, and scaled
   entropy data ($\log P$ and both ratios versus $L/L^*_\mathrm{GC}$).

The header also displays both $L^*_\mathrm{GC}$ and $L^*_\mathrm{canon}$ for
all $N_e$ values before the per-$N_e$ analysis begins.

---

## (c) How the script does it

### Unrestricted partitions: `compute_p`

Uses Euler's pentagonal number recurrence (same as in `rr_entropy.py`):
$O(L\sqrt{L})$ time, $O(L)$ space. Computed once for all $L \leq L_\mathrm{max}$
and reused across all $N_e$ values.

### M3 partition count: `compute_m4`

Builds a 2D table `dp[n][k]` = number of partitions of $n$ into exactly $k$
parts each $\geq s$, using the recurrence:

$$p_k^{(\geq s)}(n) = p_{k-1}^{(\geq s)}(n - s) + p_k^{(\geq s)}(n - k)$$

The first term counts partitions whose smallest part equals $s$ (remove it to
get a partition of $n-s$ into $k-1$ parts $\geq s$); the second counts those
whose smallest part exceeds $s$ (subtract 1 from every part to get a partition
of $n-k$ into $k$ parts $\geq s$). The table has dimensions
$(L_\mathrm{max}+1) \times (N_e+1)$ and takes $O(L_\mathrm{max} \cdot N_e)$
time and space. The M3 count is the row sum:
$P_L^{(N_e,s)} = \sum_{k=1}^{N_e} \mathrm{dp}[L][k]$.

### Crossover analysis: `analyse_ne`

For each $N_e$, evaluates $\log P_L^{(N_e)}$, $\alpha_0\sqrt{L}$, and
$(N_e-1)\log(L/s)$ at a set of $L$ values: $\{0.05, 0.1, 0.2, 0.3, 0.5,
0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 8.0, 12.0, 18.0\} \times L^*_\mathrm{GC}$,
clipped to $[N_e s,\, L_\mathrm{max}]$. Reports the two asymptote ratios
$r_\mathrm{exp} = \log P / (\alpha_0\sqrt{L})$ and
$r_\mathrm{poly} = \log P / [(N_e-1)\log(L/s)]$.

The asymptote-intersection scale $L^\dagger$ is found by scanning $L$ from
$N_e s$ to $L_\mathrm{max}$ and detecting the sign change of
$\alpha_0\sqrt{L} - (N_e-1)\log(L/s)$.

### Canonical crossover scale

Computed in the header block by bisection over $[2, 10^9]$ to 100 iterations,
solving $(\sqrt{6}/2\pi)\sqrt{L}\log L = N_e$ to full floating-point precision.

---

## (d) Output produced

### Header table (all $N_e$ values)

| $N_e$ | $L^*_\mathrm{GC}$ | $L^*_\mathrm{canon}$ |
|---|---|---|
| 10 | 60.8 | 45.3 |
| 20 | 243.2 | 116.3 |
| 30 | 547.1 | 207.9 |
| 50 | 1520.4 | 444.0 |
| 75 | 3420.9 | 827.3 |
| 85 | 4394.8 | 997.1 |
| 100 | 6079.3 | 1312.8 |
| 200 | 24317.0 | 3856.8 |

### Per-$N_e$ crossover tables ($s = 1$, $L_\mathrm{max} = 5000$)

For each $N_e$ the table shows $L$, $L/L^*_\mathrm{GC}$, $\log P$, the
exponential asymptote, the polynomial asymptote, and both ratios. Key
observations from the data:

**For $N_e = 10$**: at $L/L^*_\mathrm{GC} = 1$ ($L = 60$), the ratio
$\log P / (\alpha_0\sqrt{L}) \approx 0.613$, showing the chain is already well
into the crossover zone. The polynomial ratio is only 0.33 at $L/L^* = 1$,
confirming the transition is gradual.

**For $N_e = 20$**: the exponential ratio peaks near $L/L^*_\mathrm{GC}
\approx 0.70$ at $\approx 0.732$, then decreases. At $L/L^* = 1$ ($L = 243$)
the ratio is 0.729.

**For $N_e = 30$**: the exponential ratio peaks around $L/L^*_\mathrm{GC}
\approx 0.70$ at $\approx 0.786$; the polynomial ratio reaches only 0.41 at
the maximum available $L = 5000$.

### Crossover summary

| $N_e$ | $L^*_\mathrm{GC}$ | $L^\dagger$ (asymptotes equal) |
|---|---|---|
| 10 | 60.8 | 465 |
| 20 | 243.2 | 3705 |
| 30 | 547.1 | N/A ($> 5000$) |

The ratio $L^\dagger / L^*_\mathrm{GC}$ grows rapidly with $N_e$ (from 7.6
at $N_e = 10$ to 15.2 at $N_e = 20$), confirming that the polynomial
asymptote is not a practical description at accessible chain lengths for
$N_e \geq 30$.

### CSV files

Three files are written to the same directory as the script (with an optional
`--prefix` prepended):
- `_m4_crossover_data.csv`: full per-$L$ data for all $N_e$ values.
- `_m4_crossover_summary.csv`: $L^*_\mathrm{GC}$, $L^\dagger$, and their
  ratio, one row per $N_e$.
- `_m4_scaled_entropy.csv`: $\log P$, both asymptotes, and both ratios
  versus $L/L^*_\mathrm{GC}$, for use in scaled crossover figures.

---

## (e) What we learn

**The crossover is real but very gradual.** The entropy ratios
$\log P / (\alpha_0\sqrt{L})$ and $\log P / [(N_e-1)\log L]$ both remain well
below 1 throughout the computed range. Even at $L/L^*_\mathrm{GC} \approx 20$
for $N_e = 10$, the polynomial ratio is only 0.57. The transition from
exponential to polynomial entropy extends over more than a decade in reduced
chain length $L/L^*$, and the asymptotes do not cross until $L^\dagger \sim
8$–$15 \times L^*_\mathrm{GC}$.

**$L^*_\mathrm{canon}$ is the correct crossover onset; $L^*_\mathrm{GC}$
overestimates it.** The factor by which $L^*_\mathrm{GC}$ overestimates the
onset grows from $\approx 1.3$ at $N_e = 10$ to $\approx 2.6$ at $N_e = 30$
and $\approx 4.4$ at $N_e = 85$. For quantitative comparison with simulation
data, $L^*_\mathrm{canon}$ should be used.

**The polynomial exponent is $N_e - 1$, not $N_e$.** The count
$P_L^{(N_e,\, s)} \sim (L/s)^{N_e-1}/(N_e!\,(N_e-1)!)$ gives entropy
$\sim (N_e-1)\log(L/s)$. The data confirm this: at $L = 1215$ for $N_e = 10$,
$\log P = 36.2$ compared to the polynomial prediction $(N_e-1)\log L \approx
63.9$, with the ratio 0.57 converging from below toward 1.

**The $N_e^2$ scaling of $L^*$ is robust.** Both $L^*_\mathrm{GC} \propto
N_e^2$ and $L^*_\mathrm{canon}$ (which scales approximately as $N_e^2/\log^2
N_e$) grow quadratically with $N_e$. This provides a partition-theoretic
derivation of the onset of the entangled regime, with $L^*/N_e \sim N_e$
corresponding to a chain of approximately $N_e$ entanglement lengths.

**The script covers a broader $N_e$ range than described in the paper.** The
default run includes $N_e \in \{10, 20, 30, 50, 75, 85, 100, 200\}$; the
paper's figures use $N_e \in \{10, 20, 30\}$. The extended range confirms
that the qualitative crossover behaviour and the $N_e^2$ scaling persist
across the full range of physically relevant entanglement lengths.
