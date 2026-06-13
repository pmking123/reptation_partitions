# Appendix A Convergence Diagnostics

## Summary of `appendix_diagnostics.py`

Produces the acceptance rate and burn-in convergence diagnostics for the
uniform-measure sampler cited in Appendix A §A.2–A.3 of the paper.

**Dependencies:** `mc_uniform_fast.py`, `partition_numpy.py`. No external
dependencies beyond numpy.

---

## What it produces

### (1) Acceptance rates

Runs 200,000 steps after a short equilibration at each production $L$ and
reports the fraction of accepted proposals. Expected range 0.89–0.96 (rising
with $L$), as quoted in Appendix A §A.2.

### (2) Burn-in diagnostics

Run at each $L$ via `burnin_diagnostics(L)`, with burn-in and post-burn-in
lengths both equal to $1000L$ steps and 4 independent chains:

- **Running mean**: step at which the running mean of $k$ first enters the
  $\pm 1\%$ band around the second-half mean, expressed as a fraction of the
  burn-in length.

- **Gelman–Rubin $\hat{R}$**: computed across 4 independent chains (all
  starting from the default `init_arrays` state, different seeds). Values
  below 1.01 pass. The paper cites $\hat{R} = 1.007$ at $L = 500$.

- **Starting-point sensitivity**: runs the same burn-in from three starting
  states — under-dispersed ($k = 1$, single part), default
  ($k \approx 1.28\sqrt{L}$), and over-dispersed ($k = L$, all unit parts)
  — and reports the spread in second-half mean $k$. The paper cites a spread
  below $0.5\%$.

---

## Note on missing $L = 100$

`main()` currently calls `burnin_diagnostics` only for $L \in \{200, 500\}$.
The $L = 100$ result cited in Appendix A §A.3 requires either adding
`burnin_diagnostics(100)` to `main()` or running it manually. This should
be fixed before the repository is finalised.
