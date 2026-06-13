# Verification manifest — *Reptating Polymers and Integer Partitions*

Audit trail mapping each computational claim in the manuscript to the script
that produces it, the archived log, the value quoted in the paper, the value
recomputed from the log, and a status. Prepared to accompany the referee
response and to seed the reproducibility appendix.

Status key: **OK** = paper matches log exactly; **FIX** = correction required,
exact edit identified; **NOTE** = matches, minor housekeeping only.

---

## 1. Exact enumeration (deterministic, arbitrary-precision)

| # | Claim (paper location) | Script | Log | Paper value | Recomputed | Status |
|---|---|---|---|---|---|---|
| 1.1 | Decoupling factorisation over all compositions of L≤20 | `partition_simulations.py` | `partition_simulations_out_3.txt` | 1,048,575 compositions, 0 failures | 1,048,575; 0 failures | **OK** |
| 1.2 | Lemma (geometric = injective projections), all signed realisations L≤20 | `partition_simulations.py` | `partition_simulations_out_3.txt` | 3,486,784,400 realisations, 0 mismatches | 3,486,784,400; 0 mismatches | **OK** |
| 1.3 | Worked example w_strong(3,2,3,2,3) | `partition_simulations.py` (extremal block) | `partition_simulations_out_3.txt` | w_H=2, w_V=2, w_strong=4 | 2·2=4 confirmed by enumeration | **OK** (corrects old "8" in TeX; see edit set) |
| 1.4 | Remark weak_sa counts | direct enumeration | (this session) | w_weak=28, w_strong=4 | 28; 4 (H-counts 8,6,6,8) | **OK** (corrects old "8") |
| 1.5 | Prop. 7 attractor non-uniqueness, L≤20 | `fixed_point_verification.py` | `fpv_out_1.txt` | 679,104 divergent / 1,048,575; 0 gap-2 failures | 679,104; PASS all 4 tests | **OK** |
| 1.6 | Prop. 7 zero-failure confirmation, L≤25 | `fixed_point_verification.py` | `fpv_out_1.txt` | 33,554,431 compositions, 0 failures | 33,554,431; PASS (25,420,443 divergent) | **OK** |
| 1.7 | Entropy coefficients α0, α_RR, α_d | `rr_entropy.py` | `rre_out_1.txt`, `entropy_coefficients.csv` | 2.5651, 1.62231, 1.81380; α_RR/α0=√(2/5)=0.6325 | identical | **OK** |
| 1.8 | DP spot checks p(L), d(L) | `rr_entropy.py` | `rre_out_1.txt` | — | p(20)=627, p(50)=204226, d(50)=3658 | **OK** (also fixes provenance of Table kstar / canonical means) |
| 1.9 | O(log L/√L) convergence of α | `rr_entropy.py` | `rre_out_1.txt` | fitted k≈0.95 (RR,d), 1.26 (unrestricted) | identical | **OK** |

## 2. Monte Carlo — uniform partition measure (stochastic)

| # | Claim (paper location) | Script | Log | Paper value | Recomputed | Status |
|---|---|---|---|---|---|---|
| 2.1 | Sampler validation vs exact enumeration | `mc_uniform_fast.py` `validate()` | (validate output) | "L≤20"; TV<0.04; E[k] within 2% | validate runs L=4..20 | **FIX** parenthetical: p(16)=231 → p(20)=627 |
| 2.2 | KL table, BE vs geom, 12-chain jackknife | `rld_fast.py` (Step 2c) | `rld_out_1.txt` (rld_updated) | see table below | matches to quoted precision | **OK** |
| 2.3 | k* modal estimate at L=500 | `mc_uniform_fast.py` `run_mode_task` | `mode_500.json` | modal 53; σ_k≈20 | modal_k=53; ksig=19.81 | **OK** |
| 2.4 | k* / MC sampling extent L≤2000 | `mc_uniform_fast.py` | `mode_2000.json` | "L≤2000" | modal_k=128 (18 h run) | **OK** |
| 2.5 | BE-convergence monotone in L (burn-in adequacy at L=500) | `rld_fast.py` (Step 3) | `rld_out_1.txt` | "visible for L≥100" | KL_BE 7.7e-4→2.8e-4→1.1e-4→2e-5 (L=100/200/500/2000) | **OK** |
| 2.6 | Weighted-sampler burn-in convergence at L=500 | `burnin_diagnostics.py` | `burnin_convergence_L500.txt` | (scaling-study prerequisite) | running-mean PASS; R̂=1.0067; start spread 0.49% | **OK** |

### KL table cross-check (claim 2.2)

| L | Paper KL_BE | Log Step 2c | Paper KL_geom | Log Step 2c | Paper ratio | Log Step 2c |
|---|---|---|---|---|---|---|
| 100 | (1.0±0.3)e-3 | 0.0010±0.0003 | (8.8±0.4)e-2 | 0.0881±0.0035 | 87±27 | 87.0±27.4 |
| 200 | (3±1)e-4 | 0.0003±0.0001 | (1.14±0.04)e-1 | 0.1142±0.0044 | 410±60 | 414.5±62.4 |
| 500 | (1±1)e-4 | 0.0001±0.0001 | (1.66±0.07)e-1 | 0.1657±0.0072 | 1200±600 | 1202.9±585.7 |

## 3. Lattice reptation — Prediction 2 feasibility (stochastic)

| # | Claim (paper location) | Script | Log | Paper value | Recomputed | Status |
|---|---|---|---|---|---|---|
| 3.1 | Power-law fit exponent c (L=100–500) | `prediction2.py` + fit | `primary.csv` | c=1.87, R²=0.998 | **2.13, R²=0.999** (paper's 1.87 = full L≤1000 range, contradicts stated exclusion) | **FIX** |
| 3.2 | √L fit slope a (L=100–500) | `prediction2.py` + fit | `primary.csv` | a=0.191, a/α0=0.074, R²=0.918 | **a=0.273, a/α0=0.107, R²=0.983** | **FIX** |
| 3.3 | Power law beats √L at accessible L | — | `primary.csv` | "decisively better" | holds (R² 0.999 vs 0.983) | **OK** |
| 3.4 | L=750,1000 undersampling exclusion | `prediction2.py` | `primary.csv`, logs | T_obs/τ ≳60 (L≤500), <30 (L≥750) | exclusion now consistently applied after refit | **NOTE** verify τ thresholds if τ available |
| 3.5 | Disorder-sensitivity check (independent obstacle fields) | `prediction2.py` | `seed_check_1.csv`, `seed_check_2.csv` | (separate realisations) | files present | **NOTE** cross-check if cited quantitatively |

## 4. M4 crossover

| # | Claim (paper location) | Script | Log | Paper value | Recomputed | Status |
|---|---|---|---|---|---|---|
| 4.1 | M4 crossover L* ~ Ne² | `crossover.py` | `s1_m4_crossover_summary.csv` (+s2,s3) | crossover form | L*_obs/L*_theory rises 7.6→15.2→20.3 (Ne=10,20,30) | **NOTE** confirm paper's claim matches the observed pre-asymptotic ratio trend |

---

## Outstanding edits (this round)

**Typos (factual):**
- Conclusion: "validated ... for $L \le 16$" → "$L \le 20$".
- Section 8 validation para: "($p(16)=231$ partitions)" → "($p(20)=627$ partitions)".

**Prediction 2 refit (Section 8, tab:pred2 paragraph):**
- c: 1.87 → 2.13 (R² 0.998 → 0.999)
- a: 0.191 → 0.27, a/α0: 0.074 → 0.11 (R² 0.918 → 0.983)
- Add: "unweighted least squares over L=100–500; excluded L=750,1000 discussed below."

**Decoupling section (separate edit set, 10 edits + lemma):**
- Definitions of signed realisation + strong self-avoidance; equivalence Lemma;
  restated Theorem + proof; corrected Example (w_strong=4); Remark weak_sa (8→4);
  edge-vs-site Remark; verification subsection recast; table/discussion refs.

## Housekeeping (script self-consistency for supplementary release)
- `partition_simulations.py`: Part-2 banner still says "L=1..12"; docstring still
  says "4095 compositions". Update to describe combined check, 2^20−1 / 3^20−1.
- L=2000 canonical mean: quote exact DP 145.64 (not MC running mean 145.42) wherever
  the canonical mean appears.
- Appendix A.2 still needs per-run acceptance rate and τ_int columns if not present
  in per-job .out files.
