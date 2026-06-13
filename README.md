# Reptation Partitions

Computational supplement to

> **Reptating Polymers and Integer Partitions: a Bose-Einstein Correspondence**  
> *Journal of Physics A* (submitted 2026)

This repository contains all Python scripts, SLURM job scripts, and human-readable summaries needed to reproduce the numerical results, figures, and verification claims in the paper. Every quantitative claim in the manuscript is mapped to its generating script and output log in `MANIFEST.md`.

---

## Repository structure

```
.
├── python_scripts/       # Analysis, simulation, and verification code
├── job_scripts/          # SLURM shell scripts for HPC cluster execution
├── summaries/            # Markdown (and PDF) documentation for each script
├── MANIFEST.md           # Audit trail: paper claim → script → output value
└── README.md
```

---

## Paper–script map

The table below maps each section of the paper to the script(s) that produce its numerical content. The `summaries/` directory contains a detailed documentation file for each script.

| Paper section | Script(s) | What it produces |
|---|---|---|
| §5.2 Rogers–Ramanujan entropy | `rr_entropy.py` | Entropy coefficients $\alpha_\mathrm{RR}$, $\alpha_d$; convergence tables; CSV data for Figs. 3–4 |
| §6.1–6.4 Canonical vs GC mean | `be_saddle.py` | Exact $p_k(L)$ via DP; $k^*$, $E[k]_\mathrm{canon}$, $E[k]_\mathrm{GC}$; scaling fits; Table 2 |
| §7.3 M3 crossover | `crossover.py` | $P_L^{(N_e)}$ via DP; crossover tables and CSV for all $N_e$ values |
| §8.2 Decoupling theorem | `partition_simulations.py` | Exhaustive verification over $2^{20}-1$ compositions; $w(\lambda)$ exact enumeration |
| §8.2 Proposition 7 | `fixed_point_verification.py` | Gap-2 fixed-point verification over $2^{25}-1$ compositions |
| §8.3 BE run-length distribution | `mc_uniform_fast.py`, `rld_fast.py` | Uniform-measure MC; KL divergence jackknife table; $\chi^2$ table |
| §8.4 Prediction 2 feasibility | `prediction2.py` | Lattice reptation simulation; $\log\tau$ vs $L$ scaling fits |
| Appendix A.2–A.3 | `appendix_diagnostics.py` | Acceptance rates; Gelman–Rubin $\hat{R}$; starting-point sensitivity (uniform sampler) |
| Appendix A.2 | `burnin_diagnostics.py` | Trace plots; running-mean convergence; $\hat{R}$; starting-point sensitivity (weighted sampler) |

**Shared infrastructure** (no direct paper section, used by scripts above):

| Script | Role |
|---|---|
| `partition_numpy.py` | Numpy-accelerated partition representation and MH move operations; shared by `mc_fast.py` and `mc_uniform_fast.py` |
| `mc_fast.py` | Weighted-measure $\tilde{w}$ MH sampler; validates canonical/GC inequivalence |
| `aggregate.py` | Collects JSON output from SLURM job-array tasks and prints consolidated tables |
| `core.py` | Generates all partitions of $L$ as non-increasing tuples; used by `mc_fast.py` validation |

---

## Dependencies

All scripts require **Python 3.8+** and **NumPy**. The simulation and diagnostic scripts additionally require **Matplotlib** (`burnin_diagnostics.py`). No other external dependencies are needed.

The SLURM scripts in `job_scripts/` are configured for an AMD EPYC 9454 cluster (48 cores, 256 GB RAM) and honour `SLURM_CPUS_PER_TASK` for parallelism. They can be adapted to other HPC environments by adjusting the `#SBATCH` headers.

---

## Running the scripts

### Standalone (sequential) mode

Every script can be run directly without SLURM:

```bash
cd python_scripts
python3 be_saddle.py                      # §6 canonical/GC comparison (L up to 10000; ~30 min)
python3 rr_entropy.py                     # §5.2 entropy coefficients (~5 min)
python3 fixed_point_verification.py       # Prop. 7, L≤25 (~2 hours)
python3 partition_simulations.py          # §8.2 decoupling theorem (~10 min)
python3 crossover.py                      # §7.3 M3 crossover (~10 min)
python3 mc_uniform_fast.py                # §8.3 uniform sampler, quick run
python3 rld_fast.py                       # §8.3 KL divergences (~11 hours on 48 cores)
python3 prediction2.py                    # §8.4 Prediction 2 feasibility (~several hours)
python3 appendix_diagnostics.py           # Appendix A diagnostics (~10 min)
```

For scripts with long runtimes (`rld_fast.py`, `fixed_point_verification.py` at $L \leq 25$, `prediction2.py`), the SLURM job scripts are strongly recommended.

### SLURM job-array mode

Each script that supports parallel execution has a matching job script:

```bash
cd job_scripts
sbatch rld.sh                   # §8.3 KL table (rld_fast.py, parallel)
sbatch mc_uniform.sh            # §8.3 mode/RLD tasks (mc_uniform_fast.py, parallel)
sbatch mc.sh                    # §6 weighted sampler (mc_fast.py, parallel)
sbatch partition_simulations.sh # §8.2 decoupling (partition_simulations.py)
sbatch fixed_point_verification.sh  # Prop. 7
sbatch be_saddle.sh             # §6 DP tables
sbatch crossover.sh             # §7.3 M3 crossover
sbatch rr_entropy.sh            # §5.2 entropy coefficients
sbatch prediction2.sh           # §8.4 Prediction 2 (6 configurations)
sbatch mc_burnin_diagnostics.sh # Appendix A (weighted sampler)
```

After SLURM jobs complete, combine outputs with `aggregate.py`:

```bash
python3 aggregate.py --mode uniform --indir results_mc_uniform
python3 aggregate.py --mode weighted --indir results_mc
```

---

## Summaries

The `summaries/` directory contains a documentation file for each script, structured in five sections: (a) theory, (b) what the script does, (c) how it does it, (d) output produced, (e) what we learn. Both Markdown and PDF versions are provided.

| Summary | Script documented |
|---|---|
| `be_saddle_summary` | `be_saddle.py` |
| `rr_entropy_summary` | `rr_entropy.py` |
| `partition_simulations_summary` | `partition_simulations.py` |
| `fixed_point_verification_summary` | `fixed_point_verification.py` |
| `mc_fast_summary` | `mc_fast.py` |
| `mc_uniform_fast_summary` | `mc_uniform_fast.py` |
| `partition_numpy_summary` | `partition_numpy.py` |
| `rld_fast_summary` | `rld_fast.py` |
| `crossover_summary` | `crossover.py` |
| `prediction2_summary` | `prediction2.py` |
| `aggregate_summary` | `aggregate.py` |
| `appendix_diagnostics_summary` | `appendix_diagnostics.py` |

`burnin_diagnostics.py` and `core.py` do not have dedicated summaries. `burnin_diagnostics.py` produces burn-in trace plots, running-mean convergence diagnostics, Gelman–Rubin $\hat{R}$, and starting-point sensitivity checks for the weighted sampler (`mc_fast.py`); its results are not directly cited in the paper. `core.py` is a 15-line utility that generates all partitions of $L$; it is used only by `mc_fast.py`'s small-$L$ validation routine.

---

## Verification

`MANIFEST.md` provides a full audit trail mapping every quantitative claim in the paper to its generating script, output log, paper value, and recomputed value, with a pass/fail status. All deterministic claims (exact enumeration, dynamic programming) have status **OK**. The two stochastic claims with refitted values (Prediction 2 fit coefficients) are also **OK** after the refit applied prior to submission.

---

## License

MIT License. Copyright (c) 2026 Paul King.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. The Software is provided "as is", without warranty of any kind.
