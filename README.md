# Reproducibility Materials for Joint UQ-RRT Modeling

This repository contains the estimation, testing, simulation, and diagnostic code supporting the manuscript **A Randomized Response Framework for Joint Modeling of Latent Sensitive Behavior and Ordinal Attitudinal Responses**.

## Contents

- `R/Greenberg_EM_reviewer_sensitivity.R`: Greenberg UQ-RRT data generation, EM estimation, reduced and heterogeneous models, joint Wald and likelihood-ratio tests, numerical-Hessian standard errors, multiple-start diagnostics, checkpointing, and simulation summaries.
- `R/QQ_nonproportional_odds_simulation_formal.R`: proportional-odds misspecification study.
- `R/Read_TSCS2012_revised.R`: TSCS 2012 data preparation used by the empirical analysis.
- `R/TSCS2012_empirical_recoding_sensitivity.R`: empirical ordinal-recoding sensitivity analysis.
- `scripts/summarize_reviewer_simulations.py`: aggregation of scenario- and parameter-level simulation outputs.
- `docs/ANALYSIS_DESIGN.md`: definitions of the reviewer-requested analyses.
- `docs/DATA_ACCESS.md`: data availability and access instructions.
- `results/small_sample_diagnostics.csv`: the 12 setting-level diagnostics reported in Web Appendix Table S6.
- `results/software_session_info.txt`: tested software environment.

## Software

The release was tested with R 4.6.1 on macOS (arm64). The R scripts use base R and the recommended packages `MASS` and `foreign`; version details are recorded in `results/software_session_info.txt`. The aggregation utility uses Python 3.9 or later and only the Python standard library.

## Quick validation

Run these commands from the repository root. They create output files under ignored local result directories.

```sh
EM_SIM_KT=2 EM_N_CORES=1 EM_N_STARTS=2 EM_RUN_ID=quick_test \
  Rscript R/Greenberg_EM_reviewer_sensitivity.R

QQ_RUN=1 QQ_QUICK_TEST=1 \
  Rscript R/QQ_nonproportional_odds_simulation_formal.R
```

Formal simulations use 1,000 replications and are computationally intensive. Every run writes its settings, seed, effective denominators, convergence diagnostics, standard-error failures, extreme-solution flags, and multiple-start comparisons to the output directory.

## Example formal settings

The small-sample C1 low-prevalence setting at `n=200` can be run with:

```sh
EM_CASE=C1 EM_PREVALENCE=low EM_SIM_N=200 EM_SIM_KT=1000 \
  EM_SIM_SEED=20260916 EM_N_STARTS=10 EM_RUN_ID=C1_low_n200 \
  Rscript R/Greenberg_EM_reviewer_sensitivity.R
```

Change `EM_CASE` to `C6`, `EM_PREVALENCE` to `high`, and `EM_SIM_N` to `300` or `400` for the other Table S6 settings. Failed replications are retained and are not replaced.

## Data availability

The Taiwan Social Change Survey microdata are third-party data and are not redistributed in this repository. Authorized users can set `TSCS_DATA` to the downloaded `.sav` file and run:

```sh
TSCS_DATA=/path/to/tscs2012q2.sav \
  Rscript R/TSCS2012_empirical_recoding_sensitivity.R
```

See `docs/DATA_ACCESS.md` for details. No raw or respondent-level TSCS records are included.

## Output interpretation

Point-estimation summaries use finite converged estimates. Standard-error and coverage summaries use replications with valid observed-information covariance matrices. Extreme solutions are flagged by the prespecified rule

`max(abs(estimate)) > 10` or `max(SE) > 20`

and are retained in the primary summaries. Conditional rejection or coverage rates must be read together with their effective denominators and undefined-result rates.

## License

The code is released under the MIT License. The TSCS data remain subject to the terms imposed by their data provider and are not covered by this repository's license.
