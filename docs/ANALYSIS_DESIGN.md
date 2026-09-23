# Reviewer requested simulation design

## Existing baseline results

The previously completed baseline contains C1-C6, low and high prevalence, sample sizes 500 and 1000, and 1000 fixed Monte Carlo replications per scenario at the article RRT setting. Those results are not rerun here.

## Non proportional odds robustness

- Data-generating scenarios: proportional odds, mild non-proportional odds, and moderate non-proportional odds.
- Sample sizes: 500 and 1000.
- Replications: 1000 accepted fits per scenario, with all failed attempts retained.
- Fitted model: the article proportional-odds EM model.
- Outputs: convergence, Hessian validity, extreme solutions, prevalence error, posterior-status error, fitted ordinal-probability error, parameter summaries where a common truth exists, and joint Wald/LR rejection rates with Monte Carlo standard errors.
- Interpretation: because both latent groups use the same NPO data-generating distribution, equality across latent groups is a true null. Wald/LR rejection rates therefore assess size robustness under proportional-odds misspecification.

## Known RRT design sensitivity

The fitted RRT mechanism equals the data-generating mechanism. New settings are:

- weaker information: `p=0.30, c=0.25`;
- stronger information: `p=0.70, c=0.25`;
- different unrelated-question prevalence: `p=0.50, c=0.50`.

For each setting, four anchors are used: C1 low prevalence with `n=500`, C1 high prevalence with `n=1000`, C6 low prevalence with `n=500`, and C6 high prevalence with `n=1000`. C1 provides type-I-error evidence; C6 provides power and heterogeneous-model estimation evidence.

## RRT mechanism misspecification

Data are generated at `p=0.50, c=0.25`. The fitted values are separately changed to:

- `p=0.45, c=0.25`;
- `p=0.55, c=0.25`;
- `p=0.50, c=0.20`;
- `p=0.50, c=0.30`.

C1 and C6 are run under low prevalence at sample sizes 500 and 1000. These scenarios quantify bias, RMSE, coverage, posterior effects, test rejection behavior, convergence, and extreme solutions when the assumed device probabilities differ moderately from their true values.

## Smaller sample performance

C1 and C6 are run at sample sizes 200, 300, and 400 under both low and high prevalence. Results are interpreted as scenario-specific operating characteristics and not as a universal minimum-sample-size theorem.

## Replication and numerical controls

- Formal RRT sensitivity scenarios use 1000 fixed replications and five starting values by default.
- Formal NPO scenarios use 1000 accepted fits and five starting values; failed attempts remain in the replication-level output.
- Outputs retain point-estimation failures, invalid standard errors, Hessian diagnostics, extreme estimates, start agreement, and undefined test statistics.
- Type-I error and power must use the same joint null for Wald and LR comparisons.

## Empirical recoding analysis dependency

The available source calls `Read_Zhang2018()` but does not contain that function or the underlying empirical data. Original-category frequencies, alternative recoding, and empirical ordinal-versus-collapsed comparisons cannot be computed until the data or a runnable loader is supplied. No synthetic substitute should be reported as empirical evidence.
