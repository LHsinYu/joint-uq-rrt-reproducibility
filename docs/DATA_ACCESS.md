# Data access

The empirical analysis uses the 2012 Taiwan Social Change Survey (TSCS) gender module. These microdata are provided by a third party and may be subject to registration, citation, and redistribution conditions. This repository therefore does not contain the source `.sav` or `.csv` files or any respondent-level derivatives.

To reproduce the empirical analysis:

1. Obtain the TSCS 2012 data from the official survey archive under its current access terms.
2. Keep the downloaded file outside version control.
3. Set the `TSCS_DATA` environment variable to the local `.sav` file.
4. Run `R/TSCS2012_empirical_recoding_sensitivity.R` from the repository root.

Example:

```sh
TSCS_DATA=/absolute/path/to/tscs2012q2.sav \
  Rscript R/TSCS2012_empirical_recoding_sensitivity.R
```

The script checks the required variable names before analysis. The variables used are `id`, `a1`, `a10`, `b1`, `j1`, `j2`, `j5`, `k9`, and `wr`. Users should consult the official TSCS codebook for definitions and permitted uses.

The simulation programs do not require the TSCS data and can be run independently.
