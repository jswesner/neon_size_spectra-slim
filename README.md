
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Size spectra in freshwater streams are consistent across temperature and resource supply

This page provides data and code for *Gjoni et al. Size spectra in
freshwater streams are consistent across temperature and resource
supply*.

All figures and tables in the manuscript can be recreated by running the
R scripts in the `code` folder. The scripts are named in order (e.g., 1,
2, 3…), so that script 2 will not work without running script 1 first,
and so forth.

The scripts provide the following:

**1) choose priors.R**: Conducts prior predictive simulation.

**2) fit models.R**: Fits truncated Pareto GLMs used in the paper.
(NOTE: These models take ~10-20 hours to run. We suggest using a cluster
to fit them, as we did).

**3.1-3) check_models….R**: Posterior predictive checks of the models
from step 2.

**4-12)…**: Create main figures in the manuscript.

**13-14)…**: Refits models with only fish or only inverts.

These scripts use data that were wrangled from raw NEON data to generate
values of body sizes and Gross Primary Production. Code for that can be
found in a separate repo at
<https://github.com/jswesner/neon_size_spectra>.
