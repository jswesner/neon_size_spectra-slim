
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Temperature and resource supply drive continental variation in size structure of freshwater food webs

This page provides data and code for *Gjoni et al. Size spectra in
freshwater streams are consistent across temperature and resource
supply*.

All figures and tables in the manuscript can be recreated by running the
R scripts below. The scripts are named in order (e.g., 1, 2, 3…), so
that script 2 will not work without running script 1 first, and so
forth.

The scripts below are in the folder `code`:

**1) choose priors.R**: Conducts prior predictive simulation.

**1b) summarize data.R**: Calculates basic summary statistics of the
data. e.g., total number of body sizes, number of taxa, etc.

**2) fit models.R**: Fits truncated Pareto GLMs used in the paper.
(NOTE: These models take ~10-20 hours to run. We suggest using a cluster
to fit them, as we did).

**3.1-3) check_models….R**: Posterior predictive checks of the models
from step 2.

**4-12)…**: Creates figures in the manuscript and SI.

**13-14)…**: Refits models with only fish or only inverts.

**15) interaction plot fish and macros.R**: Creates figures for SI
interaction plot

**sampler_types.R**: Refits models with only inverts collected from
Surber samplers

**isd_plots.R**: Creates ISD plots for each site, averaged over
collections

## Raw data processing

The raw data were obtained from the \[National Ecological Observatory
Network (NEON)\] (<https://www.neonscience.org/>) using the
`neonUtilities` package. Code to download and wrangle the raw data are
in the folder `code/get_neon_data`.

The scripts below are in the folder `code/get_neon_data`:

**1) get_neon_macro_data.R**: Downloads NEON data for macroinvertebrates
using the `neonUtilities` package. Then wrangles data to get individual
body masses using length-mass regression.

**1b) get_neon_macro_data_surber_only.R**: Same as above, but limits
data to only samples from Surber samplers. This is for testing how
sampling method affects the size spectrum.

**2) get_neon_fish_data.R**: Downloads NEON data for fish using the
`neonUtilities` package. Then wrangles data to get individual body
masses and fish abundance data from NEON’s three-pass removal sampling.

**3) fit_three_pass_fish.R**: Fits three-pass removal models to estimate
fish densities.

**4) get_fish_size_density.R**: Combines NEON data on fish sizes with
results of three-pass removal estimates. The result is a data frame with
fish sizes (mgDM) and the associated density of those sizes (no_m2).
This places the fish data in the same format as the macroinvertebrate
data so they can be combined.

**5) combine_fish_macros.R**: Combines the macroinvertebrate and fish
body size data sets.

**5b) combine_fish_macros_samplertype.R**: Same as above but limits
macroinvertebrate data to only samples from Surber samplers.

**6) cull_xmin…**: Estimates $x_{min}$, the minimum body size for which
the data follow a power law. Then culls the data below that size. This
produces the primary data that is used to fit ISD models.

## Packages

    #>          package     version
    #> 1      bayesplot 1.15.0.9000
    #> 2           brms      2.23.1
    #> 3        cowplot       1.2.0
    #> 4          ggh4x       0.3.1
    #> 5        ggimage       0.3.5
    #> 6        ggplot2       4.0.3
    #> 7       ggridges       0.5.7
    #> 8         ggtext       0.1.2
    #> 9       ggthemes       5.2.0
    #> 10           hms       1.1.4
    #> 11      isdbayes       0.1.0
    #> 12       janitor       2.2.1
    #> 13         knitr        1.51
    #> 14     lubridate       1.9.5
    #> 15 neonUtilities       3.0.3
    #> 16     patchwork       1.3.2
    #> 17      poweRlaw       1.0.0
    #> 18        readxl       1.4.5
    #> 19          renv       1.2.3
    #> 20     rmarkdown        2.31
    #> 21 rnaturalearth       1.2.0
    #> 22         rstan      2.32.7
    #> 23          rsvg       2.7.0
    #> 24     tidybayes       3.0.7
    #> 25     tidyverse       2.0.0
    #> 26       viridis       0.6.5
