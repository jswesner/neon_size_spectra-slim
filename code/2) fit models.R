library(brms)
library(isdbayes)
library(tidyverse)

# Instead of re-fitting, you can load the fitted models below

# load models (these are already fit)
fit_temp = readRDS("models/fit_temp.rds")
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

# Re-fit models -----------------------------------------------------------
# !!!! Each of the models below takes up to 24 hours to run. They were run on a cluster. Here, we've added them
# for clarity, but with only a small number of iterations to prevent overwhelming memory. If you fit the models below
# they should not be used for inference (without increasing the iterations)

#1) load data (NEON body size data)
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")

#2) fit main model
# fit_temp_om_gpp = brm(dw | vreal(no_m2, xmin, xmax) ~ log_om_s*mat_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id),
#                       data = dat_2022_clauset,
#                       stanvars = stanvars,    # required for truncated Pareto via isdbayes package
#                       family = paretocounts(),# required for truncated Pareto via isdbayes package
#                       prior = c(prior(normal(-2, 0.5), class = "Intercept"),
#                                 prior(normal(0, 0.2), class = "b"),
#                                 prior(exponential(7), class = "sd")),
#                       iter = 2000,
#                       cores = 4,
#                       threads = 12, 
#                       chains = 4)
# 
# saveRDS(fit_temp_om_gpp, file = "models/fit_temp_om_gpp.rds")

#2) fit submodel using the update function

# fit_temp = update(fit_temp_om_gpp, formula = . ~ mat_s + (1 | sample_id) + (1 | year) + (1 | site_id))
# saveRDS(fit_temp, file = "models/fit_temp.rds")
