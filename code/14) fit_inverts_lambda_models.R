library(rstan)
library(tidyverse)

rstan_options(autowrite = TRUE)
rstan_options(threads_per_chain = 1)

# load data
neon_sizes_2016_2021 = readRDS(file = "data/derived_data/macro_dw-wrangled.rds") %>% 
  filter(year >= 2016 & year <= 2021)

# compile model
stan_spectra_mod_gpp_x_temp = stan_model("models/stan_spectra_mod_gpp_x_temp.stan")


# make data and fit model ---------------------------------------------------------

dat = neon_sizes_2016_2021

stan_data_interaction = list(N = nrow(dat),
                             mat_s = dat$mat_s,
                             gpp_s = dat$log_gpp_s,
                             year = dat$year_int,
                             site = dat$site_int,
                             sample = dat$sample_int,
                             n_years = length(unique(dat$year_int)),
                             n_sites = length(unique(dat$site_int)),
                             n_samples = length(unique(dat$sample_int)),
                             counts = dat$no_m2,
                             x = dat$dw,
                             xmin = dat$xmin,
                             xmax = dat$xmax)

fit_interaction = sampling(object = stan_spectra_mod_gpp_x_temp,
                           data = stan_data_interaction,
                           iter = 2000, chains = 4, cores = 4)

saveRDS(fit_interaction, file = paste0("models/stan_fishonly_gppxtemp",Sys.Date(),".rds"))










