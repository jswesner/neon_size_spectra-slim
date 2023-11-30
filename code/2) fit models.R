library(brms)
library(isdbayes)
library(tidyverse)



# Instead of re-fitting, you can load the fitted models below

# load models (these are already fit)
fit_temp = readRDS("models/fit_temp.rds")
fit_om = readRDS("models/fit_om.rds")
fit_gpp = readRDS("models/fit_gpp.rds")
fit_temp_om = readRDS("models/fit_temp_om.rds")
fit_temp_gpp = readRDS("models/fit_temp_gpp.rds")
fit_om_gpp = readRDS("models/fit_om_gpp.rds")
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")




# Re-fit models -----------------------------------------------------------
# !!!! Each of the models below takes up to 24 hours to run. They were run on a cluster. Here, we've added them
# for clarity, but with only a small number of iterations to prevent overwhelming memory. If you fit the models below
# they should not be used for inference (without increasing the iterations)

#1) load data (NEON body size data)
dat_all = readRDS("data/derived_data/dat_all.rds")

#2) fit main model
fit_temp_om_gpp = brm(dw | vreal(no_m2, xmin, xmax) ~ log_om_s*mat_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id),
                      data = dat_all,
                      stanvars = stanvars,    # required for truncated Pareto via isdbayes package
                      family = paretocounts(),# required for truncated Pareto via isdbayes package
                      prior = c(prior(normal(-1.5, 0.2), class = "Intercept"),
                                prior(normal(0, 0.1), class = "b"),
                                prior(exponential(7), class = "sd")),
                      iter = 10, 
                      chains = 1)

#2) fit submodels using the update function

fit_temp = update(fit_temp_om_gpp, formula = . ~ mat_s + (1 | sample_id) + (1 | year) + (1 | site_id))
fit_om = update(fit_temp_om_gpp, formula = . ~ log_om_s + (1 | sample_id) + (1 | year) + (1 | site_id))
fit_gpp = update(fit_temp_om_gpp, formula = . ~ log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id))
fit_temp_om = update(fit_temp_om_gpp, formula = . ~ mat_s*log_om_s + (1 | sample_id) + (1 | year) + (1 | site_id))
fit_temp_gpp = update(fit_temp_om_gpp, formula = . ~ mat_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id))
fit_om_gpp = update(fit_temp_om_gpp, formula = . ~ log_om_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id))

#3) add errors in variables
# load temperature data that has correct sd's (on the standardized scale)
temperature_mean_annual = readRDS(file = "data/derived_data/temperature_mean-annual.rds") %>% 
  rename(site_id = siteID) %>% 
  ungroup

dat_all_sd = dat_all %>% left_join(temperature_mean_annual %>% select(site_id, mat_s_mean, mat_s_sd))
saveRDS(dat_all_sd, file = "data/derived_data/dat_all_sd.rds")

fit_temp_sd = brm(dw | vreal(no_m2, xmin, xmax) ~ me(mat_s_mean, mat_s_sd) + (1 | sample_id) + (1 | year) + (1 | site_id),
                      data = dat_all_sd,
                      stanvars = stanvars,    # required for truncated Pareto via isdbayes package
                      family = paretocounts(),# required for truncated Pareto via isdbayes package
                      prior = c(prior(normal(-1.5, 0.2), class = "Intercept"),
                                prior(normal(0, 0.1), class = "b"),
                                prior(exponential(7), class = "sd")),
                      iter = 10, 
                      chains = 1)                                  
                                  
