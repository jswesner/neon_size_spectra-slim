library(brms)
library(isdbayes)
library(tidyverse)
library(janitor)

# Instead of re-fitting, you can load the fitted models below

# load models (these are already fit)
# fit_temp = readRDS("models/fit_temp.rds")
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")

# Re-fit models -----------------------------------------------------------
# !!!! Each of the models below takes up to 24 hours to run. They were run on a cluster. 

#1) load data (NEON body size data)
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")

#2) fit main model
# fit_temp_om_gpp_year = brm(dw | vreal(no_m2, xmin, xmax) ~ log_om_s*mat_s*log_gpp_s + (1 |site_id:sample_id) + (1 + log_om_s*mat_s*log_gpp_s|year),
#                            data = dat_2022_clauset,
#                            stanvars = stanvars,    # required for truncated Pareto via isdbayes package
#                            family = paretocounts(),# required for truncated Pareto via isdbayes package
#                            prior = c(prior(normal(-2, 0.5), class = "Intercept"),
#                                      prior(normal(0, 0.2), class = "b"),
#                                      prior(exponential(7), class = "sd")),
#                            iter = 2000,
#                            cores = 4,
#                            threads = 12,
#                            chains = 4)
# 
# saveRDS(fit_temp_om_gpp_year, file = "models/fit_temp_om_gpp_year.rds")

#2) fit submodel using the update function

# fit_temp_year = update(fit_temp_om_gpp, formula = . ~ mat_s + (1 | site_id:sample_id) + (1 + mat_s | year),
#                        cores = 4, threads = 12, chains = 4, iter = 2000)
# saveRDS(fit_temp_year, file = "models/fit_temp_year.rds")



# check seasonal effect ---------------------------------------------------

# first samples of the year only
dat_first = dat_2022_clauset %>%
  group_by(site_id, year) %>%
  mutate(n = length(unique(collect_date))) %>%
  filter(n > 1) %>% # only sites/years with more than 1 collection
  filter(collect_date == min(collect_date))

# last samples of the year only
dat_last = dat_2022_clauset %>%
  group_by(site_id, year) %>%
  mutate(n = length(unique(collect_date))) %>%
  filter(n > 1) %>%
  filter(collect_date == max(collect_date))

# fit_first = update(fit_temp_om_gpp, newdata = dat_first, chains = 1, iter = 1000, data2 = list(model = "first"))
# fit_last = update(fit_temp_om_gpp, newdata = dat_last, chains = 1, iter = 1000, data2 = list(model = "last"))

# saveRDS(fit_first, file = "models/fit_first.rds")
# saveRDS(fit_last, file = "models/fit_last.rds")


# fit with just macros --------------------------------------------------
# not that the macros only data set has more observations than the macro/fish dataset
# because the culling procedure keeps some smaller individuals with macros only, which leads to 
# higher n.
# dat_macros = readRDS("data/dat_2022_macros_clauset.rds")
# 
# fit_macros = update(fit_temp_om_gpp, newdata = dat_macros, chains = 1, iter = 1000,
#                     data2 = list(model = "macros"))
# 
# saveRDS(fit_macros, file = "models/fit_macros.rds")



# fit_with just fish ------------------------------------------------------

dat_fish = readRDS(file = "data/dat_2022_fish_clauset.rds")
fit_fish = update(fit_temp_om_gpp, newdata = dat_fish, 
                  chains = 1, iter = 1000, data2 = list(model = "last"))

saveRDS(fit_fish, file = "models/fit_fish.rds")


# add Gaussian process for spatial autocorrelation ---------------------------------------------

fit_temp_om_gpp_year = readRDS(file = "models/fit_temp_om_gpp_year.rds")

# load data
neon_latlong <- read_csv(file = "data/raw_data/site_lat_longs.csv") %>% distinct(siteID, lat, long) %>% 
  clean_names() 

dat_lat_long = left_join(fit_temp_om_gpp_year$data, neon_latlong)

fit_latlong = update(fit_temp_om_gpp_year, 
                     formula = . ~ log_om_s * mat_s * log_gpp_s + (1 | site_id:sample_id) + (1 + log_om_s * mat_s * log_gpp_s | year) +
                       gp(lat, long), newdata = dat_lat_long, 
                     iter = 2000, chains = 4)

saveRDS(fit_latlong, file = "models/fit_latlong.rds")

