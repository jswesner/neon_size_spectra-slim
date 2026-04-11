# read in coordinates
library(tidyverse)
library(ape)
library(sf)


site_predictors <- readRDS("./data/predictors.rds") %>% 
  rename(siteID = 'site_id')

site_df = read.csv(file = "./data/classified_coordinates.csv", header = TRUE) %>%
  distinct %>% 
  left_join(site_predictors,.) %>% 
  st_as_sf(., coords = c("long", "lat"), crs = 4326)

site_dist = site_df %>% 
  st_distance() %>% 
  units::set_units(value = "km") %>% 
  as.data.frame

site_inv_dist = as.matrix(1/site_dist)
diag(site_inv_dist) <- 0
site_inv_dist_num = as.numeric(site_inv_dist)

#  run the Morans test on
# OM
ape::Moran.I(as.vector(site_df$log_om_s), site_inv_dist)

# GPP
ape::Moran.I(as.vector(site_df$log_gpp_s), site_inv_dist)

# temperature
ape::Moran.I(as.vector(site_df$mat_s), site_inv_dist)

# residuals
## get model 
fit_temp_year = readRDS(file = "models/fit_temp_year.rds")

resids_temp_year = dat_2022_clauset %>% tidybayes::add_residual_draws(fit_temp_year)
