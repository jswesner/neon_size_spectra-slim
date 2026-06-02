library(tidyverse)
library(brms)
library(tidybayes)
library(isdbayes)

# Literature Figure Comparisons -------------------------------------------
theme_set(theme_default())

predictors = readRDS("data/predictors.rds")

predictors_scaled

predictors %>% 
  distinct(site_id,
           mean_temp,
           sd_temp,
           mean_gpp,
           sd_gpp,
           mean_om,
           sd_om
           ) %>% 
  mutate(across(where(is.numeric), round, 0))%>% 
  write_csv(., file = "tables/environmental_summary.csv")


