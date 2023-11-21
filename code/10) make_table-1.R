library(tidyverse)
library(brms)
library(tidybayes)
library(isdbayes)

# Literature Figure Comparisons -------------------------------------------
theme_set(theme_default())

dat_all = readRDS("data/derived_data/dat_all.rds")

dat_all %>% 
  distinct(temp_mean, temp_sd,
           gpp, gpp_sd,
           mean_om, sd_om) %>% 
  mutate(across(where(is.numeric), round, 0))%>% 
  write_csv(., file = "tables/environmental_summary.csv")
