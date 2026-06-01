# tests how lambda varies with sample size
# n = unique body sizes in a sample from sample_size.rds
library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
theme_set(brms::theme_default())

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")

# load data
sample_size = readRDS("data/sample_size.rds")
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
dat_all = readRDS(file = "data/dat_2022_clauset.rds")

post_medians = fit_temp_om_gpp$data %>% select(-dw, -no_m2) %>% 
  distinct() %>%
  mutate(no_m2 = 1) %>% 
  left_join(sample_size) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~ (1|sample_id)) %>% 
  group_by(sample_id, n) %>% 
  median_qi(.epred)


post_medians %>% 
  arrange(.epred) %>% 
  ggplot(aes(x = n, y = .epred)) +
  geom_point()




