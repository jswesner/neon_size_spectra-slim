library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
theme_set(brms::theme_default())

pparetocounts = function(x, xmin, xmax, lambda){(1 - (x^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1))))}

# plots will save here
directory = "plots/"

#1) Load models and data
# load data and models -----------------------------------------
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_updated12192024.rds")
dat = fit_temp_om_gpp$data

#2) Resample data
n_samples = 5000

dat_resampled_rank = dat %>% 
  group_by(sample_id) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number()),
         prob_yx = n_yx/max(row_number())) %>% 
  ungroup()

dat_resampled_rank %>%
  filter(site_id == "ARIK") %>% 
  ggplot(aes(x = dw, y = n_yx)) + 
  geom_point() +
  scale_y_log10() +
  scale_x_log10()

#3) ISD by sample ----------------------------------

epred_draws = dat_resampled_rank %>% 
  distinct(xmin, xmax, log_om_s, mat_s, log_gpp_s, sample_id, year, site_id) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)

epred_draws_summary = epred_draws %>% 
  group_by(sample_id, xmin, xmax) %>% 
  median_qi(.epred) %>% 
  pivot_longer(cols = c(.epred, .lower, .upper),
               names_to = "quantile", 
               values_to = "lambda")

epred_draws_list = epred_draws_summary %>% 
  group_by(sample_id) %>% 
  group_split()


isd_lines = epred_draws_summary %>%
  expand_grid(interval = 1:50) %>% 
  group_by(quantile, sample_id) %>% 
  mutate(x = seq(min(xmin), max(xmax), length.out = max(interval))) %>% 
  mutate(prob_yx = pparetocounts(x = x, xmin = xmin, xmax = xmax, lambda = lambda)) %>% 
  select(-.width, -.point, -.interval, -lambda) %>% 
  pivot_wider(names_from = quantile, values_from = prob_yx)


isd_lines %>% 
  filter(sample_id <= 20) %>%
  ggplot(aes(x = x, y = .epred)) + 
  geom_line(aes(group = sample_id)) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.001, 1)) +
  geom_point(data = dat_resampled_rank %>% 
               filter(sample_id <= 20), aes(y = prob_yx, x = dw)) +
  theme_void() +
  facet_wrap(~sample_id)






  
  

