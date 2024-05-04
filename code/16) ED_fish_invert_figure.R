library(rstan)
library(tidybayes)
library(brms)
library(tidyverse)


fish_dw <- readRDS("data/derived_data/fish_dw-wrangled.rds")
fit_fishonly_temp_om_gpp = readRDS(file = "models/fit_fishonly_temp_om_gpp.rds")



fish_posts = as_draws_df(fit_fishonly_temp_om_gpp) %>% select(.draw, a, starts_with("beta"))
fish_params = names(fish_posts)[1:8]


fish_lines = tibble(mat_s = seq(min(fish_dw$mat_s), max(fish_dw$mat_s), length.out = 20)) %>% 
  mutate(log_gpp_s = quantile(fish_dw$log_gpp_s, probs = 0.25),
         log_om_s = quantile(fish_dw$log_om_s, na.rm = T, probs = 0.75)) %>% 
  merge(fish_posts) %>% 
  as_tibble() %>% 
  mutate(.epred = a + mat_s*beta_mat + mat_s*log_gpp_s*beta_gpp_mat + mat_s*log_gpp_s*log_om_s*beta_om_mat_gpp)


fish_lines %>% 
  ggplot(aes(x = mat_s, y = .epred)) + 
  stat_lineribbon(.width = 0.95)

