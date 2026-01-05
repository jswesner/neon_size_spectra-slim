library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
theme_set(brms::theme_default())

# pparetocounts = function(x, xmin, xmax, lambda){(1 - (x^(lambda + 1) - (xmin^(lambda+1)))/(xmax^(lambda + 1) - (xmin^(lambda+1))))}

# plots will save here
directory = "plots/"

#1) Load models and data
# load data and models -----------------------------------------
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")
dat = fit_temp_om_gpp$data

#2) Resample data
n_samples = 5000


dat_resampled_rank = dat %>% 
  group_by(sample_id) %>% 
  slice_sample(
    prop = 1,              # sample 100% of rows in each group
    weight_by = no_m2,     # weighted by no_m2
    replace = TRUE) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number()),
         prob_yx = n_yx/max(row_number())) %>% 
  ungroup()


#3) ISD by sample ----------------------------------

epred_draws = dat_resampled_rank %>% 
  distinct(xmin, xmax, log_om_s, mat_s, log_gpp_s, sample_id, year, site_id) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)

epred_draws_summary = epred_draws %>% 
  group_by(sample_id, xmin, xmax, site_id) %>% 
  median_qi(.epred) %>% 
  pivot_longer(cols = c(.epred, .lower, .upper),
               names_to = "quantile", 
               values_to = "lambda")

epred_draws_list = epred_draws_summary %>% 
  group_by(sample_id, site_id) %>% 
  group_split()

isd_lines = epred_draws_summary %>%
  expand_grid(interval = 1:50) %>% 
  group_by(quantile, sample_id) %>% 
  mutate(x = exp(seq(min(log(xmin)), max(log(xmax)), length.out = max(interval)))) %>%
  # mutate(x = seq(min(xmin), max(xmax), length.out = max(interval))) %>% 
  mutate(prob_yx = pparetocounts(x = x, 
                                 xmin = xmin - 0.001, # this subtraction is a workaround. for some reason
                                 # the exp(log(...)) procedure produces prob_yx = 0 when x = xmin, but it doesn't do that using the linear procedure.
                                 # this fix ensures that xmin is smaller than x
                                 xmax = xmax, lambda = lambda)) %>% 
  select(-.width, -.point, -.interval, -lambda) %>% 
  pivot_wider(names_from = quantile, values_from = prob_yx)


set.seed(204)
id_pull = isd_lines %>% ungroup %>% distinct(sample_id) %>% 
  sample_n(21) %>% pull(sample_id)

isd_fits = isd_lines %>% 
  filter(sample_id %in% id_pull) %>%
  ggplot(aes(x = x, y = .epred, group = interaction(site_id, sample_id))) +
  geom_point(data = dat_resampled_rank %>%
               filter(sample_id %in% id_pull), aes(y = prob_yx, x = dw),
             size = 0.2) + 
  geom_line(color = "red") +
  geom_ribbon(aes(ymin = .lower, ymax = .upper), 
              fill = "red") +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.001, 1)) +
  facet_wrap(~interaction(site_id,sample_id), ncol = 3) +
  labs(y = "Pr(X \u2265 x)",
       x = "Individual Dry Mass (mg)") +
  theme(strip.text = element_text(size = 9))

ggsave(isd_fits, file = "plots/isd_fits.jpg",
       width = 6, height = 8, dpi = 400)






  
  

