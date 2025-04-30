library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
theme_set(brms::theme_default())

#1) Load models and data
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")
dat = as_tibble(fit_temp_om_gpp$data)

posts_sample_lambdas = fit_temp_om_gpp$data %>% 
  distinct(sample_id, mat_s, log_gpp_s, log_om_s, year, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp)

#2) Resample data
n_samples = 5000

dat_resampled_rank = dat %>% 
  group_by(sample_id, log_om_s, log_gpp_s, mat_s, site_id, year) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  select(site_id, year, sample_id, no_m2, dw) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number())) %>% 
  ungroup()

#3) ISD by sample ----------------------------------

epred_draws = dat_resampled_rank %>% 
  distinct(sample_id, log_om_s, log_gpp_s, mat_s, site_id, year, xmin, xmax) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)

epred_draws_summary = epred_draws %>% 
  group_by(sample_id, log_om_s, log_gpp_s, mat_s, site_id, year, xmin, xmax, no_m2) %>% 
  median_qi(.epred) %>% 
  pivot_longer(cols = c(.epred, .lower, .upper),
               names_to = "quantile", 
               values_to = "lambda")

epred_draws_list = epred_draws_summary %>% 
  group_by(sample_id, quantile) %>% 
  group_split()

temp = lapply(epred_draws_list, function(df) {
  df %>% expand_grid(x = 10^seq(log10(xmin), log10(xmax), length.out = 50)) %>% 
    mutate(prob_yx = (1 - (x^(lambda + 1) - (xmin^(lambda + 1))) / ((xmax)^(lambda + 1) - (xmin^(lambda + 1)))),
           n_yx = prob_yx * n_samples) 
})


site_ordered = dat_resampled_rank %>% 
  ungroup %>% 
  distinct(site_id, sample_id, year) %>%
  arrange(site_id, sample_id) %>% 
  mutate(site_order = 1:nrow(.),
         site_order_site = paste0(site_order, ") ", site_id))

isd_lines = bind_rows(temp) %>% 
  left_join(site_ordered)

sample_max = 133

dat_sliced = dat_resampled_rank %>% 
  filter(sample_id <= sample_max) %>%
  group_by(sample_id) %>% 
  # slice_sample(n = 1000) %>% 
  left_join(site_ordered)

plot_isd_sample = isd_lines %>% 
  filter(sample_id <= sample_max) %>%
  select(-prob_yx, -lambda) %>% 
  pivot_wider(names_from = quantile, values_from = n_yx) %>% 
  mutate(dw = x) %>%
  ggplot(aes(x = dw, y = .epred)) +
  geom_point(data = dat_sliced ,
             aes(y = n_yx,
                 color = site_id),
                 size = 0.05,
             shape = 1)  +
  geom_line() +
  geom_text(data = site_ordered %>% 
              filter(sample_id <= sample_max), aes(label = site_order_site),
            x = -1.1,
            y = 0.5,
            size = 1.2) +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.4) +
  scale_x_log10() +
  scale_y_log10() +
  coord_cartesian(ylim = c(0.5, NA)) +  
  facet_wrap(~site_order) +
  guides(size = "none",
         color = "none") +
  theme(axis.text = element_blank(),
        strip.text = element_blank(),
        panel.spacing = unit(0.1,'lines')) +
  labs(y = "\u03bb",
       x = "mgDM") +
  NULL

ggsave(plot_isd_sample, width = 6.5, height = 7,
       file = "plots/plot_isd_sample.jpg", dpi = 300)
saveRDS(plot_isd_sample, file = "plots/plot_isd_sample.rds")


#4) ISD by site ----------------------------------

n_site_samples = 5000

dat_resampled_rank_site = dat %>% 
  group_by(log_om_s, log_gpp_s, mat_s, site_id, year) %>% 
  sample_n(n_site_samples, weight = no_m2, replace = T) %>% 
  select(site_id, year, no_m2, dw) %>% 
  group_by(site_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  group_by(site_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number())) %>% 
  ungroup()
  
epred_draws_site = dat_resampled_rank_site %>% 
  distinct(log_om_s, log_gpp_s, mat_s, site_id, xmin, xmax) %>% 
  mutate(no_m2 = 1,
         site = "new",
         year = "new") %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, allow_new_levels = T)

epred_draws_summary_site = epred_draws_site %>% 
  group_by(log_om_s, log_gpp_s, mat_s, site_id, xmin, xmax, no_m2) %>% 
  median_qi(.epred) %>% 
  pivot_longer(cols = c(.epred, .lower, .upper),
               names_to = "quantile", 
               values_to = "lambda")

epred_draws_list_site = epred_draws_summary_site %>% 
  group_by(site_id, quantile) %>% 
  group_split()

temp_site = lapply(epred_draws_list_site, function(df) {
  df %>% expand_grid(x = 10^seq(log10(xmin), log10(xmax), length.out = 50)) %>% 
    mutate(prob_yx = (1 - (x^(lambda + 1) - (xmin^(lambda + 1))) / ((xmax)^(lambda + 1) - (xmin^(lambda + 1)))),
           n_yx = prob_yx * n_site_samples) 
})

isd_lines_site = bind_rows(temp_site)

dat_sliced_site = dat_resampled_rank_site %>% 
  group_by(site_id) %>% 
  # slice_sample(n = 5000) %>% 
  group_by(site_id) %>% 
  mutate(prob_yx = n_yx/max(n_yx))

plot_isd_site = isd_lines_site %>% 
  select(-n_yx, -lambda) %>% 
  pivot_wider(names_from = quantile, values_from = prob_yx) %>% 
  mutate(dw = x) %>%
  ggplot(aes(x = dw, y = .epred)) +
  geom_line() +
  geom_ribbon(aes(ymin = .lower, ymax = .upper),
              alpha = 0.4) +
  geom_point(data = dat_sliced_site,
             aes(y = prob_yx, size = dw,
                 color = dw),
             shape = 16) +
  scale_x_log10() +
  scale_y_log10() +
  # coord_cartesian(ylim = c(0.8, NA)) +
  facet_wrap(~site_id) +
  scale_color_viridis(trans = "log",
                      breaks = c(0.01,  1, 100, 10000),
                      labels = c("0.01", "1", "100","10,000"),
                      option = "D",
                      begin = 0,
                      end = 0.9,
                      alpha = 0.7) +
  guides(size = "none",
         color = "none") +
  NULL

ggsave(plot_isd_site, width = 6.5, height = 7,
       file = "plots/plot_isd.jpg", dpi = 300)
saveRDS(plot_isd_site, file = "plots/plot_isd.rds")

