library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")

fit_temp_om_gpp$preds = "temp*om*gpp"

# load data
dat_all = readRDS("data/derived_data/dat_all.rds") %>% mutate(temp_mean = mean)

mean_temp = mean(unique(dat_all$temp_mean))
sd_temp = sd(unique(dat_all$temp_mean))
mean_om = mean(unique(dat_all$log_om))
sd_om = sd(unique(dat_all$log_om))
mean_gpp = mean(unique(dat_all$log_gpp))
sd_gpp = sd(unique(dat_all$log_gpp))

sample_posts = fit_temp_om_gpp$data %>% glimpse %>% 
  select(-dw, -no_m2) %>% distinct() %>% 
  mutate(no_m2 = 1) %>% 
  left_join(dat_all %>% ungroup %>% distinct(date, sample_id)) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)

sample_posts_summary = sample_posts %>% group_by(sample_id, site_id, date, mat_s) 

global_mean_lambda = as_draws_df(fit_temp_om_gpp) %>% select(b_Intercept, starts_with("sd_")) %>% 
  mutate(.epred = b_Intercept +
           rnorm(nrow(.), 0, sd_sample_id__Intercept) +
           rnorm(nrow(.), 0, sd_site_id__Intercept)+
           rnorm(nrow(.), 0, sd_year__Intercept)) %>% 
  expand_grid(site_id = dat_all %>% distinct(site_id) %>% pull()) %>% 
  expand_grid(date = seq(min(dat_all$date),
                         max(dat_all$date),
                         length.out = 2))

time_series_fig = sample_posts_summary %>% 
  ggplot(aes(x = date, y = .epred, group = site_id)) +
  stat_lineribbon(data = global_mean_lambda, alpha = 0.1, .width = c(0.95), fill = "black",
                  linewidth = 0.2) +
  stat_lineribbon(data = global_mean_lambda, alpha = 0.1, .width = c(0.5), fill = "black",
                  linewidth = 0) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper),
                  size = 0.2) +
  geom_line() +
  facet_wrap(~site_id) +
  theme_default() +
  scale_x_date(date_labels = "%Y") +
  ylim(-4, 0) +
  guides(color = "none") +
  labs(y = "\u03bb",
       x = "Collection Date")

ggsave(time_series_fig, file = "plots/time_series_fig.jpg", width = 6, height = 9, dpi = 400)
