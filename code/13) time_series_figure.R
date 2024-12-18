library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggh4x)

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")

fit_temp_om_gpp$preds = "temp*om*gpp"

# load data
dat_all = readRDS("data/derived_data/dat_all.rds") %>% mutate(temp_mean = mean)
climate_zones = read_csv("data/classified_coordinates.csv") %>%
  clean_names() %>% 
  distinct(site_id, climate_zone) %>% 
  mutate(zone_order = as.integer(as.factor(climate_zone)),
         zone_letter = letters[zone_order],
         zone_title = paste0(zone_letter, ") Zone ", climate_zone))

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
  left_join(climate_zones) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)

sample_posts_summary = sample_posts %>% 
  group_by(sample_id, site_id, date, mat_s, climate_zone, zone_title) %>% 
  median_qi(.epred)

one_site_per_zone = climate_zones %>% 
  group_by(climate_zone, zone_title) %>% 
  slice(1)

global_mean_lambda = as_draws_df(fit_temp_om_gpp) %>% select(b_Intercept, starts_with("sd_")) %>% 
  mutate(.epred = b_Intercept +
           rnorm(nrow(.), 0, sd_sample_id__Intercept) +
           rnorm(nrow(.), 0, sd_site_id__Intercept)+
           rnorm(nrow(.), 0, sd_year__Intercept)) %>% 
  expand_grid(site_id = one_site_per_zone %>% distinct(site_id) %>% pull()) %>%
  left_join(climate_zones) %>%
  expand_grid(date = seq(min(dat_all$date),
                         max(dat_all$date),
                         length.out = 2)) %>% 
  group_by(site_id, date, climate_zone, zone_title) %>% 
  reframe(low50 = quantile(.epred, probs = 0.25),
          high50 = quantile(.epred, probs = 0.75),
          low95 = quantile(.epred, probs = 0.025),
          high95 = quantile(.epred, probs = 0.975),
          .epred = median(.epred))


labels = dat_all %>% group_by(site_id) %>% filter(date == max(date))

library(viridis)
time_series_fig = sample_posts_summary %>% 
  ggplot(aes(x = date, y = .epred)) +
  # stat_lineribbon(data = global_mean_lambda, alpha = 0.1, .width = c(0.95), fill = "black",
                  # linewidth = 0.2) +
  # stat_lineribbon(data = global_mean_lambda, alpha = 0.1, .width = c(0.5), fill = "black",
                  # linewidth = 0) +
  geom_ribbon(data = global_mean_lambda, aes(ymin = low50, 
                                             ymax = high50,
                                             fill = climate_zone),
                  alpha = 0.25) +
  geom_ribbon(data = global_mean_lambda, aes(ymin = low95,
                                             ymax = high95, 
                                             fill = climate_zone),
              alpha = 0.25) +
  geom_line(data = global_mean_lambda, linewidth = 0.2, aes(color = climate_zone)) +
  geom_point(aes(group = site_id), 
                  color = "black",
                  size = 0.5,
                  linewidth = 0.2) +
  geom_line(aes(group = site_id), color = "black", linewidth = 0.1) +
  facet_wrap2(zone_title~.) +
  theme_default() +
  scale_x_date(date_labels = "%Y") +
  ylim(-4, -1) +
  guides(color = "none", 
         fill = "none") +
  labs(y = "\u03bb",
       x = "Collection Date") +
  scale_color_viridis() +
  scale_fill_viridis()

ggsave(time_series_fig, file = "plots/time_series_fig.jpg", width = 6, height = 7, dpi = 400)
