library(tidybayes)
library(tidyverse)
library(isdbayes)

model_list = readRDS(file = "models/model_list.rds")

# mod = model_list$`models/fit_temp_newxmin_sumnorm.rds`
mod = model_list$`models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds`
dat = mod$data

post_dots = dat %>% 
  select(-dw, -no_m2) %>% distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(mod, re_formula = NULL)

post_lines = tibble(mat_s = seq(min(dat$mat_s),
                        max(dat$mat_s),
                        length.out = 30)) %>% 
  mutate(log_gpp_s = 0,
         log_om_s = 0,
         no_m2 = 1,
         xmin = min(dat$dw),
         xmax = max(dat$dw)) %>% 
  add_epred_draws(mod, re_formula = NA)

post_dots %>% 
  # mutate(cari = case_when(site_id == "CARI" ~ "CARI",
                          # TRUE ~ "Other")) %>% 
  ggplot(aes(x = mat_s, y = .epred)) +
  stat_pointinterval(aes(group = sample_id), size = 1) +
  stat_lineribbon(data = post_lines, .width = 0.95, alpha = 0.4)

# factorial plot ----------------------------------------------------------
# make grid
mod = model_list$`models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds`
dat = mod$data

preds_raw = mod$data %>% select(ends_with("_s"), starts_with("log"))
quantiles = c(0.25, 0.5, 0.75)

grid_a = expand_grid(mat_s = seq(min(preds_raw$mat_s), max(preds_raw$mat_s), length.out = 30),
                     log_gpp_s = quantile(preds_raw$log_gpp_s, probs = quantiles),
                     log_om_s = quantile(preds_raw$log_om_s, probs = quantiles)) %>% 
                       mutate(predictor = "mat_s",
                              gpp_group = case_when(log_gpp_s == min(log_gpp_s) ~ "gpp25",
                                                    log_gpp_s == max(log_gpp_s) ~ "gpp75",
                                                    T ~ "gpp50"),
                              om_group = case_when(log_om_s == min(log_om_s) ~ "om25",
                                                    log_om_s == max(log_om_s) ~ "om75",
                                                    T ~ "om50")) %>% 
  mutate(xmin = min(dat$xmin),
         xmax = max(dat$xmax),
         no_m2 = 1)

grid_b = expand_grid(log_gpp_s = seq(min(preds_raw$log_gpp_s), max(preds_raw$log_gpp_s), length.out = 30),
                     mat_s = quantile(preds_raw$mat_s, probs = quantiles),
                     log_om_s = quantile(preds_raw$log_om_s, probs = quantiles)) %>% 
  mutate(predictor = "log_gpp_s",
         mat_group = case_when(mat_s == min(mat_s) ~ "mat25",
                               mat_s == max(mat_s) ~ "mat75",
                               T ~ "mat50"),
         om_group = case_when(log_om_s == min(log_om_s) ~ "om25",
                              log_om_s == max(log_om_s) ~ "om75",
                              T ~ "om50")) %>% 
  mutate(xmin = min(dat$xmin),
         xmax = max(dat$xmax),
         no_m2 = 1)

grid_c = expand_grid(log_om_s = seq(min(preds_raw$log_om_s), max(preds_raw$log_om_s), length.out = 30),
                     mat_s = quantile(preds_raw$mat_s, probs = quantiles),
                     log_gpp_s = quantile(preds_raw$log_gpp_s, probs = quantiles)) %>% 
  mutate(predictor = "log_om_s",
         mat_group = case_when(mat_s == min(mat_s) ~ "mat25",
                               mat_s == max(mat_s) ~ "mat75",
                               T ~ "mat50"),
         gpp_group = case_when(log_gpp_s == min(log_gpp_s) ~ "gpp25",
                              log_gpp_s == max(log_gpp_s) ~ "gpp75",
                              T ~ "gpp50"))%>% 
  mutate(xmin = min(dat$xmin),
         xmax = max(dat$xmax),
         no_m2 = 1)

post_a = grid_a %>% mutate(group = paste(parse_number(gpp_group), 
                                         parse_number(om_group), sep = "_")) %>% add_epred_draws(mod, re_formula = NA) %>% mutate(x = mat_s)
post_b = grid_b %>% mutate(group = paste(parse_number(mat_group),
                                         parse_number(om_group), sep = "_")) %>% add_epred_draws(mod, re_formula = NA) %>% mutate(x = log_gpp_s)
post_c = grid_c %>% mutate(group = paste(parse_number(mat_group), 
                                         parse_number(gpp_group), sep = "_")) %>% add_epred_draws(mod, re_formula = NA) %>% mutate(x = log_om_s)

bind_rows(post_a, post_b, post_c) %>% 
  filter(.draw <= 100) %>% 
  ggplot(aes(x = x, y = .epred, fill = group)) + 
  stat_lineribbon(.width = 0.95, alpha = 0.4) +
  ggh4x::facet_grid2(group~predictor) +
  guides(fill = "none")
    
post_a %>% 
  # filter(.draw <= 100) %>% 
  ggplot(aes(x = mat_s, y = .epred, fill = group)) + 
  stat_lineribbon(.width = 0.95, alpha = 0.4) +
  ggh4x::facet_grid2(gpp_group ~ om_group) +
  guides(fill = "none")

post_b %>% 
  # filter(.draw <= 100) %>% 
  ggplot(aes(x = log_gpp_s, y = .epred, fill = group)) + 
  stat_lineribbon(.width = 0.95, alpha = 0.4) +
  ggh4x::facet_grid2(mat_group ~ om_group) +
  guides(fill = "none")

post_c %>% 
  # filter(.draw <= 100) %>% 
  ggplot(aes(x = log_om_s, y = .epred, fill = group)) + 
  stat_lineribbon(.width = 0.95, alpha = 0.4) +
  ggh4x::facet_grid2(gpp_group ~ mat_group) +
  guides(fill = "none")

slopes = bind_rows(post_a, post_b, post_c) %>% 
  group_by(predictor, group) %>% 
  filter(x == min(x) | x == max(x)) %>% 
  select(group, predictor, .draw, x, .epred) %>% 
  mutate(x_range = max(x) - min(x),
         x_minmax = case_when(x == min(x) ~ "min",
                              TRUE ~ "max")) %>% 
  select(-x) %>% 
  pivot_wider(names_from = x_minmax, values_from = .epred) %>% 
  mutate(slope = (max-min)/x_range)

slopes %>% 
  group_by(predictor, .draw) %>% 
  reframe(slope = mean(slope)) %>% 
  ggplot(aes(x = slope, y = predictor)) +
  stat_slab(alpha = 0.4) +
  geom_vline(aes(xintercept = 0))

slopes %>% 
  group_by(predictor, .draw) %>% 
  reframe(slope = mean(slope)) %>% 
  group_by(predictor) %>% 
  reframe(prob_pos = sum(slope>0)/max(.draw))
