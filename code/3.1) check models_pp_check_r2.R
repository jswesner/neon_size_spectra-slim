library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(patchwork)
theme_set(brms::theme_default())

# 1) load models
mod = readRDS("models/fit_temp_om_gpp.rds")

# 2) get data
dat = as_tibble(mod$data)

# resampling before running pp_check to remove influence of no_m2 
dat_resampled = dat %>% ungroup %>% sample_n(10000, weight = no_m2, replace = T)
dat_resampled_site = dat %>% group_by(site_id) %>% sample_n(5, weight = no_m2, replace = T)

pp_check_global = pp_check(mod, newdata = dat_resampled) + scale_x_log10()
pp_check_sites = pp_check(mod, newdata = dat_resampled_site, group = "site_id", type = "dens_overlay_grouped") + scale_x_log10()

pp_check_plot = (pp_check_global + theme(legend.position = c(0.8, 0.9)) +
    labs(subtitle = "a)"))/
  (pp_check_sites + guides(color = "none") + labs(subtitle = "b)") +
     theme(axis.text.x = element_text(angle =90, vjust = 0.5, hjust = 0))) 

ggsave(pp_check_plot, file = "plots/fig_s6_pp_check_plot.jpg", width = 6.5, height = 9, dpi = 400)

# 
# # 3) extract posteriors 
# posts_sample_lambdas = dat %>% 
#   distinct(sample_id, .keep_all = T) %>% 
#   select(-dw) %>% 
#   mutate(no_m2 = 1) %>% 
#   add_epred_draws(mod, re_formula = NULL) %>% 
#   rename(lambda = .epred) %>% 
#   ungroup()
# 
# saveRDS(posts_sample_lambdas, file = "posteriors/posts_sample_lambdas.rds")
# 
# # 4) merge posts and raw data
# n_samples = 5000
# 
# dat_resampled = dat %>% 
#   group_by(sample_id, log_om_s, log_gpp_s, mat_s, site_id, year) %>% 
#   sample_n(n_samples, weight = no_m2, replace = T) %>% 
#   select(site_id, year, sample_id, no_m2, dw) %>% 
#   group_by(sample_id) %>% 
#   mutate(xmin = min(dw),
#          xmax = max(dw),
#          data = "y_raw") %>% 
#   ungroup()
# 
# posts_raw = posts_sample_lambdas %>% 
#   filter(.draw <= 10) %>% 
#   select(lambda, sample_id, .draw) %>% 
#   right_join(dat_resampled %>% 
#                group_by(sample_id) %>% 
#                slice_sample(n = 500) %>% 
#                select(sample_id, dw, no_m2, xmin, xmax, site_id), multiple = "all",
#              relationship = "many-to-many") 
# 
# # 5) sample posterior preds
# sim_posts = posts_raw %>% 
#   mutate(x = rparetocounts(lambda = lambda, xmin = xmin, xmax = xmax, n = nrow(.))) %>% 
#   mutate(data = "y_rep") %>% 
#   # mutate(x = round(x,2)) %>% 
#   bind_rows(dat_resampled %>% 
#               mutate(data = "y") %>% 
#               mutate(.draw = 0,
#                      x = dw)) %>% 
#   rename(sim = x)
# 
# # 6) pick a sample to plot
# id = as.integer(runif(15, 1, length(unique(dat$sample_id))))
# 
# # 7) Make Plots
# # violin
# post_pred_si = sim_posts %>%
#   # filter(sample_id %in% id) %>%
#   ggplot(aes(x = as.integer(.draw), y = sim, color = data, group = .draw)) + 
#   geom_violin() +
#   scale_y_log10() +
#   scale_x_continuous(breaks = c(0,1, 2,3, 4,5, 6,7, 8,9, 10),
#                      labels = c(expression(italic("y")), 
#                                 "1", "2","3", "4","5", "6","7", "8","9", "10")) +
#   facet_wrap(~ site_id, ncol = 4) +
#   scale_color_colorblind(labels = expression(italic("y"[rep]) ~ "rep")) +
#   labs(y = "Individual Body Size (mgDM)",
#        x = "Draw",
#        color = "") +
#   guides(color = "none") +
#   theme(legend.position = c(0.8, 0.08)) +
#   NULL
# 
# # ggview::ggview(post_pred_si, width = 6.5, height = 9)
# ggsave(post_pred_si, file = "plots/post_pred_i.jpg", width = 6.5, height = 9)
# 
# sim_posts %>%
#   filter(site_id == "ARIK") %>% 
#   # filter(sample_id %in% id) %>%
#   ggplot(aes(x = .draw, y = sim, color = data, group = .draw)) + 
#   geom_violin() +
#   scale_y_log10() +
#   facet_wrap(~ sample_id) +
#   NULL
# 
# # density
# sim_posts %>%
#   ggplot(aes(x = sim, color = data, group = .draw)) + 
#   geom_density() +
#   scale_x_log10() +
#   facet_wrap(~site_id, scales = "free_y") +
#   scale_color_colorblind() +
#   NULL

# r2 ----------------------------------------------------------------------
# from: Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for Bayesian regression models. The American Statistician.
n_samples = 5000
dat_resampled = dat %>% 
  group_by(sample_id, log_gpp_s, log_om_s, mat_s, site_id, year) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  select(site_id, year, sample_id, no_m2, dw) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  ungroup()

set.seed = 23234

preds = dat_resampled %>% 
  slice_sample(n = 5000) %>% 
  mutate(no_m2 = 1) %>% 
  add_predicted_draws(mod, re_formula = NA, ndraws = 1000) 

preds %>% 
  group_by(.draw) %>% 
  reframe(v_fit = var(.prediction),
          v_pred = var(dw - .prediction)) %>% 
  mutate(r2 = v_fit/(v_fit + v_pred)) %>% 
  reframe(mean = mean(r2),
          sd = sd(r2))


# confirm r2 - the bayes_R2() function in brms does not work with isd models. The code below checks that our calculation of R2 above is correct by
# comparing the hand-calculated R2 from a basic stan model to the bayes_R2(). There is an exact match, suggesting that our calculation above is correct.

# mod_test = brm(mpg ~ hp, data = mtcars,
#                family = "Gaussian")

# as_tibble(t(fitted(mod_test, summary = F))) %>%
#   mutate(y = mod_test$data$mpg) %>%
#   pivot_longer(cols = -y) %>%
#   mutate(.draw = parse_number(name)) %>%
#   group_by(.draw) %>%
#   reframe(v_fit = var(value),
#           v_pred = var(y - value)) %>%
#   mutate(r2 = v_fit/(v_fit + v_pred)) %>%
#   mean_qi(r2)
# 
# bayes_R2(mod_test)

# compare to model with temperature only ----------------------------------

# check for temp model only

mod_temp = readRDS("models/fit_temp_nosyca_wideprior.rds")
dat_temp = mod_temp$data

posts_sample_lambdas_temp = dat_temp %>% 
  distinct(sample_id, .keep_all = T) %>% 
  select(-dw) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(mod_temp, re_formula = NULL) %>% 
  rename(lambda = .epred) %>% 
  ungroup()

# 4) merge posts and raw data
n_samples = 5000

dat_resampled_temp = dat_temp %>% 
  group_by(sample_id, mat_s, site_id, year) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  select(site_id, year, sample_id, no_m2, dw) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  ungroup()

posts_raw_temp = posts_sample_lambdas_temp %>% 
  filter(.draw <= 10) %>% 
  select(lambda, sample_id, .draw) %>% 
  right_join(dat_resampled_temp %>% 
               group_by(sample_id) %>% 
               slice_sample(n = 500) %>% 
               select(sample_id, dw, no_m2, xmin, xmax, site_id), multiple = "all",
             relationship = "many-to-many") 

# 5) sample posterior preds
sim_posts_temp = posts_raw_temp %>% 
  mutate(x = rparetocounts(lambda = lambda, xmin = xmin, xmax = xmax, n = nrow(.))) %>% 
  mutate(data = "y_rep") %>% 
  # mutate(x = round(x,2)) %>% 
  bind_rows(dat_resampled_temp %>% 
              mutate(data = "y") %>% 
              mutate(.draw = 0,
                     x = dw)) %>% 
  rename(sim = x)

# 6) pick a sample to plot
id = as.integer(runif(15, 1, length(unique(dat_temp$sample_id))))

# 7) Make Plots
# violin
post_pred_si_temp = sim_posts_temp %>%
  # filter(sample_id %in% id) %>%
  ggplot(aes(x = as.integer(.draw), y = sim, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  scale_x_continuous(breaks = c(0,1, 2,3, 4,5, 6,7, 8,9, 10),
                     labels = c(expression(italic("y")), 
                                "1", "2","3", "4","5", "6","7", "8","9", "10")) +
  facet_wrap(~ site_id, ncol = 4) +
  scale_color_colorblind(labels = expression(italic("y"[rep]) ~ "rep")) +
  labs(y = "Individual Body Size (mgDM)",
       x = "Draw",
       color = "") +
  guides(color = "none") +
  theme(legend.position = c(0.8, 0.08)) +
  NULL

# ggview::ggview(post_pred_si, width = 6.5, height = 9)
# ggsave(post_pred_si, file = "plots/post_pred_i.jpg", width = 6.5, height = 9)

sim_posts_temp %>%
  filter(site_id == "ARIK") %>% 
  # filter(sample_id %in% id) %>%
  ggplot(aes(x = .draw, y = sim, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~ sample_id) +
  NULL

# density
sim_posts_temp %>%
  ggplot(aes(x = sim, color = data, group = .draw)) + 
  geom_density() +
  scale_x_log10() +
  facet_wrap(~site_id, scales = "free_y") +
  scale_color_colorblind() +
  NULL


# r2
dat_resampled %>% 
  slice_sample(n = 5000) %>% 
  mutate(no_m2 = 1) %>% 
  add_predicted_draws(mod_temp, re_formula = NA, ndraws = 1000) %>% 
  group_by(.draw) %>% 
  reframe(v_fit = var(.prediction),
          v_pred = var(dw - .prediction)) %>% 
  mutate(r2 = v_fit/(v_fit + v_pred)) %>% 
  reframe(mean = mean(r2),
          sd = sd(r2))
