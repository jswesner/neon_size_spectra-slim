library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
theme_set(brms::theme_default())

# For these data and models, pp_check() and other functions in brms will not work for model checking. This code does the same thing as pp_check, but it first 
# accounts for the no_m2 for each body size. This is necessary, because the raw size data "dw" will not match the posterior predictions without first correcting 
# for no_m2. Because the range of no_m2 is so large, we have to resample dw with replacement, and then plot that resampled data against model predictions.

# 1) load models
# load models 
fit_temp = readRDS("models/fit_temp.rds")
fit_om = readRDS("models/fit_om.rds")
fit_gpp = readRDS("models/fit_gpp.rds")
fit_temp_om = readRDS("models/fit_temp_om.rds")
fit_temp_gpp = readRDS("models/fit_temp_gpp.rds")
fit_om_gpp = readRDS("models/fit_om_gpp.rds")
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

# 2) get data
dat = as_tibble(fit_temp_om_gpp$data)

# 3) extract posteriors 
posts_sample_lambdas = dat %>% 
  distinct(sample_id, .keep_all = T) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  rename(lambda = .epred) %>% 
  ungroup()

# 4) merge posts and raw data
n_samples = 5000

dat_resampled = dat %>% 
  group_by(sample_id) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  select(site_id, year, sample_id, no_m2, dw) %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         data = "y_raw") %>% 
  ungroup()

posts_raw = posts_sample_lambdas %>% 
  filter(.draw <= 10) %>% 
  select(lambda, sample_id, .draw) %>% 
  right_join(dat_resampled %>% select(sample_id, dw, no_m2, xmin, xmax, site_id), multiple = "all",
             relationship = "many-to-many") 

# 5) sample posterior preds
sim_posts = posts_raw %>% 
  mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
  mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  mutate(data = "y_rep") %>% 
  # mutate(x = round(x,2)) %>% 
  bind_rows(dat_resampled %>% 
              mutate(data = "y_raw") %>% 
              mutate(.draw = 0,
                     x = dw)) %>% 
  rename(sim = x)

# 6) pick a sample to plot
id = as.integer(runif(15, 1, length(unique(dat$sample_id))))

#7) Make plots
# violin
sim_posts %>%
  # filter(sample_id %in% id) %>%
  ggplot(aes(x = .draw, y = sim, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~ site_id) +
  NULL

# density
sim_posts %>%
  ggplot(aes(x = sim, color = data, group = .draw)) + 
  geom_density() +
  scale_x_log10() +
  facet_wrap(~site_id, scales = "free_y") +
  scale_color_colorblind() +
  NULL

# bayes pvalue ------------------------------------------------------------
id = as.integer(runif(15, 1, length(unique(dat$sample_id))))

posts_raw_p = posts_sample_lambdas %>% 
  # filter(sample_id %in% id) %>% 
  filter(.draw <= 100) %>% 
  select(lambda, sample_id, .draw) %>% 
  right_join(dat_resampled %>% distinct(sample_id, dw, no_m2, xmin, xmax, site_id), 
             by = "sample_id", multiple = "all",
             relationship = "many-to-many") 

# sample posterior preds
n_samples = 1000
posts_raw_preds = posts_raw_p %>% 
  # filter(.draw <= 4) %>%
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
  mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  mutate(data = "y_rep") %>% 
  # mutate(x = round(x,2)) %>% 
  rename(sim = x) %>% 
  group_by(sample_id, .draw) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  ungroup

# by site
t_stat_post_site = posts_raw_preds %>% 
  # filter(site_id %in% id) %>% 
  group_by(.draw, site_id) %>% 
  reframe(gm = exp(mean(log(sim))),
          median = median(sim))

t_stat_raw_site = dat_resampled %>% 
  group_by(site_id) %>%
  reframe(gm_raw = exp(mean(log(dw))),
          median_raw = median(dw))


t_stat_post_site %>% 
  ggplot(aes(x = gm)) + 
  geom_histogram() + 
  facet_wrap(~site_id) + 
  geom_vline(data = t_stat_raw_site, aes(xintercept = gm_raw))


# by sample
t_stat_post_sample = posts_raw_preds %>% 
  # filter(sample_id %in% id) %>% 
  group_by(.draw, sample_id) %>% 
  reframe(gm = exp(mean(log(sim))),
          median = median(sim))

t_stat_raw_sample = dat_resampled %>% 
  group_by(sample_id) %>%
  reframe(gm_raw = exp(mean(log(dw))),
          median_raw = median(dw))

t_stat_summary = t_stat_post_sample %>% 
  # filter(sample_id == 2) %>% 
  # group_by(sample_id) %>% 
  # reframe(mean_gm = mean(gm),
  # sd_gm = sd(gm)) %>% 
  left_join(t_stat_raw_sample) %>% 
  left_join(dat %>% ungroup %>% distinct(sample_id, site_id)) %>% 
  mutate(diff_gm = gm - gm_raw,
         diff_median = median - median_raw, 
         higher_gm = case_when(diff_gm > 0 ~ "higher", TRUE ~ "lower"),
         higher_f = as.integer(as.factor(higher_gm)) - 1, 
         higher_median = case_when(diff_median > 0 ~ "higher", TRUE ~ "lower"),
         higher_f_median = as.integer(as.factor(higher_median)) - 1) %>% 
  ungroup() %>% 
  mutate(p_value_overall = sum(higher_f >0)/nrow(.)) %>% 
  group_by(site_id) %>% 
  mutate(p_value_site = sum(higher_f >0)/max(row_number())) %>% 
  group_by(sample_id) %>% 
  mutate(p_value_sample = sum(higher_f >0)/max(row_number())) %>% 
  ungroup


bayes_pvalue = t_stat_summary %>% 
  ggplot(aes(x = gm, y = gm_raw)) + 
  geom_point(aes(color = higher_gm), size = 0.2, shape = 21) + 
  ggthemes::scale_color_colorblind() +
  geom_abline() + 
  scale_y_log10() + 
  scale_x_log10() +
  annotate(geom = "text", x = 0.02, y = 10, label = paste("Bayesian P = ", 
                                                          round(unique(t_stat_summary$p_value_overall), 2))) +
  coord_cartesian(ylim = c(0.005, 10),
                  xlim = c(0.005, 10)) + 
  labs(x = "Predicted",
       y = "Observed",
       subtitle = "Geometric Mean body size") +
  theme_default() +
  guides(color = "none")

ggview::ggview(bayes_pvalue, width = 5, height = 5)


bayesian_p_median = as.character(round(sum(t_stat_summary$higher_f_median > 0)/nrow(t_stat_summary), 2))

t_stat_summary %>% 
  ggplot(aes(x = median, y = median_raw)) + 
  geom_point(aes(color = higher_median), size = 0.2, shape = 21) + 
  ggthemes::scale_color_colorblind() +
  geom_abline() + 
  scale_y_log10() + 
  scale_x_log10() +
  annotate(geom = "text", x = 0.02, y = 10, label = paste("Bayesian P = ", bayesian_p_median)) +
  coord_cartesian(ylim = c(0.005, 10),
                  xlim = c(0.005, 10)) + 
  labs(x = "Predicted",
       y = "Observed") +
  theme_default()


t_stat_summary %>% 
  group_by(sample_id) %>% 
  reframe(bayes_p = sum(higher_f)/length(row_number())) %>% 
  filter(bayes_p <= 0.025 | bayes_p >= 0.975)

t_stat_summary %>% ungroup() %>% distinct(sample_id, site_id, 
                                          p_value_sample, p_value_overall,
                                          p_value_site) %>% 
  pivot_longer(cols = c(p_value_sample, p_value_site)) %>% 
  mutate(level = str_sub(name, 9, 20),
         higher = case_when(value > 0.975 | value < 0.025 ~ "no", 
                            T ~ "yes"),
         site_no = as.integer(as.factor(site_id)),
         unique = case_when(level == "sample" ~ sample_id,
                            TRUE ~ site_no)) %>% 
  distinct(unique, level, p_value_overall, value, higher) %>% 
  ggplot(aes(x = level, y = value, color = higher)) + 
  geom_jitter(height = 0, width = 0.05) + 
  scale_color_colorblind() + 
  geom_hline(aes(yintercept = p_value_overall),
             linetype = "dotted") + 
  labs(y = "Bayesian p-value",
       color = "Good fit") +
  geom_hline(aes(yintercept = 0.025)) +
  geom_hline(aes(yintercept = 0.975)) +
  annotate(geom = "text", x = 2.3, y = 0.6, label = "Overall p-value") +
  annotate(geom = "text", x = 2.4, y = 0.95, label = "0.975") +
  annotate(geom = "text", x = 2.4, y = 0.05, label = "0.025")



# among sites -------------------------------------------------------------

t_stat_summary %>% 
  ggplot(aes(x = gm, y = gm_raw)) + 
  geom_point(aes(color = higher_gm), size = 0.2, shape = 21) + 
  ggthemes::scale_color_colorblind() +
  geom_abline() + 
  scale_y_log10() + 
  scale_x_log10() +
  # annotate(geom = "text", x = 0.02, y = 10, label = paste("Bayesian P = ", bayesian_p_gm)) +
  # coord_cartesian(ylim = c(0.005, 10),
  # xlim = c(0.005, 10)) + 
  labs(x = "Predicted",
       y = "Observed") +
  theme_default() +
  facet_wrap(~site_id, scales = "free")

# r2 ----------------------------------------------------------------------
# from: Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for Bayesian regression models. The American Statistician.
posts_raw_p = posts_sample_lambdas %>% 
  # filter(sample_id %in% id) %>% 
  filter(.draw <= 100) %>% 
  select(lambda, sample_id, .draw) %>% 
  right_join(dat_resampled %>% distinct(sample_id, dw, no_m2, xmin, xmax, site_id), 
             by = "sample_id", multiple = "all",
             relationship = "many-to-many") 

# 5) sample posterior preds
n_samples = 500
posts_raw_preds = posts_raw_p %>% 
  # filter(.draw <= 4) %>%
  group_by(sample_id) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  ungroup %>% 
  mutate(u = runif(nrow(.), 0, 1)) %>% # uniform draw
  mutate(x = (u*xmax^(lambda+1) +  (1-u) * xmin^(lambda+1) ) ^ (1/(lambda+1))) %>% 
  mutate(data = "y_rep") %>% 
  # mutate(x = round(x,2)) %>% 
  rename(sim = x) %>% 
  group_by(sample_id, .draw) %>% 
  sample_n(n_samples, weight = no_m2, replace = T) %>% 
  ungroup

posts_raw_preds %>% 
  group_by(.draw) %>% 
  reframe(v_fit = var(sim),
          v_pred = var(dw - sim)) %>% 
  mutate(r2 = v_fit/(v_fit + v_pred)) %>% 
  mean_qi(r2)


# confirm r2 - the bayes_R2() function in brms does not work with isd models. The code below checks that our calculation of R2 above is correct. It does this by
# comparing the hand-calculated R2 from a basic stan model to the bayes_R2(). There is an exact match, suggesting that our calculation above is correct.

mod_test = brm(mpg ~ hp, data = mtcars,
               family = "Gaussian")

as_tibble(t(fitted(mod_test, summary = F))) %>%
  mutate(y = mod_test$data$mpg) %>%
  pivot_longer(cols = -y) %>%
  mutate(.draw = parse_number(name)) %>%
  group_by(.draw) %>%
  reframe(v_fit = var(value),
          v_pred = var(y - value)) %>%
  mutate(r2 = v_fit/(v_fit + v_pred)) %>%
  mean_qi(r2)

bayes_R2(mod_test)


