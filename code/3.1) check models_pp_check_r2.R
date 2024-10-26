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
model_list = readRDS(file = "models/model_list.rds")
mod = model_list$`models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds`

# 2) get data
dat = as_tibble(mod$data)

# 3) extract posteriors 
posts_sample_lambdas = dat %>% 
  distinct(sample_id, .keep_all = T) %>% 
  select(-dw) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(mod, re_formula = NULL) %>% 
  rename(lambda = .epred) %>% 
  ungroup()

saveRDS(posts_sample_lambdas, file = "posteriors/posts_sample_lambdas.rds")

# 4) merge posts and raw data
n_samples = 5000

dat_resampled = dat %>% 
  group_by(sample_id, log_om_s, log_gpp_s, mat_s, site_id, year) %>% 
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
  right_join(dat_resampled %>% 
               group_by(sample_id) %>% 
               slice_sample(n = 500) %>% 
               select(sample_id, dw, no_m2, xmin, xmax, site_id), multiple = "all",
             relationship = "many-to-many") 

# 5) sample posterior preds
sim_posts = posts_raw %>% 
  mutate(x = rparetocounts(lambda = lambda, xmin = xmin, xmax = xmax, n = nrow(.))) %>% 
  mutate(data = "y_rep") %>% 
  # mutate(x = round(x,2)) %>% 
  bind_rows(dat_resampled %>% 
              mutate(data = "y") %>% 
              mutate(.draw = 0,
                     x = dw)) %>% 
  rename(sim = x)

# 6) pick a sample to plot
id = as.integer(runif(15, 1, length(unique(dat$sample_id))))

# 7) Make Plots
# violin
post_pred_si = sim_posts %>%
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
ggsave(post_pred_si, file = "plots/post_pred_i.jpg", width = 6.5, height = 9)

sim_posts %>%
  filter(site_id == "ARIK") %>% 
  # filter(sample_id %in% id) %>%
  ggplot(aes(x = .draw, y = sim, color = data, group = .draw)) + 
  geom_violin() +
  scale_y_log10() +
  facet_wrap(~ sample_id) +
  NULL

# density
sim_posts %>%
  ggplot(aes(x = sim, color = data, group = .draw)) + 
  geom_density() +
  scale_x_log10() +
  facet_wrap(~site_id, scales = "free_y") +
  scale_color_colorblind() +
  NULL

# r2 ----------------------------------------------------------------------
# from: Gelman, A., Goodrich, B., Gabry, J., & Vehtari, A. (2019). R-squared for Bayesian regression models. The American Statistician.
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


# confirm r2 - the bayes_R2() function in brms does not work with isd models. The code below checks that our calculation of R2 above is correct. It does this by
# comparing the hand-calculated R2 from a basic stan model to the bayes_R2(). There is an exact match, suggesting that our calculation above is correct.

# mod_test = brm(mpg ~ hp, data = mtcars,
#                family = "Gaussian")

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
