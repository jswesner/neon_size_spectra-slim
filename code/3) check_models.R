library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(ggh4x)
theme_set(brms::theme_default())


model_list = readRDS(file = "models/model_list.rds")
mod = model_list$`models/fit_temp_newxmin_sumnorm_clauset.rds`
dat_all = readRDS("data/derived_data/dat_all.rds")
dat_info = dat_all %>% select(-dw, -no_m2) %>% distinct()

# chains ------------------------------------------------------------------
rhat_list = list()
for(i in 1:length(model_list)){
  rhat_list[[i]] = tibble(parameter = names(model_list[[i]]$fit),
               rhat = rhat(model_list[[i]])) %>% 
  mutate(model = names(model_list)[i])
}

bind_rows(rhat_list) %>% 
  ggplot(aes(y = model, x = rhat)) + 
  geom_jitter(height = 0.1, width = 0, size = 0.2) +
  xlim(1, 1.1)

# global pp_check -----------------------

dat_sims = mod$data %>% 
  sample_n(2000, weight = no_m2, replace = T)

post_sims = mod$data %>%
  distinct(mat_s, xmin, xmax) %>% 
  mutate(no_m2 = 1) %>% 
  expand_grid(rep = 1:100) %>% 
  add_predicted_draws(mod, 
                      re_formula = NA, 
                      ndraws = 300)

post_sims %>% 
  filter(rep <= 10) %>% 
  ggplot(aes(.prediction)) + 
  geom_density(aes(group = rep)) +
  scale_x_log10() +
  geom_density(data = dat_sims, aes(x = dw), color = "dodgerblue")

post_sims %>% 
  group_by(rep) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  ggplot(aes(x = gm)) + 
  geom_histogram(bins = 50) +
  geom_vline(data = dat_sims %>% ungroup %>% reframe(median = median(dw),
                                                     gm = exp(mean(log(dw)))),
             aes(xintercept = gm))


# stream pp_check ---------------------------

dat_sims_stream = mod$data %>% 
  group_by(site_id) %>% 
  sample_n(2000, weight = no_m2, replace = T)

post_sims_stream = mod$data %>%
  distinct(mat_s, xmin, xmax, site_id) %>% 
  mutate(no_m2 = 1) %>% 
  expand_grid(rep = 1:100) %>% 
  add_predicted_draws(mod, 
                      re_formula = ~ (1|site_id),
                      ndraws = 300)

post_sims_stream %>% 
  filter(rep <= 10) %>% 
  ggplot(aes(x = .prediction)) + 
  geom_density(aes(group = rep)) + 
  facet_wrap(~site_id, scales = "free") + 
  scale_x_log10() +
  geom_density(data = dat_sims_stream, aes(x = dw), color = "dodgerblue")

post_sims_stream %>% 
  group_by(site_id, rep) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  ggplot(aes(x = gm)) + 
  facet_wrap(~site_id, scales = "free_x") +
  geom_histogram(bins = 50) + 
  scale_x_log10() +
  geom_vline(data = dat_sims_stream %>% group_by(site_id) %>% reframe(median = median(dw),
                                                     gm = exp(mean(log(dw)))),
             aes(xintercept = gm))

# sample pp_check ---------------------------
dat_sims_sample = mod$data %>% 
  group_by(sample_id, site_id) %>% 
  sample_n(1000, weight = no_m2, replace = T)

post_sims_sample = mod$data %>%
  distinct(log_om_s, mat_s, xmin, xmax, sample_id, site_id) %>% 
  mutate(no_m2 = 1) %>% 
  expand_grid(rep = 1:100) %>% 
  add_predicted_draws(mod, 
                      re_formula = ~ (1|sample_id), 
                      ndraws = 300)

post_sims_sample %>% 
  filter(site_id == "LEWI") %>% 
  filter(rep <= 10) %>% 
  ggplot(aes(x = .prediction)) + 
  geom_density(aes(group = interaction(sample_id, rep))) + 
  facet_wrap(~site_id) + 
  scale_x_log10() +
  geom_density(data = dat_sims_sample %>% filter(site_id == "LEWI"),
               aes(x = dw, group = sample_id), color = "dodgerblue") +
  facet_wrap(~sample_id)


samples = sample(unique(mod$data$sample_id), 10)
post_sims_sample %>% 
  group_by(sample_id, rep, site_id) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  filter(sample_id %in% samples) %>% 
  ggplot(aes(x = gm, group = sample_id)) + 
  geom_histogram() + 
  facet_wrap(site_id ~ sample_id) +
  scale_x_log10() +
  geom_vline(data = dat_sims_sample %>% filter(sample_id %in% samples) %>% 
               group_by(sample_id, site_id) %>% reframe(gm = exp(mean(log(dw))),
                                               median = median(dw)),
             aes(xintercept = gm))

sample_id_new = tibble(sample_id = sort(unique(mod$data$sample_id))) %>% 
  mutate(sample_id_new = 1:nrow(.))

sample_post_preds = post_sims_sample %>% 
  group_by(sample_id, rep, site_id) %>% 
  reframe(median = median(.prediction),
          gm = exp(mean(log(.prediction)))) %>% 
  left_join(dat_sims_sample %>% 
              group_by(sample_id, site_id) %>% reframe(raw_gm = exp(mean(log(dw))),
                                                       raw_median = median(dw))) %>% 
  left_join(sample_id_new)


sample_post_preds_median = sample_post_preds %>% mutate(diff_gm = raw_gm - gm,
                             diff_median = raw_median - median) %>% 
  group_by(sample_id) %>% 
  mutate(median_diff_gm = median(diff_gm),
         median_diff_median = median(diff_median)) %>% 
  pivot_longer(starts_with("diff")) %>% 
  group_by(sample_id, site_id, sample_id_new, name, median) 

sample_post_preds_median %>% 
  ggplot(aes(x = reorder(sample_id_new, median_diff_gm))) + 
  stat_pointinterval(aes(y = value)) + 
  facet_wrap(~name)

sample_post_preds_median %>% 
  group_by(sample_id, name, median_diff_gm, raw_gm, raw_median) %>% glimpse() %>% 
  median_qi(gm) %>% 
  filter(name == "diff_gm") %>% 
  mutate(cred_includes_true = case_when(raw_gm >= .lower & raw_gm <= .upper ~ 'yes', 
                                        TRUE ~ 'no')) %>% 
  ggplot(aes(x = raw_gm, color = cred_includes_true)) +
  geom_pointrange(aes(y = gm, ymin = .lower, ymax = .upper)) +
  geom_abline() +
  scale_x_log10() +
  scale_y_log10()

sample_post_preds_median %>% 
  group_by(sample_id, name, median_diff_median, raw_gm, raw_median) %>% glimpse() %>% 
  median_qi(median) %>% 
  filter(name == "diff_median") %>% 
  mutate(cred_includes_true = case_when(raw_median >= .lower & raw_median <= .upper ~ 'yes', 
                                        TRUE ~ 'no')) %>% 
  ggplot(aes(x = raw_median, color = cred_includes_true)) +
  geom_pointrange(aes(y = median, ymin = .lower, ymax = .upper)) +
  geom_abline() +
  scale_x_log10() +
  scale_y_log10()


# get lambdas -------------------------------------------------------------

mod_dat = mod$data

post_lambda_lines = tibble(mat_s = seq(min(mod$data$mat_s),
                   max(mod$data$mat_s),
                   length.out = 30)) %>% 
  expand_grid(log_gpp_s = round(quantile(mod$data$log_gpp_s, probs = c(0.25, 0.5, 0.75)), 1),
              log_om_s = round(quantile(mod$data$log_om_s, probs = c(0.25, 0.5, 0.75)), 1),
              xmin = quantile(mod$data$xmin, probs = c(0.25, 0.5, 0.75)),
              xmax = quantile(mod$data$xmax, probs = c(0.25, 0.5, 0.75))
              ) %>% 
  # mutate(xmin = min(mod$data$xmin),
         # xmax = max(mod$data$xmax)) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(mod, re_formula = NA, ndraws = 500)

post_lambda_dots = mod_dat %>% 
  select(-dw, -no_m2) %>% 
  distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(mod, re_formula = NULL, ndraws = 500)

post_lambda_lines %>% 
  group_by(mat_s, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  ggplot(aes(x = mat_s, y = .epred)) + 
  stat_lineribbon(.width = 0.95, alpha = 0.2) + 
  stat_pointinterval(data = post_lambda_dots, aes(group = sample_id))


# play around -------------------------------------------------------------

post_lambda_dots %>% 
  group_by(sample_id, site_id) %>% 
  median_qi(.epred) %>% 
  left_join(dat_info %>% ungroup %>% select(sample_id, year, date)) %>% 
  arrange(.epred) %>% 
  print(n = 20)


post_lambda_dots %>% 
  group_by(site_id) %>% 
  mutate(median = median(.epred)) %>% 
  ggplot(aes(x = .epred, y = reorder(site_id, median),  group = sample_id)) +
  stat_slab()
