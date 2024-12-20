library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(ggview)
theme_set(brms::theme_default())

# load these
model_list = readRDS(file = "models/model_list.rds")
fitdat2 = readRDS(file = "data/fitdat2.rds") # raw re-sampled data
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_updated12192024.rds")

# 1) Wrangle  -------------------------------------------------------------

# resample raw data
n_pred = 500 # number of individuals to sample. This is 5000 for the ms. But it works ok with 500 as well.

fitdat_temp = fit_temp_om_gpp$data %>% 
  # filter(sample_id %in% 1:50) %>% 
  group_by(sample_id) %>% 
  sample_n(size = n_pred, replace = T, weight = no_m2) %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number()),
         n = max(row_number()),
         ismax = dw - xmax) 

# get max and min to add
# get_maxmin = fitdat_temp %>% 
#   select(-dw, -no_m2, -n_yx) %>% 
#   distinct() %>% 
#   mutate(dwxmin = xmin,
#          dwxmax = xmax) %>% 
#   pivot_longer(cols = starts_with("dw"),
#                values_to = "dw") %>% 
#   select(-name)

# # get sequential sample ids within site (for plotting)
sequential_sample_ids = fitdat_temp %>%
  ungroup %>%
  distinct(site_id, sample_id) %>%
  arrange(site_id, sample_id) %>%
  group_by(site_id) %>%
  mutate(sample_id_seq = row_number())

# add max and min
fitdat2 = bind_rows(fitdat_temp) %>% 
  left_join(sequential_sample_ids)

saveRDS(fitdat2, file = "data/fitdat2.rds")

# 2) Sample Posterior  -------------------------------------------------------------

# sample yrep from posterior
sample_posts = fit_temp_om_gpp$data %>% 
  select(-dw, -no_m2) %>% 
  distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, ndraws = 300) %>% 
  ungroup %>% 
  distinct(sample_id, .draw, .epred)

post_preds2 = sample_posts %>% filter(.draw <= 100) %>% 
  left_join(sequential_sample_ids) %>% 
  left_join(fitdat2 %>% ungroup %>% distinct(sample_id, xmin, xmax)) %>%
  expand_grid(individual = 1:n_pred) %>% 
  ungroup  %>% 
  mutate(yrep = rparetocounts(nrow(.), lambda = .epred, xmin = xmin, xmax = xmax)) 

saveRDS(post_preds2, file = "posteriors/post_preds2.rds")

# 3) Plot Summary Statistics  -------------------------------------------------------------

# calculate summary statistic
sample_gm = fitdat2 %>% 
  group_by(sample_id,sample_id_seq, site_id) %>% 
  reframe(gm = exp(mean(log(dw))),
          median_dw = median(dw)) 

# plot posts vs raw summary statistics
# randomly select 1 sample from each site
sample_ids = fitdat2 %>% ungroup %>% distinct(sample_id, site_id) %>% 
  group_by(site_id) %>% 
  sample_n(1) %>% 
  pull(sample_id)

reorder_ids = post_preds2 %>% 
  ungroup %>% 
  group_by(sample_id_seq, site_id) %>% 
  reframe(rank_gm = exp(mean(log(yrep)))) %>% 
  arrange(rank_gm) %>% 
  mutate(sample_id_ordered = 1:nrow(.))

post_pred_gm_summary = post_preds2 %>% 
  # filter(sample_id_seq == 2) %>%
  group_by(sample_id_seq, site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep))),
          median_dw = median(yrep)) %>% 
  left_join(reorder_ids)

post_pred_gm = post_pred_gm_summary %>% 
  ggplot(aes(x = sample_id_ordered, y = gm)) +
  stat_pointinterval(alpha = 0.5,
                     size = 0.2,
                     aes(color = "Posterior Predictive", shape = "Posterior Predictive")) +
  scale_y_log10() +
  geom_point(data = sample_gm %>% 
               left_join(reorder_ids), 
             aes(y = gm, color = "Raw Data", shape = "Raw Data"),
             size = 1) +
  guides(fill = "none") +
  scale_color_manual(values = c("#56b4e9", "black")) +
  # scale_shape_manual(values = c(16, 17)) +
  labs(x = "NEON Sample (n = 133)",
       y = "Geometric mean body size (mgDM)") +
  guides(shape = "none") +
  guides(color = guide_legend(override.aes = list(shape = c(16, 17)),
                              title = "")) +
  theme(legend.position = c(0.2, 0.95))

# ggview(post_pred_gm, width = 6.5, height = 5.5)
saveRDS(post_pred_gm, file = "plots/post_pred_gm.rds")
ggsave(post_pred_gm, file = "plots/post_pred_gm.jpg", width = 6.5, height = 5.5)       
       
     
post_pred_median = post_preds2 %>% 
  # filter(sample_id_seq == 2) %>%
  group_by(sample_id_seq, site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep))),
          median_dw = median(yrep)) %>% 
  left_join(reorder_ids) %>% 
  ggplot(aes(x = sample_id_ordered, y = median_dw)) +
  stat_pointinterval(alpha = 0.5, aes(color = "Posterior Predictive", shape = "Posterior Predictive")) +
  scale_y_log10() +
  geom_point(data = sample_gm %>% 
               left_join(reorder_ids), 
             aes(color = "Raw Data", shape = "Raw Data"),
             size = 2.3) +
  guides(fill = "none") +
  scale_color_manual(values = c("#56b4e9", "black")) +
  scale_shape_manual(values = c(16, 17)) +
  labs(x = "NEON Sample (n = 133)",
       y = "Geometric mean body size (mgDM)") 

# 4) Bayesian P-Values  -------------------------------------------------------------

bayes_p = post_preds2 %>% 
  # filter(sample_id_seq == 2) %>%
  group_by(sample_id_seq, site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep)))) %>% 
  left_join(reorder_ids) %>% 
  left_join(sample_gm %>% 
              left_join(reorder_ids) %>% 
              rename(gm_raw = gm)) %>% 
  mutate(diff = gm - gm_raw,
         maxdraws = max(.draw)) %>% 
  group_by(sample_id, maxdraws) %>% 
  reframe(sumdiff = sum(diff>0)) %>% 
  mutate(bayes_p = sumdiff/maxdraws)

range(bayes_p$bayes_p)
mean(bayes_p$bayes_p)
sd(bayes_p$bayes_p)
bayes_p %>% 
  filter(bayes_p > 0.9 | bayes_p <0.1)

