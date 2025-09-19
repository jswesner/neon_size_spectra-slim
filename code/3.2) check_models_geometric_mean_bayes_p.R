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
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")

# 1) Wrangle  -------------------------------------------------------------

# resample raw data
n_pred = 500 # number of individuals to sample. This is 5000 for the ms. But it works ok with 500 as well.

# resampled data ordered by dw with raw "frequency" = n_yx
fitdat_temp = fit_temp_om_gpp$data %>% 
  # filter(sample_id %in% 1:50) %>% 
  group_by(sample_id) %>% 
  sample_n(size = n_pred, replace = T, weight = no_m2) %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number()),
         n = max(row_number()),
         ismax = dw - xmax) 

# get sequential sample ids within site (for plotting)
sequential_sample_ids = fitdat_temp %>%
  ungroup %>%
  distinct(site_id, sample_id) %>%
  arrange(site_id, sample_id) %>%
  group_by(site_id) %>%
  mutate(sample_id_seq = row_number())


# add max and min
fitdat2 = bind_rows(fitdat_temp) %>% 
  left_join(sequential_sample_ids)

# 2) Sample Posterior  -------------------------------------------------------------

# sample yrep from posterior
sample_posts = fit_temp_om_gpp$data %>% 
  select(-dw, -no_m2) %>% 
  distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, ndraws = 300) %>% 
  ungroup %>% 
  distinct(sample_id, .draw, .epred)

site_posts = fit_temp_om_gpp$data %>% 
  select(-dw, -no_m2, -sample_id, -year) %>% 
  distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~(1|site_id), ndraws = 300) %>% 
  ungroup %>% 
  distinct(site_id, .draw, .epred)

# posterior y_reps
post_preds_sample = sample_posts %>% filter(.draw <= 100) %>% 
  left_join(sequential_sample_ids) %>% 
  left_join(fitdat2 %>% ungroup %>% distinct(sample_id, xmin, xmax)) %>%
  expand_grid(individual = 1:n_pred) %>% 
  ungroup  %>% 
  mutate(yrep = rparetocounts(nrow(.), lambda = .epred, xmin = xmin, xmax = xmax)) 

post_preds_site = site_posts %>% filter(.draw <= 100)  %>% 
  left_join(fitdat2 %>% ungroup %>% distinct(sample_id, site_id, xmin, xmax) %>% 
              group_by(site_id) %>% reframe(xmin = median(xmin), xmax = median(xmax))) %>%
  expand_grid(individual = 1:n_pred) %>% 
  ungroup  %>% 
  mutate(yrep = rparetocounts(nrow(.), lambda = .epred, xmin = xmin, xmax = xmax)) 

# id ordered by gm for plotting
reorder_ids = post_preds_sample %>% 
  ungroup %>% 
  group_by(sample_id_seq, site_id) %>% 
  reframe(rank_gm = exp(mean(log(yrep)))) %>% 
  arrange(rank_gm) %>% 
  mutate(sample_id_ordered = 1:nrow(.))

reorder_ids_site = post_preds_site %>% 
  ungroup %>% 
  group_by(site_id) %>% 
  reframe(rank_gm = exp(mean(log(yrep)))) %>% 
  arrange(rank_gm) %>% 
  mutate(site_id_ordered = 1:nrow(.))

# 3) Bayesian P-Values by sample  -------------------------------------------------------------

# calculate summary statistic
sample_gm = fitdat2 %>% 
  group_by(sample_id,sample_id_seq, site_id) %>% 
  reframe(gm_raw = exp(mean(log(dw))),
          median_dw_raw = median(dw)) 

post_gm_sample = post_preds_sample %>% 
  # filter(sample_id_seq == 2) %>%
  group_by(sample_id_seq, sample_id, site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep)))) %>% 
  left_join(reorder_ids) %>%  
  left_join(sample_gm) 

bayes_p_sample = post_gm_sample %>% 
  mutate(diff = gm - gm_raw,
         maxdraws = max(.draw)) %>% 
  group_by(sample_id_seq, sample_id, sample_id_ordered, maxdraws) %>% 
  reframe(sumdiff = sum(diff>0)) %>% 
  mutate(bayes_p = sumdiff/maxdraws,
         outside = case_when(bayes_p > 0.9 | bayes_p < 0.1 ~ "no", TRUE ~ "yes"))

range(bayes_p_sample$bayes_p)
mean(bayes_p_sample$bayes_p)
sd(bayes_p_sample$bayes_p)

# 3b) Bayesian P-Values by site  -------------------------------------------------------------

# calculate summary statistic
site_gm = fitdat2 %>% 
  group_by(site_id) %>% 
  reframe(gm_raw = exp(mean(log(dw))),
          median_dw_raw = median(dw)) 

post_gm_site = post_preds_site %>% 
  group_by(site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep)))) %>% 
  left_join(reorder_ids_site) %>%  
  left_join(site_gm) 

bayes_p_site = post_gm_site %>% 
  mutate(diff = gm - gm_raw,
         maxdraws = max(.draw)) %>% 
  group_by(site_id, maxdraws) %>% 
  reframe(sumdiff = sum(diff>0)) %>% 
  mutate(bayes_p = sumdiff/maxdraws,
         outside = case_when(bayes_p > 0.9 | bayes_p < 0.1 ~ "no", TRUE ~ "yes"))

range(bayes_p_site$bayes_p)
mean(bayes_p_site$bayes_p)
sd(bayes_p_site$bayes_p)

# 4) Plot Summary Statistics - sample  -------------------------------------------------------------

# plot posts vs raw summary statistics
post_pred_gm_summary = post_preds_sample %>% 
  group_by(sample_id_seq, sample_id, site_id, xmin, xmax, .draw) %>% 
  reframe(gm = exp(mean(log(yrep))),
          median_dw = median(yrep)) %>% 
  left_join(reorder_ids) %>% 
  group_by(sample_id, sample_id_seq, sample_id_ordered) 


post_pred_gm = post_gm_sample %>%
  group_by(sample_id, rank_gm, gm_raw) %>% 
  median_qi(gm) %>% 
  left_join(bayes_p_sample) %>% 
  ungroup %>% 
  arrange(gm) %>% 
  mutate(rowname = 1:nrow(.)) %>% 
  ggplot(aes(x = rowname, y = gm)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper, color = "y[rep]"), 
                  alpha = 0.5,
                  size = 0.3) +
  geom_point(data = . %>% ungroup %>% distinct(sample_id, gm_raw, rank_gm, outside, rowname), 
             aes(y = gm_raw, color = "y"), shape = 20, size = 0.12) +
  scale_y_log10() +
  labs(x = "NEON Sample (n = 133)",
       y = "Geometric mean body size (mgDM)") +
  scale_color_manual(values = c("y[rep]" = "#3399ff",
                                "y" = "black"), labels = c(expression(italic(y)),
                                                           expression(italic(y)[rep]))) +
  theme(legend.title = element_blank())

# ggview(post_pred_gm, width = 6.5, height = 5.5)
saveRDS(post_pred_gm, file = "plots/fig_s7_post_pred_gm.rds")
ggsave(post_pred_gm, file = "plots/fig_s7_post_pred_gm.jpg", width = 6.5, height = 5.5)       
       
# 4) Plot Summary Statistics - site  -------------------------------------------------------------

post_gm_site %>%
  group_by(site_id, rank_gm, gm_raw) %>% 
  median_qi(gm) %>% 
  left_join(bayes_p_site) %>% 
  ggplot(aes(x = reorder(site_id, rank_gm), y = gm)) + 
  geom_pointrange(aes(ymin = .lower, ymax = .upper)) +
  geom_point(data = . %>% ungroup %>% distinct(site_id, gm_raw, rank_gm, outside), 
             aes(y = gm_raw, color = outside))

