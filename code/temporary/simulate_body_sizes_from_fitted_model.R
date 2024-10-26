library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(ggh4x)
theme_set(brms::theme_default())

# use the posterior lambdas to simulate body sizes at each site/sample/etc. Use the global xmin and site-level xmax values to do this instead of
# the culled values.

mod = readRDS("models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")
dat_all = readRDS("data/derived_data/dat_all.rds")
dat_minmax = dat_all %>% 
  group_by(site_id) %>%
  # ungroup %>% 
  mutate(xmax = max(dw)) %>%  # site xmax
  ungroup %>%
  mutate(xmin = min(dw)) %>%  # global xmin
  select(-dw, -no_m2) %>% distinct() %>% 
  left_join(mod$data %>% distinct(sample_id, xmin, xmax) %>% rename(xmin_clauset = xmin, xmax_clauset = xmax))


# body sizes per site -----------------------------------------------------

posts_site = dat_minmax %>%
  ungroup %>% 
  distinct(site_id, xmin, xmax, log_om_s, log_gpp_s, mat_s) %>% 
  mutate(no_m2 = 1)

post_individuals = posts_site %>% 
  add_predicted_draws(mod, re_formula = ~ (1|site_id), ndraws = 4000)

site_lambdas = posts_site %>% 
  add_epred_draws(mod, re_formula = ~ (1|site_id)) %>% 
  group_by(site_id) %>% 
  median_qi(.epred)

post_individuals %>% 
  left_join(site_lambdas) %>% 
  ggplot(aes(x = reorder(site_id, .epred), y = .prediction)) + 
  geom_jitter(aes(size = .prediction),
              width = 0.2,
              alpha = 0.06) +
  # geom_point(shape = 1, size = 0.2) +
  scale_y_log10() +
  scale_size_continuous(breaks = c(0.01, 0.1, 1, 10, 100, 1000, 
                                   1e3, 1e4)) +
  geom_boxplot(aes(group = .epred), outlier.shape = NA) +
  # stat_halfeye(alpha = 0.7) +
  guides(size = "none") +
  coord_flip() +
  labs(x = "Streams (ordered by \u3bb)",
       y = "Individual dry mass (mg)") +
  NULL

# same graph but with raw data
dat_all_sims = dat_all %>% 
  group_by(site_id) %>% 
  sample_n(size = 4000, weight = no_m2, replace = T) %>% 
  mutate(.prediction = dw,
         data = "raw") 

post_individuals %>% 
  left_join(site_lambdas) %>% 
  mutate(data = "simulated") %>% 
  bind_rows(dat_all_sims %>% left_join(site_lambdas)) %>% 
  ggplot(aes(x = reorder(site_id, .epred), y = .prediction)) + 
  geom_jitter(aes(size = .prediction),
              width = 0.2,
              alpha = 0.06) +
  # geom_point(shape = 1, size = 0.2) +
  scale_y_log10() +
  scale_size_continuous(breaks = c(0.01, 0.1, 1, 10, 100, 1000, 
                                   1e3, 1e4)) +
  geom_boxplot(aes(group = .epred), outlier.shape = NA) +
  # stat_halfeye(alpha = 0.7) +
  guides(size = "none") +
  facet_wrap(~data) +
  coord_flip() +
  labs(x = "Streams (ordered by \u3bb)",
       y = "Individual dry mass (mg)") +
  NULL


# create walk through with one site -------------------------
posts_ind_clauset_xmins = dat_minmax %>%
  ungroup %>% 
  distinct(site_id, xmin_clauset, xmax_clauset, log_om_s, log_gpp_s, mat_s) %>% 
  mutate(no_m2 = 1, 
         data = "simulated_clauset") %>% 
  left_join(site_lambdas) %>% 
  mutate(xmin = xmin_clauset, xmax = xmax_clauset) %>% # changes xmin/xmax to be the "culled" versions 
  add_predicted_draws(mod, re_formula = ~ (1|site_id), ndraws = 4000)

post_ind_allsteps = post_individuals %>% 
  left_join(site_lambdas) %>% 
  mutate(data = "e) simulated from ISD, empirical xmin") %>% 
  bind_rows(dat_all_sims %>% left_join(site_lambdas) %>% mutate(data = "b) raw corrected")) %>% 
  bind_rows(dat_all %>% left_join(site_lambdas) %>% mutate(data = "a) raw uncorrected",
                                                           .prediction = dw)) %>% 
  bind_rows(posts_ind_clauset_xmins %>% mutate(data = "d) simulated from ISD, wrong xmin")) %>% 
  bind_rows(mod$data %>% group_by(site_id) %>% sample_n(4000, replace = T, weight = no_m2) %>% 
              mutate(.prediction = dw, 
                     data = "c) raw corrected and culled"))

post_ind_allsteps_summarized = 
  post_ind_allsteps %>% 
  group_by(site_id, data) %>% 
  reframe(arithmetic_mean = mean(.prediction),
          geometric_mean = exp(mean(log(.prediction))),
          median = median(.prediction)) %>% 
  pivot_longer(cols = c(arithmetic_mean, geometric_mean, median),
               names_to = "summary", values_to = "summary_value") %>% 
  filter(site_id == "WALK") 

shape_label = post_ind_allsteps_summarized %>% filter(data == "a) raw uncorrected") %>% 
  mutate(label = c("Arithmetic Mean" , "Geometric Mean", "Median"))

plot_walkthrough = post_ind_allsteps %>% 
  filter(site_id == "WALK") %>%  
  ggplot(aes(y = reorder(data, desc(data)), x = .prediction)) +
  geom_jitter(width = 0, height = 0.01, alpha = 0.1, shape = 1, size = 0.2) +
  scale_x_log10() +
  geom_point(data = post_ind_allsteps_summarized,
             aes(color = summary,
                 shape = summary, 
                 x = summary_value),
             size = 3) +
  scale_color_brewer(type = "qual", palette = 3) +
  # scale_color_colorblind() +
  theme(legend.position = "top",
        legend.text = element_text(size = 8)) +
  labs(color = "",
       y = "",
       x = "mg dry mass",
       subtitle = "What is the 'mean' body size?") +
  ggrepel::geom_text_repel(data = shape_label, aes(label = label, x = summary_value),
                           size = 2.6, nudge_y = 0.5) +
  guides(color = "none",
         shape = "none")

ggsave(plot_walkthrough, file = "plots/plot_walkthrough.jpg", width = 6.5, height = 4)


# compare medians
simulated_medians = post_individuals %>% 
  left_join(site_lambdas) %>%
  group_by(site_id) %>% 
  reframe(sim_lower = quantile(.prediction, probs = 0.025),
          sim_median = median(.prediction),
          sim_gm = exp(mean(log(.prediction))),
          sim_upper = quantile(.prediction, probs = 0.975)) %>% 
  mutate(data = "simulated")

raw_medians = fit_temp_om_gpp_newxmin_sumnorm$data %>% 
  group_by(site_id) %>% 
  sample_n(4000, weight = no_m2, replace = T) %>% 
  left_join(site_lambdas) %>%
  group_by(site_id) %>% 
  reframe(raw_lower = quantile(dw, probs = 0.025),
          raw_median = median(dw),
          raw_gm = exp(mean(log(dw))),
          raw_upper = quantile(dw, probs = 0.975)) %>% 
  mutate(data = "raw")

simulated_medians %>% left_join(raw_medians %>% select(raw_lower, raw_upper, raw_median, site_id)) %>% 
  ggplot(aes(x = raw_median, y = sim_median)) + 
  geom_pointrange(aes(xmin = raw_lower, xmax = raw_upper)) +
  geom_pointrange(aes(ymin = sim_lower, ymax = sim_upper)) +
  scale_x_log10() +
  scale_y_log10()

simulated_medians %>% left_join(raw_medians %>% select(raw_lower, raw_upper, raw_gm,
                                                       raw_median, site_id)) %>% 
  ggplot(aes(x = raw_gm, y = sim_gm)) + 
  geom_point() +
  scale_x_log10() +
  scale_y_log10()


simulated_medians %>% 
  arrange(-sim_gm) %>% 
  print(n = Inf)

simulated_medians %>% 
  ggplot(aes(x = sim_gm, y = sim_median)) + 
  geom_point()

raw_medians %>% 
  arrange(-raw_gm) %>% 
  print(n = Inf)

raw_medians %>% 
  ggplot(aes(x = raw_gm, y = raw_median)) + 
  geom_point()
