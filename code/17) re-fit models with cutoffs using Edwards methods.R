library(brms)
library(isdbayes)
library(tidyverse)
library(tidybayes)
library(lubridate)

# Instead of re-fitting, you can load the fitted models below

# load models (these are already fit)
fit_temp = readRDS("models/fit_temp.rds")
fit_om = readRDS("models/fit_om.rds")
fit_gpp = readRDS("models/fit_gpp.rds")
fit_temp_om = readRDS("models/fit_temp_om.rds")
fit_temp_gpp = readRDS("models/fit_temp_gpp.rds")
fit_om_gpp = readRDS("models/fit_om_gpp.rds")
fit_temp_om_gpp = readRDS("models/temporary/fit_temp_om_gpp_newxmin.rds")

# load data
neon_cutoff = read_csv("data/NEON_cutoff.csv") %>% mutate(site_id = site)

dat = fit_temp_om_gpp$data
dat_with_cutoff = dat_all %>% left_join(neon_cutoff %>% select(site_id, sum_norm, count_norm, date))

dat_sumnorm = dat_with_cutoff %>% filter(dw >= sum_norm) %>% mutate(xmin = min(dw))
dat_countnorm = dat_with_cutoff %>% filter(dw >= count_norm) %>% mutate(xmin = min(dw))

fit_temp_om_gpp_sumnorm = update(fit_temp_om_gpp, newdata = dat_sumnorm, iter = 1000, chains = 1)
saveRDS(fit_temp_om_gpp_sumnorm, file = "models/fit_temp_om_gpp_sumnorm.rds")

fit_temp_om_gpp_countnorm = update(fit_temp_om_gpp, newdata = dat_countnorm, iter = 100, chains = 1)
saveRDS(fit_temp_om_gpp_countnorm, file = "models/fit_temp_om_gpp_countnorm.rds")


test = dat_sumnorm %>% filter(sample_id == 82)

# fit1 = brm(dw | vreal(no_m2, xmin, xmax) ~ 1, 
#            data = test,
#            stanvars = stanvars,    # required for truncated Pareto
#            family = paretocounts(),# required for truncated Pareto
#            chains = 1, iter = 1000)
# 
# saveRDS(fit1, file = "models/fit1.rds")

fit1 = readRDS(file = "models/fit1.rds")

fit2 = update(fit1, newdata = dat_sumnorm %>% filter(sample_id == 19))
fit2_count = update(fit1, newdata = dat_countnorm %>% filter(sample_id == 19))
fit2_original = update(fit1, newdata = dat %>% filter(sample_id == 19))


fit_newxmin_list = NULL

dat_cutoff_tomodel = dat_with_cutoff %>% mutate(original_cutoff = xmin) %>% 
  pivot_longer(cols = c(original_cutoff, sum_norm, count_norm)) %>% 
  mutate(above_cutoff = case_when(dw >= value ~ "keep")) %>% 
  filter(!is.na(above_cutoff)) %>% 
  group_by(sample_id, name) %>% 
  mutate(xmin = case_when(name == "original_cutoff" ~ xmin,
                          TRUE ~ min(dw))) %>% 
  group_split()

for(i in 1:length(dat_cutoff_tomodel)){
  fit_newxmin_list[[i]] = update(fit1, newdata = dat_cutoff_tomodel[[i]], 
                                 data2 = list(sample_id = unique(dat_cutoff_tomodel[[i]]$sample_id),
                                              cutoff_name = unique(dat_cutoff_tomodel[[i]]$name),
                                              cutoff_value = unique(dat_cutoff_tomodel[[i]]$value),
                                              xmin = unique(dat_cutoff_tomodel[[i]]$xmin),
                                              n = nrow(dat_cutoff_tomodel[[i]])))
}

saveRDS(fit_newxmin_list, file = "models/temporary/fit_newxmin_list.rds")

summary_newxmin_list = NULL

for(i in 1:length(fit_newxmin_list)){
  summary_newxmin_list[[i]] = fit_newxmin_list[[i]]$data %>% 
    distinct(xmin, xmax) %>% 
    mutate(no_m2 = 1) %>% 
    tidybayes::add_epred_draws(fit_newxmin_list[[i]]) %>% 
    median_qi(.epred) %>% 
    mutate(sample_id = fit_newxmin_list[[i]]$data2$sample_id,
           cutoff_name = fit_newxmin_list[[i]]$data2$cutoff_name,
           cutoff_value = fit_newxmin_list[[i]]$data2$cutoff_value,
           xmin = fit_newxmin_list[[i]]$data2$xmin,
           n = fit_newxmin_list[[i]]$data2$n)
}

saveRDS(summary_newxmin_list, file = "posteriors/temporary/summary_newxmin_list.rds")


bind_rows(summary_newxmin_list) %>% 
  left_join(dat %>% distinct(sample_id, mat_s, log_om_s, log_gpp_s)) %>% 
  ggplot(aes(x = mat_s, y = .epred, ymin = .lower, ymax = .upper,
             color = cutoff_name)) + 
  geom_pointrange() +
  geom_smooth(method = lm) +
  facet_wrap(~cutoff_name) +
  NULL




# new xmin
fit_sumnorm_list_xmin = NULL

dat_cutoff_tomodel = dat_with_cutoff %>% mutate(original_cutoff = xmin) %>% 
  pivot_longer(cols = c(original_cutoff, sum_norm, count_norm)) %>% 
  mutate(above_cutoff = case_when(dw >= value ~ "keep")) %>% 
  filter(!is.na(above_cutoff)) %>% 
  group_by(sample_id, name) %>% 
  group_split()

for(i in 1:length(dat_cutoff_tomodel)){
  fit_sumnorm_list_xmin[[i]] = update(fit1, newdata = dat_cutoff_tomodel[[i]] %>% 
                                        mutate(xmin = min(dw)), 
                                 data2 = list(sample_id = unique(dat_cutoff_tomodel[[i]]$sample_id),
                                              cutoff_name = unique(dat_cutoff_tomodel[[i]]$name),
                                              cutoff_value = unique(dat_cutoff_tomodel[[i]]$value),
                                              xmin = unique(dat_cutoff_tomodel[[i]]$xmin),
                                              n = nrow(dat_cutoff_tomodel[[i]])))
}

saveRDS(fit_sumnorm_list_xmin, file = "models/temporary/fit_sumnorm_list_xmin.rds")

summary_sumnorm_list_xmin = NULL

for(i in 1:length(fit_sumnorm_list)){
  summary_sumnorm_list_xmin[[i]] = fit_sumnorm_list_xmin[[i]]$data %>% 
    distinct(xmin, xmax) %>% 
    mutate(no_m2 = 1) %>% 
    tidybayes::add_epred_draws(fit_sumnorm_list_xmin[[i]]) %>% 
    median_qi(.epred) %>% 
    mutate(sample_id = fit_sumnorm_list_xmin[[i]]$data2$sample_id,
           cutoff_name = fit_sumnorm_list_xmin[[i]]$data2$cutoff_name,
           cutoff_value = fit_sumnorm_list_xmin[[i]]$data2$cutoff_value,
           xmin = fit_sumnorm_list_xmin[[i]]$data2$xmin,
           n = fit_sumnorm_list_xmin[[i]]$data2$n)
}

saveRDS(summary_sumnorm_list_xmin, file = "posteriors/temporary/summary_sumnorm_list_xmin.rds")


bind_rows(summary_sumnorm_list_xmin) %>% 
  ggplot(aes(x = sample_id, y = .epred, ymin = .lower, ymax = .upper,
             color = cutoff_name)) + 
  geom_pointrange() +
  NULL

# new xmin
fit_sumnorm_list_xmin = NULL

dat_cutoff_tomodel = dat_with_cutoff %>% mutate(original_cutoff = xmin) %>% 
  pivot_longer(cols = c(original_cutoff, sum_norm, count_norm)) %>% 
  mutate(above_cutoff = case_when(dw >= value ~ "keep")) %>% 
  filter(!is.na(above_cutoff)) %>% 
  group_by(sample_id, name) %>% 
  group_split()

for(i in 1:length(dat_cutoff_tomodel)){
  fit_sumnorm_list_xmin[[i]] = update(fit1, newdata = dat_cutoff_tomodel[[i]] %>% 
                                        mutate(xmin = min(dw)), 
                                      data2 = list(sample_id = unique(dat_cutoff_tomodel[[i]]$sample_id),
                                                   cutoff_name = unique(dat_cutoff_tomodel[[i]]$name),
                                                   cutoff_value = unique(dat_cutoff_tomodel[[i]]$value),
                                                   xmin = unique(dat_cutoff_tomodel[[i]]$xmin),
                                                   n = nrow(dat_cutoff_tomodel[[i]])))
}

saveRDS(fit_sumnorm_list_xmin, file = "models/temporary/fit_sumnorm_list_xmin.rds")

summary_sumnorm_list_xmin = NULL

for(i in 1:length(fit_sumnorm_list)){
  summary_sumnorm_list_xmin[[i]] = fit_sumnorm_list_xmin[[i]]$data %>% 
    distinct(xmin, xmax) %>% 
    mutate(no_m2 = 1) %>% 
    tidybayes::add_epred_draws(fit_sumnorm_list_xmin[[i]]) %>% 
    median_qi(.epred) %>% 
    mutate(sample_id = fit_sumnorm_list_xmin[[i]]$data2$sample_id,
           cutoff_name = fit_sumnorm_list_xmin[[i]]$data2$cutoff_name,
           cutoff_value = fit_sumnorm_list_xmin[[i]]$data2$cutoff_value,
           xmin = fit_sumnorm_list_xmin[[i]]$data2$xmin,
           n = fit_sumnorm_list_xmin[[i]]$data2$n)
}

saveRDS(summary_sumnorm_list_xmin, file = "posteriors/temporary/summary_sumnorm_list_xmin.rds")
 

# 
cutoff_plot = bind_rows(bind_rows(summary_sumnorm_list) %>% mutate(xmin_name = "original_xmin"), 
                        bind_rows(summary_sumnorm_list_xmin) %>% mutate(xmin_name = "new_xmin")) %>% 
  left_join(dat %>% distinct(sample_id, mat_s)) %>%
  mutate(cutoff_name = as.factor(cutoff_name),
         cutoff_name = fct_relevel(cutoff_name, "original_cutoff"),
         xmin_name = fct_relevel(xmin_name, "original_xmin")) %>% 
  ggplot(aes(x = sample_id,
             y = .epred, ymin = .lower, ymax = .upper,
             color = cutoff_name, xmin_name)) + 
  geom_pointrange(size = 0.1, shape = 1) + 
  facet_grid(xmin_name ~ cutoff_name) +
  # geom_smooth(method = lm) +
  ggthemes::scale_color_colorblind() +
  labs(y = "lambda") +
  guides(color = "none") +
  NULL


ggsave(cutoff_plot, file = "plots/temporary/cutoff_plot.jpg", width = 6, height = 6)

bind_rows(dat_cutoff_tomodel) %>% 
  mutate(dw_no_m2 = dw*no_m2) %>% 
  filter(sample_id <= 20) %>% 
  ggplot(aes(x = dw_no_m2, fill = name)) + 
  geom_histogram() +
  scale_x_log10() +
  facet_grid(sample_id~name)




bind_rows(dat_cutoff_tomodel) %>% 
  group_by(sample_id, name) %>% 
  mutate(dw_no_m2 = dw*no_m2) %>% 
  sample_n(500, replace = T, weith = no_m2) %>% 
  arrange(-dw_no_m2) %>% 
  mutate(order = row_number()) %>% 
  filter(sample_id == 50) %>% 
  ggplot(aes(y = order, x = dw_no_m2, color = name)) + 
  geom_point() +
  geom_line(aes(group = sample_id)) +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~name)


# plot --------------------------------------------------------------------


mod_dat = fit_temp_om_gpp_newxmin$data
post_newxmin = fit_temp_om_gpp_newxmin$data %>% 
  select(-dw, -no_m2) %>% 
  distinct() %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp_newxmin, re_formula = NULL, ndraw = 1000)

post_newxmin %>% 
  ggplot(aes(x = mat_s, y = .epred, group = interaction(mat_s, year, site_id))) +
  stat_pointinterval()


post_lines_newxmin = tibble(mat_s = seq(min(mod_dat$mat_s),
                                        max(mod_dat$mat_s),
                                        length.out = 50)) %>% 
  mutate(log_gpp_s = 0,
         log_om_s = 0,
         xmin = min(mod_dat$xmin),
         xmax = max(mod_dat$xmax),
         no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp_newxmin, re_formula = NA)


post_lines_newxmin %>% 
  ggplot(aes(x = mat_s, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.3) +
  stat_pointinterval(data = post_newxmin, aes(group = interaction(mat_s, year, site_id)))


