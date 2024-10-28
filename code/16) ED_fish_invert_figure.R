library(rstan)
library(tidybayes)
library(brms)
library(tidyverse)
library(isdbayes)
theme_set(theme_default())

mod_all = readRDS("models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")
mod_fish = readRDS("models/fish_or_invert_models/fit_fish_clauset.rds")
mod_inverts = readRDS("models/fish_or_invert_models/fit_inverts_clauset.rds")

mod_all$data2 = list(model = "all")
mod_fish$data2 = list(model = "fish")
mod_inverts$data2 = list(model = "inverts")



mod_list = list(mod_all, mod_fish, mod_inverts)


# get post lines ----------------------------------------------------------

post_lines_list = list()

for(i in 1:length(mod_list)){
  preds = tibble(mat_s = seq(min(mod_list[[i]]$data$mat_s),
                             max(mod_list[[i]]$data$mat_s),
                             length.out = 30)) %>% 
    mutate(log_gpp_s = 0,
           log_om_s = 0,
           xmin = min(mod_list[[i]]$data$xmin),
           xmax = max(mod_list[[i]]$data$xmax)) %>%  
    mutate(no_m2 = 1)
  
  post_lines_list[[i]] = preds %>% add_epred_draws(mod_list[[i]], re_formula = NA) %>% 
    mutate(model = mod_list[[i]]$data2$model)
}

# get post dots -----------------------------------------------------------

post_dots_list = list()

for(i in 1:length(mod_list)){
  preds = mod_list[[i]]$data %>% distinct(sample_id, mat_s, log_gpp_s, log_om_s,
                                          xmin, xmax) %>% 
    mutate(no_m2 = 1)
  
  post_dots_list[[i]] = preds %>% add_epred_draws(mod_list[[i]], re_formula = ~ (1|sample_id)) %>% 
    mutate(model = mod_list[[i]]$data2$model)
}


# plot --------------------------------------------------------------------

dat_all = readRDS("data/derived_data/dat_all.rds")
temp_mean_sd = dat_all %>% distinct(site_id, mean) %>% 
  ungroup %>% reframe(temp_mean = mean(mean),
                      temp_sd = sd(mean))

post_lines_summary = bind_rows(post_lines_list) %>% 
  group_by(model, mat_s) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*temp_mean_sd$temp_sd) + temp_mean_sd$temp_mean) %>% 
  mutate(model = case_when(model == "all" ~ "a) Inverts+Fish",
                           model == "fish" ~ "c) Fish",
                           TRUE ~ "b) Invertebrates"))
  
post_dots_summary = bind_rows(post_dots_list) %>% 
  group_by(model, mat_s, sample_id) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*temp_mean_sd$temp_sd) + temp_mean_sd$temp_mean) %>% 
  mutate(model = case_when(model == "all" ~ "a) Inverts+Fish",
                           model == "fish" ~ "c) Fish",
                           TRUE ~ "b) Invertebrates"))


figure_s2_all_isd = post_dots_summary %>% 
  ggplot(aes(x = mat, y = .epred, ymin = .lower, ymax = .upper)) + 
  facet_wrap(~model) +
  geom_lineribbon(data = post_lines_summary, fill = "grey40", 
                  linewidth = 0.2) + 
  geom_pointrange(shape = ".",
                  linewidth = 0.2,
                  position = position_jitter(width = 0.05)) +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0 C)") +
  guides(fill = "none")


saveRDS(figure_s2_all_isd, file = "plots/ms_plots/figure_s2_all_isd.rds")
ggsave(figure_s2_all_isd, file = "plots/ms_plots/figure_s2_all_isd.jpg", width = 6.5, height = 2.3)










