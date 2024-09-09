library(brms)
library(isdbayes)
library(tidyverse)
library(tidybayes)
library(lubridate)

lambda = c(-2, -1.7, -1.5)
xmin = c(0.001, 0.01, 0.1, 1)
xmax = 1000
counts = 1
cull_prop = c(1, 1.5, 3)

dat = tibble(lambda = lambda,
             # xmin = xmin,
             xmax = xmax,
             counts = counts) %>% 
  expand_grid(xmin) %>% 
  expand_grid(cull_prop = cull_prop) %>% 
  expand_grid(xmin_choice = c("original", "new")) %>% 
  mutate(cull = xmin*cull_prop)

dat_list = dat %>% group_by(lambda, xmin, xmax, counts, cull, xmin_choice) %>% group_split()

d_list = NULL

for(i in 1:length(dat_list)){
  d_list[[i]] = tibble(dw = rparetocounts(n = 300, lambda = dat_list[[i]]$lambda, xmin = dat_list[[i]]$xmin,
              xmax = dat_list[[i]]$xmax)) %>% 
  mutate(xmin = dat_list[[i]]$xmin, xmax = dat_list[[i]]$xmax, counts = dat_list[[i]]$counts,
         cull = dat_list[[i]]$cull, xmin_choice = dat_list[[i]]$xmin_choice,
         true_lambda = dat_list[[i]]$lambda, cull_prop = dat_list[[i]]$cull_prop) %>% 
  filter(dw >= cull) %>% 
    mutate(xmin = case_when(xmin_choice == "original" ~ xmin, TRUE ~ min(dw)),
           culled = case_when(xmin == cull ~ "not culled", TRUE ~ "culled"),
           data_cull_xmin = paste0(culled, "_", xmin_choice),
           no_m2 = counts)
}


brm_dummy <- readRDS("models/fit1.rds")

brm_xmins = NULL

for(i in 1:length(d_list)){
  brm_xmins[[i]] = update(brm_dummy, newdata = d_list[[i]],
                          data2 = list(data = unique(d_list[[i]]$data_cull_xmin),
                                       true_lambda = unique(d_list[[i]]$true_lambda),
                                       cull = unique(dat_list[[i]]$cull)),
                          iter = 1000)
}

saveRDS(brm_xmins, file = "models/temporary/brm_xmins.rds")


brm_xmins = readRDS(file = "models/temporary/brm_xmins.rds")

brm_xmins_post = NULL

for(i in 1:length(brm_xmins)){
  brm_xmins_post[[i]] = as_draws_df(brm_xmins[[i]]) %>% 
    mutate(data = unique(brm_xmins[[i]]$data2$data),
           true_lambda = unique(brm_xmins[[i]]$data2$true_lambda),
           cull = unique(brm_xmins[[i]]$data2$cull),
           fit = i)
}

bind_rows(brm_xmins_post) %>%
  mutate(error = b_Intercept - true_lambda) %>% 
  mutate(shift = case_when(data == "not culled_original" ~ -0.02,
                           data == "culled_new" ~ 0,
                           TRUE ~ 0.02),
         true_lambda = true_lambda + shift,
         data = as.factor(data),
         data = fct_relevel(data, "not culled_original", "culled_new")) %>% 
  ggplot(aes(x = true_lambda, y = b_Intercept, color = data, group = interaction(data, cull, i, true_lambda))) +
  stat_pointinterval() +
  ggthemes::scale_color_colorblind() +
  # geom_hline(yintercept = 0) +
  # facet_wrap(~cull) +
  geom_abline() +
  labs(y = "\u03bb")
