library(brms)
library(tidyverse)
library(tidybayes)
library(ggview)
library(janitor)
library(isdbayes)

#1) load data
dat_all = readRDS("data/derived_data/dat_all.rds")

mean_temp = mean(unique(dat_all$temp_mean))
sd_temp = sd(unique(dat_all$temp_mean))

# fancy facets
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")

#2) load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

#3) get quantiles for om and gpp
qlog_om_s = quantile(unique(dat_all$log_om_s), probs = c(0.25, 0.5, 0.75), na.rm = T) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), values_to = "log_om_s", names_to = "quantile_om") %>% 
  mutate(quantile_om = c("Low OM", "Median OM", "High OM"))

qlog_gpp_s = quantile(unique(dat_all$log_gpp_s), probs = c(0.25, 0.5, 0.75), na.rm = T) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), values_to = "log_gpp_s", names_to = "quantile_gpp") %>% 
  mutate(quantile_gpp = c("Low GPP", "Median GPP", "High GPP"))

#4) extract posteriors across data grid.
post_lines_heat = tibble(mat_s = seq(min(dat_all$mat_s), max(dat_all$mat_s), length.out = 10)) %>% 
  expand_grid(log_gpp_s = seq(min(dat_all$log_gpp_s), max(dat_all$log_gpp_s), length.out = 10)) %>%
  # expand_grid(log_om_s = seq(min(dat_all$log_om_s), max(dat_all$log_om_s), length.out = 30)) %>%
  expand_grid(qlog_om_s) %>% 
  mutate(no_m2 = 1, xmin = 0.003, xmax = 20000) %>%  # placeholder values. They do not affect the lambda predictions
  add_epred_draws(fit_temp_om_gpp, re_formula = NA) %>% 
  median_qi(.epred) %>%
  mutate(temp_mean = (mat_s*sd_temp) + mean_temp)  %>% 
  mutate(quantile_om = as.factor(quantile_om)) %>%
  mutate(quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM"),
         log_gpp = (log_gpp_s*sd(unique(dat_all$log_gpp))) + mean(unique(dat_all$log_gpp))) %>% 
  left_join(facet_om)

#5) Make and plot
(isd_heat_plot = post_lines_heat %>% 
    ggplot(aes(x = temp_mean, y = log_gpp)) +
    geom_tile(aes(fill = .epred)) +
    facet_wrap(~facet_om, labeller = "label_parsed") +
    scale_fill_viridis_c(direction = -1, na.value="white") +
    geom_point(data = dat_all %>% ungroup %>% distinct(mat_s, log_gpp) %>% 
                 mutate(temp_mean = (mat_s*sd_temp) + mean_temp), shape = 21,
               col = 'black', fill = "white", size = 1.5) +
    theme_default() +
    labs(fill = "\u03bb",
         x = "Mean Annual Temperature (\u00b0C)",
         y = expression(paste("GPP ln(",gC/m ^ 2/yr,")")))+
    theme(legend.key.height= unit(0.4, 'cm'),
          legend.key.width= unit(0.4, 'cm')))

ggview::ggview(isd_heat_plot, width = 6.5, height = 2.2)
ggsave(isd_heat_plot, width = 6, height = 2, file = "plots/ms_plots/isd_heat_plot.jpg", dpi = 500)
saveRDS(isd_heat_plot, file = "plots/ms_plots/isd_heat_plot.rds")
