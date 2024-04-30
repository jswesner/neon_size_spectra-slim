library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

fit_temp_om_gpp$preds = "temp*om*gpp"

# load data
dat_all = readRDS("data/derived_data/dat_all.rds") %>% mutate(temp_mean = mean)

mean_temp = mean(unique(dat_all$temp_mean))
sd_temp = sd(unique(dat_all$temp_mean))
mean_om = mean(unique(dat_all$log_om))
sd_om = sd(unique(dat_all$log_om))
mean_gpp = mean(unique(dat_all$log_gpp))
sd_gpp = sd(unique(dat_all$log_gpp))

# labels
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")

# interaction_plot ---------------------------------
#1) set up data grid/conditions for conditional_effects()
log_gpp_s = quantile(fit_temp_om_gpp$data$log_gpp_s, probs = c(0.25, 0.5, 0.75)) %>% 
  as_tibble() %>% rename(log_gpp_s = value)

#2) Use conditional_effects to get posterior summaries
int_plot = conditional_effects(fit_temp_om_gpp, effects = "mat_s:log_om_s", conditions = log_gpp_s)

#3) Wrangle conditional_effects/backtransform values from scaled to unscaled
int_plot_data = int_plot$`mat_s:log_om_s` %>% as_tibble() %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp,
         gpp = (log_gpp_s*sd_gpp) + mean_gpp,
         om = (log_om_s*sd_om) + mean_om) %>% 
  mutate(quantile_om = case_when(om == min(om) ~ "Low OM",
                                 om == max(om) ~ "High OM",
                                 TRUE ~ "Median OM"),
         quantile_gpp = case_when(gpp == min(gpp) ~ "Low GPP",
                                  gpp == max(gpp) ~ "High GPP",
                                  TRUE ~ "Median GPP")) %>% 
  left_join(facet_gpp) %>% 
  left_join(facet_om)

#4) make interaction plot showing lambda across different values of temp, om, and gpp. Then save it
(interaction_plot = int_plot_data %>% 
    mutate(fill_color = case_when(quantile_om == "Low OM" ~ -0.25,
                                  quantile_om == "High OM" ~ -0.75,
                                  TRUE ~ -0.5)) %>% 
    ggplot(aes(x = mat, y = estimate__)) + 
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.7) + 
    geom_line() +
    facet_grid(facet_gpp~facet_om, labeller = "label_parsed") +
    theme_default() + 
    # scale_fill_viridis() +
    guides(fill = "none",
           color = "none") + 
    labs(y = "\u03bb (ISD exponent)",
         x = "Mean Annual Temperature (\u00b0C)") 
)


ggview::ggview(interaction_plot, width = 4, height = 4)
ggsave(interaction_plot, width = 4, height = 4,
       file = "plots/interaction_plot.jpg")
saveRDS(interaction_plot,  file = "plots/interaction_plot.rds")


# Make univariate plot ---------------------------------------------------------
#1) Get posterior summaries with conditional_effects()
uni_plot = plot(conditional_effects(fit_temp_om_gpp, effects = "mat_s"))

#2) Get individual lambda posterior summaries with add_epred_draws()
sample_dots = fit_temp_om_gpp$data %>% 
  distinct(sample_id, mat_s, log_gpp_s, log_om_s, year, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp)

#3) Make the marginal plot of temperature and lambda
uni_plot_dots = uni_plot$mat_s$data %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  ggplot(aes(x = mat, y = estimate__)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.6) +
  stat_pointinterval(data = sample_dots, aes(y = .epred, group = sample_id), size = 0.05, shape = 1,
                     geom = "pointrange", color = "grey30") +
  theme_default() + 
  labs(y = "\u03bb (ISD exponent)",
       x = "Mean Annual Temperature (\u00b0C)")

saveRDS(uni_plot_dots, file = "plots/uni_plot_dots.rds")
# Make Heat Map plot ------------------------------------------------------

library(brms)
library(tidyverse)
library(tidybayes)
library(ggview)
library(janitor)
library(isdbayes)

#1) load data
dat_all = readRDS("data/derived_data/dat_all.rds") %>% mutate(temp_mean = mean)

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
    scale_fill_viridis_c(direction = -1, na.value="white", option = "F") +
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


# get marginal slopes
marginal_epreds = tibble(mat_s = c(0,1)) %>% 
  # expand_grid(log_gpp_s = seq(min(dat_all$log_gpp_s), max(dat_all$log_gpp_s), length.out = 10)) %>%
  # expand_grid(log_om_s = seq(min(dat_all$log_om_s), max(dat_all$log_om_s), length.out = 30)) %>%
  expand_grid(qlog_om_s) %>% 
  expand_grid(qlog_gpp_s) %>% 
  mutate(no_m2 = 1, xmin = 0.003, xmax = 20000) %>%  # placeholder values. They do not affect the lambda predictions
  add_epred_draws(fit_temp_om_gpp, re_formula = NA)


marginal_epreds %>% 
  ungroup %>% select(mat_s, quantile_om, quantile_gpp, log_gpp_s, log_om_s, .draw, .epred) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(quantile_om, quantile_gpp) %>% 
  median_qi(slope)



# Combine plots -----------------------------------------------------------
# Combine the two plots above into a two-panel plot
library(cowplot)

uni_plot_dots = readRDS("plots/uni_plot_dots.rds") + labs(subtitle = "a)")
interaction_plot = readRDS("plots/interaction_plot.rds") + labs(subtitle = "b)") + scale_x_continuous(breaks = c(0, 11, 22))
isd_heat_plot = readRDS("plots/ms_plots/isd_heat_plot.rds") + labs(subtitle = "c)") 

top = plot_grid(uni_plot_dots, interaction_plot)
temp_twopanel = plot_grid(top, isd_heat_plot, ncol = 1, rel_heights = c(0.7, 0.65))

ggview::ggview(temp_twopanel, width = 6.5, height = 5.5)
ggsave(temp_twopanel, width = 6.5, height = 5.5,
       file = "plots/ms_plots/fig_3_temp_twopanel.jpg", dpi = 500)
saveRDS(temp_twopanel, file = "plots/ms_plots/fig_3_temp_twopanel.rds")


# Combine to 3-panel plot with old Fig. 4
sim_metab_plot = readRDS(file = "plots/ms_plots/sim_metab_plot.rds")

