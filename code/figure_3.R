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
dat_all = readRDS("data/derived_data/dat_all.rds")

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
    ggplot(aes(x = mat, y = estimate__, fill = fill_color)) + 
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.8) + 
    geom_line() +
    facet_grid(facet_gpp~facet_om, labeller = "label_parsed") +
    theme_default() + 
    scale_fill_viridis() +
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

# Combine plots -----------------------------------------------------------
# Combine the two plots above into a two-panel plot
library(patchwork)

uni_plot_dots = readRDS("plots/uni_plot_dots.rds")
interaction_plot = readRDS("plots/interaction_plot.rds")

temp_twopanel = uni_plot_dots / interaction_plot

ggview::ggview(temp_twopanel, width = 4, height = 6.5)
ggsave(temp_twopanel, width = 4, height = 6.5,
       file = "plots/ms_plots/fig_3_temp_twopanel.jpg", dpi = 500)
saveRDS(temp_twopanel, file = "plots/ms_plots/fig_3_temp_twopanel.rds")