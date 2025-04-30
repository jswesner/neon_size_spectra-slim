library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggtext)

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")
fit_temp_om_gpp$preds = "temp*om*gpp"

# load data
predictors = readRDS("data/predictors_scaled.rds") %>% 
  mutate(log_om = log(om))

mean_temp = attributes(predictors$mat_s)$`scaled:center`
sd_temp = attributes(predictors$mat_s)$`scaled:scale`
mean_om = attributes(predictors$log_om_s)$`scaled:center`
sd_om = attributes(predictors$log_om_s)$`scaled:scale`
mean_gpp = attributes(predictors$log_gpp_s)$`scaled:center`
sd_gpp = attributes(predictors$log_gpp_s)$`scaled:scale`



# labels
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")

# interaction_plot ---------------------------------
#1) set up data grid/conditions for conditional_effects()
log_gpp_s = tibble(log_gpp_s = c(-1, 0, 1))

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

# regression with dots and interaction ------------------------------------

int_data = int_plot_data %>% 
  mutate(fill_color = case_when(quantile_om == "Low OM" ~ -0.25,
                                quantile_om == "High OM" ~ -0.75,
                                TRUE ~ -0.5),
         facet_gpp = fct_relevel(facet_gpp, "paste(GPP[0.75])", "paste(GPP[0.5])")) %>% 
  mutate(om_gpp = paste(quantile_om, quantile_gpp, sep = "|")) %>% 
  filter(om_gpp %in% c("Low OM|Low GPP", "Median OM|Median GPP", "High OM|High GPP")) %>%
  mutate(facet_omgpp = case_when(grepl("Low", om_gpp) ~ "a) Low GPP|Low OM",
                                 grepl("Median", om_gpp) ~ "b) Median GPP|Median OM",
                                 TRUE ~ "c) High GPP|High OM")) %>% 
  mutate(om_gpp = as_factor(om_gpp),
         om_gpp = fct_relevel(om_gpp, "Low OM|Low GPP", "Median OM|Median GPP"),
         quantile_om = as.factor(quantile_om),
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM"))

#2) Get individual lambda posterior summaries with add_epred_draws()
sample_dots = fit_temp_om_gpp$data %>% 
  distinct(sample_id, mat_s, log_gpp_s, log_om_s, year, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp)

sample_dots_summary = sample_dots %>% 
  group_by(sample_id, mat_s) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) 

site_dots = fit_temp_om_gpp$data %>% 
  distinct(mat_s, log_gpp_s, log_om_s, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~ (1|site_id)) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  group_by(mat, site_id) %>% 
  median_qi(.epred)

interaction_with_data = int_data %>% 
  ggplot(aes(x = mat, y = estimate__)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.4) + 
  geom_line(aes(group = interaction(quantile_om, quantile_gpp)),
            color = "white") +
  facet_wrap(~facet_omgpp) +
  theme_default() + 
  guides(color = "none",
         fill = "none") +
  # scale_fill_gradientn(colors = RColorBrewer::brewer.pal(9, "Greens")[3:9]) +
  scale_fill_brewer(palette = "Greens") +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  geom_point(data = sample_dots_summary, aes(y = .epred),
             size = 0.02, shape = ".", position = position_jitter(width = 0.2,
                                                                height = 0)) +
  geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
                                        group = site_id),
                  size = 0.03, linewidth = 0.2, shape = 20, position = position_jitter(width = 0.2,
                                                                                      height = 0),
                  color = 'grey20') +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))

ggsave(interaction_with_data, file = "plots/interaction_with_data.jpg",
       width = 6.5, height = 2.3)
saveRDS(interaction_with_data, file = "plots/interaction_with_data.rds")

# Make Heat Map plot ------------------------------------------------------
# fancy facets
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")

#2) load model
fit_temp_om_gpp

#3) get quantiles for om and gpp
qlog_om_s = quantile(unique(predictors$log_om_s), probs = c(0.25, 0.5, 0.75), na.rm = T) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), values_to = "log_om_s", names_to = "quantile_om") %>% 
  mutate(quantile_om = c("Low OM", "Median OM", "High OM"))

qlog_gpp_s = quantile(unique(predictors$log_gpp_s), probs = c(0.25, 0.5, 0.75), na.rm = T) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), values_to = "log_gpp_s", names_to = "quantile_gpp") %>% 
  mutate(quantile_gpp = c("Low GPP", "Median GPP", "High GPP"))

qlog_temp_s = quantile(unique(predictors$mat_s), probs = c(0.25, 0.5, 0.75), na.rm = T) %>% 
  as_tibble() %>% 
  pivot_longer(cols = everything(), values_to = "log_temp_s", names_to = "quantile_temp") %>% 
  mutate(quantile_temp = c("Low temp", "Median temp", "High temp"))

site_cutoffs = predictors %>% distinct(site_id, log_gpp_s, log_om_s, mat_s) %>% 
  mutate(low_gpp = min(qlog_gpp_s$log_gpp_s),
         high_gpp = max(qlog_gpp_s$log_gpp_s)) %>% 
  mutate(quantile_gpp = case_when(log_gpp_s <= low_gpp ~ "Low GPP",
                                  log_gpp_s >= high_gpp ~ "High GPP", 
                                  TRUE ~ "Median GPP")) %>% 
  mutate(low_om = min(qlog_om_s$log_om_s),
         high_om = max(qlog_om_s$log_om_s)) %>% 
  mutate(quantile_om = case_when(log_om_s <= low_om ~ "Low OM",
                                  log_om_s >= high_om ~ "High OM", 
                                  TRUE ~ "Median OM"))  %>% 
  mutate(low_temp = min(qlog_temp_s$log_temp_s),
         high_temp = max(qlog_temp_s$log_temp_s)) %>% 
  mutate(quantile_temp = case_when(mat_s <= low_temp ~ "Low temp",
                                 mat_s >= high_temp ~ "High temp", 
                                 TRUE ~ "Median temp")) 

#4) extract posteriors across data grid.
post_lines_heat = tibble(mat_s = seq(min(predictors$mat_s), max(predictors$mat_s), length.out = 10)) %>% 
  # expand_grid(log_gpp_s = seq(min(predictors$log_gpp_s), max(predictors$log_gpp_s), length.out = 10)) %>%
  expand_grid(log_om_s = seq(min(predictors$log_om_s), max(predictors$log_om_s), length.out = 30)) %>%
  expand_grid(qlog_gpp_s) %>%
  mutate(no_m2 = 1, xmin = 0.003, xmax = 20000) %>%  # placeholder values. They do not affect the lambda predictions
  add_epred_draws(fit_temp_om_gpp, re_formula = NA) %>% 
  median_qi(.epred) %>%
  mutate(temp_mean = (mat_s*sd_temp) + mean_temp)  %>% 
  mutate(quantile_gpp = as.factor(quantile_gpp)) %>%
  mutate(quantile_gpp = fct_relevel(quantile_gpp, "Low GPP", "Median GPP"),
         log_om = (log_om_s*sd_om) + mean_om) %>%
  left_join(facet_gpp) %>% 
  mutate(panels = case_when(quantile_gpp == "Low GPP" ~ "d) Low GPP",
                            quantile_gpp == "Median GPP" ~ "e) Median GPP",
                            quantile_gpp == "High GPP" ~ "f) High GPP"))

#5) Make and plot
(isd_heat_plot = post_lines_heat %>% 
    # filter(quantile_gpp == "Median GPP") %>% 
    ggplot(aes(x = temp_mean, y = log_om)) +
    geom_tile(aes(fill = .epred)) +
    facet_wrap(~panels) +
    scale_fill_viridis_c(direction = -1, na.value="white", option = "D") +
    geom_point(data = predictors %>% ungroup %>% 
                 distinct(mat_s, log_om, site_id) %>% 
                 mutate(temp_mean = (mat_s*sd_temp) + mean_temp) , shape = 21,
               col = 'black', fill = "white", size = 1.5) +
    theme_default() +
    labs(fill = "\u03bb",
         x = "Mean Annual Temperature (\u00b0C)",
         # y = expression(paste("GPP ln(",gC/m ^ 2/yr,")")),
         y = "OM ln(AFDM/m<sup>2</sup>)") +
    theme(legend.key.height= unit(0.4, 'cm'),
          legend.key.width= unit(0.4, 'cm'),
          axis.title.y = ggtext::element_markdown()))

ggsave(isd_heat_plot, width = 6, height = 2, file = "plots/ms_plots/isd_heat_plot.jpg", dpi = 500)
saveRDS(isd_heat_plot, file = "plots/ms_plots/isd_heat_plot.rds")


# get marginal slopes
marginal_epreds = tibble(mat_s = c(0,1)) %>% 
  # expand_grid(log_gpp_s = seq(min(predictors$log_gpp_s), max(predictors$log_gpp_s), length.out = 10)) %>%
  # expand_grid(log_om_s = seq(min(predictors$log_om_s), max(predictors$log_om_s), length.out = 30)) %>%
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

# uni_plot_dots = readRDS("plots/uni_plot_dots_color_mat.rds") + labs(subtitle = "a)")
# interaction_plot = readRDS("plots/interaction_plot.rds") + labs(subtitle = "a)") + scale_x_continuous(breaks = c(0, 11, 22))
interaction_with_data = readRDS(file = "plots/interaction_with_data.rds") + 
  theme(strip.text = element_text(size = 9),
        text = element_text(size = 9))

(isd_heat_plot = readRDS("plots/ms_plots/isd_heat_plot.rds") + 
    theme(
    # legend.position = "bottom",
          # legend.direction = "horizontal",
          legend.text = element_text(size = 7),
          strip.text = element_text(size = 9, hjust = 0),
          text = element_text(size = 9)
    ) +
    NULL
)

# top = plot_grid(uni_plot_dots, interaction_plot)
(temp_twopanel = plot_grid(interaction_with_data, isd_heat_plot, ncol = 1,
                           align = "v"))

# ggview::ggview(temp_twopanel, width = 6.5, height = 5.5)
ggsave(temp_twopanel, width = 6.5, height = 4.2,
       file = "plots/ms_plots/fig_3_temp_twopanel_color_mat.jpg", dpi = 500)
saveRDS(temp_twopanel, file = "plots/ms_plots/fig_3_temp_twopanel.rds")


