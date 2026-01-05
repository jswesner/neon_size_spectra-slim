library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggtext)
library(ggh4x)


# fish_only ------------------------------------

# load model
fit_temp_om_gpp = readRDS("models/fit_fish.rds")
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


int_data_fish = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Fish")

fish_interaction = int_data_fish %>% 
  ggplot(aes(x = mat, y = estimate__)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.4) + 
  geom_line(aes(group = interaction(quantile_om, quantile_gpp)),
            color = "white") +
  facet_wrap(~facet_omgpp) +
  theme_default() + 
  guides(color = "none",
         fill = "none") +  scale_fill_brewer(palette = "Greens") +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  # geom_point(data = sample_dots_summary, aes(y = .epred),
  #            size = 0.02, shape = ".", position = position_jitter(width = 0.2,
  #                                                                 height = 0)) +
  # geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
  #                                       group = site_id),
  #                 size = 0.03, linewidth = 0.2, shape = 20, position = position_jitter(width = 0.2,
  #                                                                                      height = 0),
  #                 color = 'grey20') +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))

# macros_only ------------------------------------

# load model
fit_temp_om_gpp = readRDS("models/fit_macros.rds")
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


int_data_macros = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Macroinvertebrates only")



fit_macros = readRDS("models/fit_macros.rds")

post_macros = tibble(mat_s = seq(0, 1)) %>% 
  expand_grid(log_gpp_s = c(-1, 0, 1),
              log_om_s = c(-1, 0, 1)) %>% 
  mutate(no_m2 = 1,
         xmin = 0.01, xmax = 1000) %>% 
  add_epred_draws(fit_macros, re_formula = NA)

post_macros %>% 
  ungroup() %>% 
  select(mat_s, log_gpp_s, log_om_s, .draw, .epred) %>% 
  pivot_wider(values_from = .epred, names_from = mat_s) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(log_gpp_s, log_om_s) %>% 
  median_qi(slope)



# first_only --------------------------------------------------------------

# load model
fit_temp_om_gpp = readRDS("models/fit_first.rds")
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

int_data_first = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "First dates only")

first_interaction = int_data_first %>% 
  ggplot(aes(x = mat, y = estimate__)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.4) + 
  geom_line(aes(group = interaction(quantile_om, quantile_gpp)),
            color = "white") +
  facet_wrap(~facet_omgpp) +
  theme_default() + 
  guides(color = "none",
         fill = "none") +  scale_fill_brewer(palette = "Greens") +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  # geom_point(data = sample_dots_summary, aes(y = .epred),
  #            size = 0.02, shape = ".", position = position_jitter(width = 0.2,
  #                                                                 height = 0)) +
  # geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
  #                                       group = site_id),
  #                 size = 0.03, linewidth = 0.2, shape = 20, position = position_jitter(width = 0.2,
  #                                                                                      height = 0),
  #                 color = 'grey20') +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))

# last_only ------------------------------------
# load model
fit_temp_om_gpp = readRDS("models/fit_last.rds")
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

int_data_last = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Last dates only")

last_interaction = int_data_last %>% 
  ggplot(aes(x = mat, y = estimate__)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.4) + 
  geom_line(aes(group = interaction(quantile_om, quantile_gpp)),
            color = "white") +
  facet_wrap(~facet_omgpp) +
  theme_default() + 
  guides(color = "none",
         fill = "none") +  scale_fill_brewer(palette = "Greens") +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  # geom_point(data = sample_dots_summary, aes(y = .epred),
  #            size = 0.02, shape = ".", position = position_jitter(width = 0.2,
  #                                                                 height = 0)) +
  # geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
  #                                       group = site_id),
  #                 size = 0.03, linewidth = 0.2, shape = 20, position = position_jitter(width = 0.2,
  #                                                                                      height = 0),
  #                 color = 'grey20') +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))


# spatial gp model --------------------------------------------------------------
# load model
fit_temp_om_gpp = readRDS("models/fit_latlong.rds")
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


int_data_gp = int_plot_data %>% 
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

int_data_gp = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Spatial")


# surber ------------------------------------------------------------------
# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_surber.rds")
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


int_data_surber = int_plot_data %>% 
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

int_data_surber = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Surber")

# full model --------------------------------------------------------------
# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")
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

int_data_full = int_plot_data %>% 
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
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
  mutate(model = "Full model")

full_interaction = int_data_full %>% 
  ggplot(aes(x = mat, y = estimate__)) + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.4) + 
  geom_line(aes(group = interaction(quantile_om, quantile_gpp)),
            color = "white") +
  facet_wrap(~facet_omgpp) +
  theme_default() + 
  guides(color = "none",
         fill = "none") +  scale_fill_brewer(palette = "Greens") +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  # geom_point(data = sample_dots_summary, aes(y = .epred),
  #            size = 0.02, shape = ".", position = position_jitter(width = 0.2,
  #                                                                 height = 0)) +
  # geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
  #                                       group = site_id),
  #                 size = 0.03, linewidth = 0.2, shape = 20, position = position_jitter(width = 0.2,
  #                                                                                      height = 0),
  #                 color = 'grey20') +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))



# plot --------------------------------------------------------------------

macro_first_last = bind_rows(int_data_macros,
                             int_data_first,
                             int_data_last,
                             int_data_fish, 
                             int_data_gp,
                             int_data_surber) %>% 
  mutate(model = as.factor(model),
         model = fct_relevel(model, "First dates only",
                             "Last dates only", "Macroinvertebrates only", "Fish", "Surber"))

full_formatted = int_data_full %>% select(-model) %>% 
  expand_grid(model = unique(macro_first_last$model))

macro_first_last_plot = macro_first_last %>% 
  ggplot(aes(x = mat, y = estimate__, fill = model, group = model)) + 
  geom_ribbon(data = full_formatted,
              aes(ymin = lower__, ymax = upper__), alpha = 0.9,
              fill = "black") + 
  geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.6) + 
  geom_line(aes(group = interaction(model,quantile_om, quantile_gpp)),
            color = "white") +
  facet_grid2(model~facet_omgpp, scales = "free_y") + 
  theme_default() + 
  guides(color = "none",
         fill = "none") +  
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  theme(strip.text = element_text(hjust = 0, size = 8)) + 
  coord_cartesian(ylim = c(-3.5, -0.5))

ggsave(macro_first_last_plot, file = "plots/macro_first_last_plot.jpg",
       width = 6, height = 8, dpi = 400)
