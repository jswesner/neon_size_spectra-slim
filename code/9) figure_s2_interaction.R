library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggtext)


# Figure S2 ---------------------------------------------------------------

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")
fit_temp_om_gpp$preds = "temp*om*gpp"

# load data
predictors = readRDS("data/predictors_scaled.rds")

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
log_gpp_s_conds = quantile(predictors$log_gpp_s,probs = c(0.25, 0.5, 0.75))

log_om_s_conds = quantile(predictors$log_om_s, probs = c(0.25, 0.5, 0.75))


#2) Use conditional_effects to get posterior summaries
int_plot_conds = tibble(mat_s = seq(min(predictors$mat_s), max(predictors$mat_s),
                               length.out = 30)) %>% 
  expand_grid(log_gpp_s = log_gpp_s_conds) %>% 
  expand_grid(log_om_s = log_om_s_conds) %>% 
  mutate(no_m2 = 1,
         xmin = min(fit_temp_om_gpp$data$xmin),
         xmax = max(fit_temp_om_gpp$data$xmax)) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NA)

#3) Wrangle conditional_effects/backtransform values from scaled to unscaled
int_plot_data = int_plot_conds %>% 
  group_by(mat_s, log_gpp_s, log_om_s) %>% 
  median_qi(.epred) %>% 
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
                                  TRUE ~ -0.5),
           quantile_gpp = fct_relevel(quantile_gpp, "Low GPP", "Median GPP"),
           quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
    ggplot(aes(x = mat, y = .epred)) + 
    geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.7) + 
    geom_line() +
    facet_grid(quantile_gpp~quantile_om) +
    theme_default() + 
    # scale_fill_viridis() +
    guides(fill = "none",
           color = "none") + 
    labs(y = "\u03bb (ISD exponent)",
         x = "Mean Annual Temperature (\u00b0C)") +
    theme(text = element_text(size = 7))
)


ggsave(interaction_plot, width = 4, height = 4,
       file = "plots/ms_plots/interaction_plot.jpg")
saveRDS(interaction_plot,  file = "plots/ms_plots/interaction_plot.rds")



# effect size across om and gpp -------------------------------------------

slope_conds = int_plot_conds = tibble(mat_s = c(0, 1)) %>% 
  expand_grid(log_gpp_s = log_gpp_s_conds) %>% 
  expand_grid(log_om_s = log_om_s_conds) %>% 
  mutate(no_m2 = 1,
         xmin = min(fit_temp_om_gpp$data$xmin),
         xmax = max(fit_temp_om_gpp$data$xmax)) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NA)


slope_conds %>% 
  ungroup %>% 
  select(ends_with("_s"), .draw, .epred) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(log_gpp_s, log_om_s) %>% 
  median_qi(slope)

slope_conds %>% 
  ungroup %>% 
  select(ends_with("_s"), .draw, .epred) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(log_gpp_s, log_om_s)  %>% 
  reframe(prob_pos = sum(slope > 0)/max(.draw))

slope_conds %>% 
  ungroup %>% 
  select(ends_with("_s"), .draw, .epred) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(log_gpp_s, log_om_s) %>% 
  select(log_gpp_s, log_om_s, .draw, slope) %>% 
  mutate(log_gpp_s = round(log_gpp_s, 2),
         log_om_s = round(log_om_s, 2)) %>% 
  pivot_wider(names_from = log_om_s, values_from = slope) %>% 
  mutate(diff = `0.62` - `-0.55`) %>% 
  group_by(log_gpp_s) %>% 
  median_qi(diff)


slope_conds %>% 
  ungroup %>% 
  select(ends_with("_s"), .draw, .epred) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(log_gpp_s, log_om_s) %>% 
  select(log_gpp_s, log_om_s, .draw, slope) %>% 
  mutate(log_gpp_s = round(log_gpp_s, 2),
         log_om_s = round(log_om_s, 2)) %>% 
  pivot_wider(names_from = log_gpp_s, values_from = slope) %>% 
  mutate(diff = `0.56` - `-0.98`) %>% 
  group_by(log_om_s) %>% 
  median_qi(diff)


# interaction by year -----------------------------------------------------

int_plot_data_year = tibble(mat_s = seq(min(fit_temp_om_gpp$data$mat_s), 
                                        max(fit_temp_om_gpp$data$mat_s),
                                        length.out = 20)) %>% 
  expand_grid(log_om_s = c(quantile(unique(fit_temp_om_gpp$data$log_om_s), probs = 0.25)[[1]], 
                           quantile(unique(fit_temp_om_gpp$data$log_om_s), probs = 0.5)[[1]], 
                           quantile(unique(fit_temp_om_gpp$data$log_om_s), probs = 0.75)[[1]]),
              log_gpp_s = c(quantile(unique(fit_temp_om_gpp$data$log_gpp_s), probs = 0.25)[[1]],
                            quantile(unique(fit_temp_om_gpp$data$log_gpp_s), probs = 0.5)[[1]], 
                            quantile(unique(fit_temp_om_gpp$data$log_gpp_s), probs = 0.75)[[1]])) %>% 
  expand_grid(year = unique(fit_temp_om_gpp$data$year)) %>% 
  mutate(site_id = "new",
         no_m2 = 1,
         xmin = min(fit_temp_om_gpp$data$xmin),
         xmax = max(fit_temp_om_gpp$data$xmax)) %>% 
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
  left_join(facet_om) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, allow_new_levels = T,
                  ndraws = 500) 

int_plot_data_summary = int_plot_data_year %>% 
  group_by(mat_s, log_om_s, log_gpp_s, mat, om, gpp, year, quantile_om, quantile_gpp) %>% 
  median_qi(.epred)

int_list = int_plot_data_summary %>% group_by(year) %>% group_split()
interaction_plot_list = NULL

years <- sort(unique(fit_temp_om_gpp$data$year))  # Ensure years are in order

for(i in seq_along(years)){
  interaction_plot_list[[i]] = int_list[[i]] %>%
    mutate(fill_color = case_when(quantile_om == "Low OM" ~ -0.25,
                                  quantile_om == "High OM" ~ -0.75,
                                  TRUE ~ -0.5),
           quantile_gpp = fct_relevel(quantile_gpp, "Low GPP", "Median GPP"),
           quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>%
    ggplot(aes(x = mat, y = .epred)) + 
    geom_ribbon(aes(ymin = .lower,
                    ymax = .upper,
                    group = year),
                alpha = 0.7) +
    geom_line() +
    facet_grid(quantile_gpp ~ quantile_om) +
    theme_default() + 
    guides(fill = "none",
           color = "none") + 
    labs(y = "\u03bb (ISD exponent)",
         x = "Mean Annual Temperature (\u00b0C)",
         subtitle = paste0(letters[i], ") ", years[i])) +
    theme(text = element_text(size = 7)) 
}

library(cowplot)
interaction_by_year = plot_grid(
  interaction_plot_list[[1]] + theme(text = element_text(size = 7), axis.title.x = element_blank(), strip.text.y = element_blank()) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
  interaction_plot_list[[2]] + theme(text = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                     axis.text.y = element_blank(), strip.text.y  =element_text(size = 5)) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
  interaction_plot_list[[3]] + theme(text = element_text(size = 7), axis.title.x = element_blank(), strip.text.y = element_blank(),
                                     strip.text.x = element_blank()) + 
    scale_y_continuous(breaks = c(-3, -2, -1)) ,
  interaction_plot_list[[4]] + theme(text = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                    axis.text.y = element_blank(), strip.text.x = element_blank(), strip.text.y  =element_text(size = 5)) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
  interaction_plot_list[[5]]+ theme(text = element_text(size = 7), axis.title.x = element_blank(), strip.text.y = element_blank(),
                                    strip.text.x = element_blank()) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
  interaction_plot_list[[6]] + theme(text = element_text(size = 7), axis.title.x = element_blank(), axis.title.y = element_blank(),
                                    axis.text.y = element_blank(), strip.text.x = element_blank(), strip.text.y  =element_text(size = 5)) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
  interaction_plot_list[[7]]+ theme(text = element_text(size = 7), 
                                    strip.text.x = element_blank(), strip.text.y  =element_text(size = 5)) + 
    scale_y_continuous(breaks = c(-3, -2, -1)),
          ncol = 2
)

ggsave(interaction_by_year, file = "plots/ms_plots/interaction_by_year.jpg",
       width = 5, height = 8, dpi = 400)
