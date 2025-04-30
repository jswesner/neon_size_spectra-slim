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
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")
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

#4) make interaction plot showing lambda across different values of temp, om, and gpp. Then save it
(interaction_plot = int_plot_data %>% 
    mutate(fill_color = case_when(quantile_om == "Low OM" ~ -0.25,
                                  quantile_om == "High OM" ~ -0.75,
                                  TRUE ~ -0.5),
           quantile_gpp = fct_relevel(quantile_gpp, "Low GPP", "Median GPP"),
           quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM")) %>% 
    ggplot(aes(x = mat, y = estimate__)) + 
    geom_ribbon(aes(ymin = lower__, ymax = upper__), alpha = 0.7) + 
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


# ggview::ggview(interaction_plot, width = 4, height = 4)
ggsave(interaction_plot, width = 4, height = 4,
       file = "plots/interaction_plot.jpg")
saveRDS(interaction_plot,  file = "plots/interaction_plot.rds")


