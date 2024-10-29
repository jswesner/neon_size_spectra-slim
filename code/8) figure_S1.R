library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
theme_set(brms::theme_default())

#1) load models
model_list = readRDS(file = "models/model_list.rds")
fit_temp = model_list$`models/fit_temp_newxmin_sumnorm_clauset.rds`
fit_om = model_list$`models/fit_om_newxmin_sumnorm_clauset.rds`
fit_gpp = model_list$`models/fit_gpp_newxmin_sumnorm_clauset.rds`
fit_temp_om = model_list$`models/fit_om_temp_newxmin_sumnorm_clauset.rds`
fit_temp_gpp = model_list$`models/fit_temp_gpp_newxmin_sumnorm_clauset.rds`
fit_om_gpp = model_list$`models/fit_om_gpp_newxmin_sumnorm_clauset.rds`
fit_temp_om_gpp = model_list$`models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds`


#2) add a name to the model files (for faceting later), combine to a list
fit_temp$preds = "temp" 
fit_om$preds = "om"
fit_gpp$preds = "gpp"
fit_temp_om$preds = "temp*om" 
fit_temp_gpp$preds = "temp*gpp"
fit_om_gpp$preds = "om*gpp"
fit_temp_om_gpp$preds = "temp*om*gpp"


all_mods = list(fit_temp, 
                fit_om, 
                fit_gpp, 
                fit_temp_om, 
                fit_temp_gpp,
                fit_om_gpp,
                fit_temp_om_gpp)


#3) load data and get means/sds for backtransforming from scaled to unscaled
dat_all = readRDS("data/derived_data/dat_all.rds") %>% mutate(temp_mean = mean)
mean_temp = mean(unique(dat_all$temp_mean))
sd_temp = sd(unique(dat_all$temp_mean))
mean_om = mean(unique(dat_all$log_om))
sd_om = sd(unique(dat_all$log_om))
mean_gpp = mean(unique(dat_all$log_gpp))
sd_gpp = sd(unique(dat_all$log_gpp))


#4) Get parameter posteriors from each mode
# custom function
get_draws_with_preds = function(model = NA){
  as_draws_df(model) %>% 
    mutate(preds = model$preds)
}

# use custom function to get parameter posteriors. Add labels.
all_draws = bind_rows(lapply(all_mods, get_draws_with_preds)) %>% 
  select(starts_with(c("b_", ".draw", "preds"))) %>%
  pivot_longer(cols = starts_with("b_")) %>% 
  filter(!is.na(value)) %>% 
  mutate(temp = case_when(grepl("mat_s", name) ~ "temperature",
                          TRUE ~ "other")) %>% 
  mutate(preds = paste0("~ ", preds, " ..."),
         preds = case_when(preds == "~ temp ..." ~ "a) ~temp ...",
                           preds == "~ gpp ..." ~ "c) ~gpp ...",
                           preds == "~ om ..." ~ "b) ~om ...",
                           preds == "~ temp*gpp ..." ~ "d) ~temp*gpp ...",
                           preds == "~ temp*om ..." ~ "e) ~temp*om ...",
                           preds == "~ om*gpp ..." ~ "f) ~om*gpp ...",
                           TRUE ~ "g) ~temp*om*gpp ...")) %>% 
  mutate(name_length = str_length(name)) %>% 
  mutate(name = str_replace(name, "mat", "temp"),
         name = str_replace(name, "b_", ""))


#5) Make plot
parameter_plot = all_draws  %>% 
  filter(name != "Intercept") %>% 
  ggplot(aes(x = reorder(name, -name_length), y = value, color = temp)) + 
  stat_pointinterval(position = position_dodge(width = 0.4)) +
  geom_hline(yintercept = 0) +
  scale_color_colorblind() + 
  theme_default() +
  labs(x = "Parameter",
       y = "Parameter Value") +
  coord_flip() +
  facet_wrap(~preds, ncol = 3) +
  guides(color = "none") +
  # scale_y_continuous(breaks = c(-0.04, -0.02, 0, 0.02, 0.04)) +
  NULL

# ggview::ggview(parameter_plot, width = 6.5, height = 7)
ggsave(parameter_plot, file = "plots/ms_plots/parameter_plot.jpg", 
       width = 6.5, height = 7, units = "in", dpi = 600 )
saveRDS(parameter_plot, file = "plots/ms_plots/parameter_plot.rds")
