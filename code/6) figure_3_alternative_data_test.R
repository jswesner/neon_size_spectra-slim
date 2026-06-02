# remake Fig 3 with alternate data, like only fish, only macors, only first dates, etc.
# See how much the patterns change relative to the main fig. 3
library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggtext)
library(ggh4x)

# load alternative models ---------------------------------
fit_fish = readRDS("models/fit_fish.rds")
fit_fish$model = "Fish"
fit_fish$preds = "temp*om*gpp"

fit_macros = readRDS("models/fit_macros.rds")
fit_macros$model = "Macroinvertebrates"
fit_macros$preds = "temp*om*gpp"

fit_first = readRDS("models/fit_first.rds")
fit_first$model = "First dates only"
fit_first$preds = "temp*om*gpp"

fit_last = readRDS("models/fit_last.rds")
fit_last$model = "Last dates only"
fit_last$preds = "temp*om*gpp"

fit_latlong = readRDS("models/fit_latlong.rds")
fit_latlong$model = "Spatial"
fit_latlong$preds = "temp*om*gpp"

fit_surber = readRDS("models/fit_temp_om_gpp_surber.rds")
fit_surber$model = "Surber"
fit_surber$preds = "temp*om*gpp"


fit_list = list(fit_fish, fit_macros, fit_first, fit_last, fit_latlong, fit_surber)


# load wrangling data -----------------------------------------------------

lat_long = read_csv("data/classified_coordinates.csv") %>% clean_names() %>% distinct(site_id, lat, long)

predictors = read_csv("data/predictors.csv") %>% 
  mutate(log_om = log(om)) %>% 
  left_join(lat_long)

mean_temp = attributes(predictors$mat_s)$`scaled:center`
sd_temp = attributes(predictors$mat_s)$`scaled:scale`
mean_om = attributes(predictors$log_om_s)$`scaled:center`
sd_om = attributes(predictors$log_om_s)$`scaled:scale`
mean_gpp = attributes(predictors$log_gpp_s)$`scaled:center`
sd_gpp = attributes(predictors$log_gpp_s)$`scaled:scale`

# labels
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")

log_gpp_s = tibble(log_gpp_s = quantile(predictors$log_gpp_s,
                                        probs = c(0.25, 0.5, 0.75)))

log_om_s = tibble(log_om_s = quantile(predictors$log_om_s,
                                      probs = c(0.25, 0.5, 0.75)))

# grid to predict over
data_list = list()

for(i in 1:length(fit_list)){
  data_list[[i]] = tibble(mat_s = seq(min(fit_list[[i]]$data$mat_s), max(fit_list[[i]]$data$mat_s),
                                 length.out = 30)) %>% 
    expand_grid(log_gpp_s) %>% 
    expand_grid(log_om_s) %>% 
    mutate(no_m2 = 1,
           xmin = min(fit_list[[i]]$data$xmin),
           xmax = max(fit_list[[i]]$data$xmax),
           lat = median(fit_list[[i]]$data$lat),
           long = median(fit_list[[i]]$data$long)) 
}



# get posteriors ----------------------------------------------------------
post_list = list()

for(i in 1:length(fit_list)){
  post_list[[i]] = data_list[[i]] %>% 
    add_epred_draws(fit_list[[i]], re_formula = NA) %>% 
    group_by(mat_s, log_gpp_s, log_om_s) %>% 
    mutate(model = unique(fit_list[[i]]$model))
}

post_fits = bind_rows(post_list)

# saveRDS(post_fits, file = "posteriors/post_fits.rds")

# get posts for full model to plot behind the alternative fits
mods = sapply(fit_list, function(x) unique(x$model))

full_model = data_list[[1]] %>% 
  add_epred_draws(readRDS("models/fit_temp_om_gpp_year.rds"), re_formula = NA, ndraws = 1000) %>%
  expand_grid(model = mods)

# wrangle posteriors --------------------------------------------------------------------

# post_fits = readRDS(file = "posteriors/post_fits.rds")

post_fits_summarized = post_fits %>% 
  group_by(mat_s, log_gpp_s, log_om_s, model) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp,
         gpp = (log_gpp_s*sd_gpp) + mean_gpp,
         om = (log_om_s*sd_om) + mean_om) %>% 
  ungroup %>% 
  mutate(quantile_om = case_when(om == min(om) ~ "Low OM",
                                 om == max(om) ~ "High OM",
                                 TRUE ~ "Median OM"),
         quantile_gpp = case_when(gpp == min(gpp) ~ "Low GPP",
                                  gpp == max(gpp) ~ "High GPP",
                                  TRUE ~ "Median GPP")) %>% 
  left_join(facet_gpp) %>% 
  left_join(facet_om) %>% 
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
  mutate(model = as.factor(model),
         model = fct_relevel(model, "First dates only", "Last dates only", "Macroinvertebrates", "Fish", "Surber"))


full_posts_summarized = full_model %>% 
  group_by(mat_s, log_gpp_s, log_om_s, model) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp,
         gpp = (log_gpp_s*sd_gpp) + mean_gpp,
         om = (log_om_s*sd_om) + mean_om) %>% 
  ungroup %>% 
  mutate(quantile_om = case_when(om == min(om) ~ "Low OM",
                                 om == max(om) ~ "High OM",
                                 TRUE ~ "Median OM"),
         quantile_gpp = case_when(gpp == min(gpp) ~ "Low GPP",
                                  gpp == max(gpp) ~ "High GPP",
                                  TRUE ~ "Median GPP")) %>% 
  left_join(facet_gpp) %>% 
  left_join(facet_om) %>% 
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
  mutate(model = as.factor(model),
         model = fct_relevel(model, "First dates only", "Last dates only", "Macroinvertebrates", "Fish", "Surber"))

# plot --------------------------------------------------------------------

macro_first_last_plot = post_fits_summarized %>% 
  ggplot(aes(x = mat, y = .epred)) +
  geom_ribbon(data = full_posts_summarized, aes(ymin = .lower, ymax = .upper),
              fill = "black", alpha = 0.8) + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper, fill = model), alpha = 0.4) + 
  geom_line(color = "white") +
  theme_default() + 
  guides(color = "none",
         fill = "none") +  
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (\u00b0C)") +
  theme(strip.text = element_text(hjust = 0)) + 
  coord_cartesian(ylim = c(-3.5, -0.5)) +
  facet_grid2(model ~ facet_omgpp, scales = "free_y") +
  theme(strip.text = element_text(size = 9))


ggsave(macro_first_last_plot, file = "plots/macro_first_last_plot.jpg",
       width = 6, height = 8.5, dpi = 400)



# summarize ---------------------------------------------------------------

# slopes
slope_grid = list()

for(i in 1:length(fit_list)){
  slope_grid[[i]] = tibble(mat_s = c(0, 1)) %>% 
    expand_grid(log_gpp_s) %>% 
    expand_grid(log_om_s) %>% 
    mutate(no_m2 = 1,
           xmin = min(fit_list[[i]]$data$xmin),
           xmax = max(fit_list[[i]]$data$xmax),
           lat = median(fit_list[[i]]$data$lat),
           long = median(fit_list[[i]]$data$long)) 
}


slope_posts = list()

for(i in 1:length(fit_list)){
  slope_posts[[i]] = add_epred_draws(slope_grid[[i]], fit_list[[i]], re_formula = NA) %>% 
    mutate(model = fit_list[[i]]$model) %>% 
    mutate(mat = (mat_s*sd_temp) + mean_temp,
           gpp = (log_gpp_s*sd_gpp) + mean_gpp,
           om = (log_om_s*sd_om) + mean_om) %>% 
    ungroup %>% 
    mutate(quantile_om = case_when(om == min(om) ~ "Low OM",
                                   om == max(om) ~ "High OM",
                                   TRUE ~ "Median OM"),
           quantile_gpp = case_when(gpp == min(gpp) ~ "Low GPP",
                                    gpp == max(gpp) ~ "High GPP",
                                    TRUE ~ "Median GPP")) %>% 
    left_join(facet_gpp) %>% 
    left_join(facet_om) %>% 
    mutate(om_gpp = paste(quantile_om, quantile_gpp, sep = "|")) %>% 
    filter(om_gpp %in% c("Low OM|Low GPP", "Median OM|Median GPP", "High OM|High GPP")) %>%
    mutate(facet_omgpp = case_when(grepl("Low", om_gpp) ~ "a) Low GPP|Low OM",
                                   grepl("Median", om_gpp) ~ "b) Median GPP|Median OM",
                                   TRUE ~ "c) High GPP|High OM")) 
}

bind_rows(slope_posts) %>% 
  ungroup %>% 
  select(mat_s, facet_omgpp, .draw, .epred, model) %>% 
  pivot_wider(names_from = mat_s, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  group_by(facet_omgpp, model) %>% 
  reframe(prob_positive = mean(slope>0)) %>% 
  arrange(model)




