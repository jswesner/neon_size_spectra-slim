library(tidyverse)
library(janitor)
library(viridis)
library(rnaturalearth)
library(patchwork)
theme_set(brms::theme_default())

# load data
neon_latlong <- read_csv(file = "data/raw_data/site_lat_longs.csv") %>% distinct(siteID, lat, long) %>% 
  clean_names() %>% 
  mutate(long = case_when(site_id == "GUIL" ~ long + 0.5,
                          TRUE ~ long))

predictors = readRDS("data/predictors_scaled.rds") 

temp_gpp_om = predictors %>% 
  ungroup %>% distinct(site_id, temp_deg_c,
                       gpp, om) %>% 
  left_join(neon_latlong) %>% 
  pivot_longer(cols = c(gpp, om, temp_deg_c))

# load map data
world <- map_data("world") %>% expand_grid(name = temp_gpp_om %>% distinct(name))
states <- map_data("state") %>% expand_grid(name = temp_gpp_om %>% distinct(name))
usa <- ne_countries(scale='medium',returnclass = 'sf') 

#1) Make Map
(map_temp <- usa %>% 
    filter(sovereignt == "United States of America") %>% 
    ggplot() + 
    # coord_sf() + 
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey70") +
    geom_sf(color = "white", fill = "grey70") +
    geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "white", fill = "grey70")  +
    geom_point(data = temp_gpp_om %>% filter(name == "temp_deg_c"), 
               aes(x = long, y = lat, fill = value),
               size = 2,
               alpha = 0.9,
               color = "black", shape = 21,
               # position = position_jitter(width = 2, height = 1, seed = 2323)
               ) +
    theme_void() +
    coord_sf(ylim = c(10, 68), xlim = c(-160, -68)) +
    labs(fill = "\u00b0C",
         title = "a) Site Map") +
    scale_fill_viridis() + 
    theme(legend.position = c(0.25, 0.4),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          text = element_text(family = "serif")) +
    NULL)

#2) Get correlations
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

abiotic_wide = fit_temp_om_gpp$data %>% 
  distinct(mat_s, log_gpp_s, log_om_s) 

correlations = tibble(gppom = cor(abiotic_wide$log_gpp_s, abiotic_wide$log_om_s),
                      gpptemp = cor(abiotic_wide$log_gpp_s, abiotic_wide$mat_s),
                      omtemp = cor(abiotic_wide$mat_s, abiotic_wide$log_om_s)) %>% 
  mutate(across(where(is.numeric), round, 2))  %>% 
  pivot_longer(cols = everything()) %>% 
  mutate(labels = paste0("R = ", value))

#3) Make correlation plots

omtemp_plot = abiotic_wide %>% 
  ggplot(aes(x = mat_s, y = log_om_s)) + 
  geom_point(size = 1) +
  labs(y = "log OM\n(z-score)",
       x = "Temperature (z-score)") +
  labs(title = "b) Correlation Plots") +
  annotate(geom = "text", x= -0.9, y = 2, label = correlations %>% filter(name == "omtemp") %>% pull(labels),
           size = 3)

gpptemp_plot = abiotic_wide %>% 
  ggplot(aes(x = mat_s, y = log_gpp_s)) + 
  geom_point(size = 1) +
  labs(y = "log GPP\n(z-score)",
       x = "Temperature (z-score)") +
  annotate(geom = "text", x= -0.9, y = 2, label = correlations %>% filter(name == "gpptemp") %>% pull(labels),
           size = 3)

gppom_plot = abiotic_wide %>% 
  ggplot(aes(x = log_gpp_s, y = log_om_s)) + 
  geom_point(size = 1) +
  labs(y = "log OM\n(z-score)",
       x = "log GPP (z-score)") +
  annotate(geom = "text", x= -0.9, y = 2, label = correlations %>% filter(name == "gppom") %>% pull(labels),
           size = 3)

correlation_plots = omtemp_plot/gpptemp_plot/gppom_plot 


#4) Combine map and correlation plots

map_correlation_plot = map_temp + correlation_plots + 
  plot_layout(widths = c(3, 1), heights = c(1, 1)) 

ggsave(map_correlation_plot, file = "plots/ms_plots/map_correlation_plot.jpg", width = 6.5, height = 8)
saveRDS(map_correlation_plot, file = "plots/ms_plots/map_correlation_plot.rds")




