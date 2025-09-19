library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggh4x)
library(viridis)
library(rnaturalearth)

fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")

dat_2022_clauset = readRDS("data/dat_2022_clauset.rds") 

# wrangle posteriors --------------------------------------------------------

data_grid_time = fit_temp_om_gpp$data %>% glimpse %>% 
  select(-dw, -no_m2) %>% distinct() %>% 
  mutate(no_m2 = 1) %>% 
  left_join(dat_2022_clauset %>% ungroup %>% distinct(collect_date, sample_id)) %>% 
  group_by(sample_id) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  mutate(site_id = as.factor(site_id),
         site_id = fct_reorder(site_id, mat)) %>% 
  group_by(year, site_id) 

site_mean_lambda = fit_temp_om_gpp$data %>% 
  select(-dw, -no_m2, xmin, xmax) %>% 
  distinct(site_id, xmin, xmax, log_om_s, log_gpp_s, mat_s) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~(1|site_id)) %>%
  expand_grid(date = seq(min(dat_2022_clauset$collect_date),
                         max(dat_2022_clauset$collect_date),
                         length.out = 2)) %>% 
  group_by(site_id, date) %>% 
  reframe(low50 = quantile(.epred, probs = 0.25),
          high50 = quantile(.epred, probs = 0.75),
          low95 = quantile(.epred, probs = 0.025),
          high95 = quantile(.epred, probs = 0.975),
          .epred = median(.epred)) %>% 
  mutate(year = year(date)) %>% 
  mutate(median = .epred) 

sample_posts = data_grid_time %>% 
  left_join(site_mean_lambda %>% ungroup %>% distinct(site_id, median)) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL)


# get variation among years within sites ----------------------------------
sample_posts %>% 
  group_by(site_id, year, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(site_id, .draw) %>% 
  reframe(sd = sd(.epred)) %>% 
  median_qi(sd, na.rm = T)

# make time series plot ---------------------------------------------------

time_series_lines = sample_posts %>% 
  # filter(site_id != "WLOU") %>%
  group_by(site_id, year, .draw, median) %>%
  reframe(.epred = mean(.epred, na.rm = T)) %>%
  ggplot(aes(x = year, y = .epred, fill = median)) + 
  geom_lineribbon(data = site_mean_lambda, aes(ymin = low95, ymax = high95),
                 alpha = 0.8,
                 color = "white", 
                 linewidth = 0.2) + 
  stat_pointinterval(size = 0.2) + 
  # geom_line(data = . %>% filter(.draw <= 300), 
  #           aes(group = .draw),
  #           alpha = 0.1) +
  facet_wrap(~reorder(site_id, -median), ncol = 5) +
  theme_default() +
  labs(y = "\u03bb",
       x = "Collection Date", 
       fill = "\u03bb") +
  scale_fill_viridis() +
  scale_x_continuous(breaks = c(2016, 2018, 2020, 2022)) +
  coord_cartesian(ylim = c(-3, -1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank()) +
  guides(fill = "none")

ggsave(time_series_lines, file = "plots/time_series_lines.jpg", width = 7, height = 7)

# make map ----------------------------------------------------------------

# load data
neon_latlong <- read_csv(file = "data/raw_data/site_lat_longs.csv") %>% distinct(siteID, lat, long) %>% 
  clean_names()

lambda_site_medians = site_mean_lambda %>% 
  filter(date == min(date)) %>% 
  left_join(neon_latlong)

# load map data
world <- map_data("world") 
states <- map_data("state") 
usa <- ne_countries(scale='medium',returnclass = 'sf') 

#1) Make Map
(map_lambda <- usa %>% 
    filter(sovereignt == "United States of America") %>% 
    ggplot() + 
    # coord_sf() + 
    geom_polygon(data = world, aes(x = long, y = lat, group = group), fill = "grey70") +
    geom_sf(color = "white", fill = "grey70") +
    geom_polygon(data = states, aes(x = long, y = lat, group = group), color = "white", fill = "grey70")  +
    geom_point(data = lambda_site_medians, 
               aes(x = long, y = lat, fill = .epred),
               size = 2,
               alpha = 0.9,
               color = "black", shape = 21,
               position = position_jitter(width = 2, height = 1, seed = 2323)) +
    theme_void() +
    coord_sf(ylim = c(10, 68), xlim = c(-160, -68)) +
    scale_fill_viridis(option = "D", begin = 0, end = 1) + 
    theme(legend.position = c(0.25, 0.4),
          legend.key.size = unit(0.4, "cm"),
          legend.title = element_text(size = 8),
          legend.text = element_text(size = 8),
          text = element_text(family = "serif")) +
    guides(fill = "none") +
    NULL)

ggsave(map_lambda, file = "plots/map_lambda.jpg", width = 6.5, height = 8)


# combine plot and map ----------------------------------------------------
time_series_lines_a = sample_posts %>% 
  # filter(site_id != "WLOU") %>%
  group_by(site_id, year, .draw, median) %>%
  reframe(.epred = mean(.epred, na.rm = T)) %>%
  ggplot(aes(x = year, y = .epred, fill = median)) + 
  geom_lineribbon(data = site_mean_lambda, aes(ymin = low95, ymax = high95),
                  alpha = 0.8,
                  color = "white", 
                  linewidth = 0.2) + 
  stat_pointinterval(size = 0.2) + 
  # geom_line(data = . %>% filter(.draw <= 300), 
  #           aes(group = .draw),
  #           alpha = 0.1) +
  facet_wrap(~reorder(site_id, -median), ncol = 6) +
  theme_default() +
  labs(y = "\u03bb",
       x = "Collection Date", 
       fill = "\u03bb") +
  scale_fill_viridis() +
  scale_x_continuous(breaks = c(2016, 2018, 2020, 2022)) +
  coord_cartesian(ylim = c(-3, -1)) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5),
        strip.text = element_text(size = 8),
        legend.text = element_text(size = 8),
        axis.title.x = element_blank()) +
  guides(fill = "none")


layout <- c(
  area(t = 3, l = 1, b = 5, r = 3),
  area(t = 1, l = 3, b = 5, r = 6)
)
# Show the layout to make sure it looks as it should
plot(layout)

time_series_lines_a + map_lambda + plot_layout(design = layout)

