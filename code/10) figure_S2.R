library(tidyverse)
library(brms)
library(ggthemes)
library(patchwork)

### Plots with fish, inverts, and fish + inverts combined

# isd_by_temp -------------------------------------------------------------

isd_by_temp_fish = readRDS(file = "plots/isd_by_temp-fishonly.rds")
isd_by_temp_inverts = readRDS(file = "plots/isd_by_temp-invertsonly.rds")
isd_by_temp = readRDS(file = "plots/isd_by_temp.rds")

points_fish = layer_data(isd_by_temp_fish, 1) %>% mutate(animal_type = "fish")
points_inverts = layer_data(isd_by_temp_inverts, 1) %>% mutate(animal_type = "inverts")
points_all = layer_data(isd_by_temp, 1) %>% mutate(animal_type = "inverts + fish", shape = 21)

lines_fish = layer_data(isd_by_temp_fish, 2) %>% mutate(animal_type = "fish")
lines_inverts = layer_data(isd_by_temp_inverts, 2) %>% mutate(animal_type = "inverts")
lines_all = layer_data(isd_by_temp, 2) %>% mutate(animal_type = "inverts + fish")

ribbon_fish = layer_data(isd_by_temp_fish, 3) %>% mutate(animal_type = "fish")
ribbon_inverts = layer_data(isd_by_temp_inverts, 3) %>% mutate(animal_type = "inverts")
ribbon_all = layer_data(isd_by_temp, 3) %>% mutate(animal_type = "inverts + fish")

points = bind_rows(points_fish,
                   points_inverts,
                   points_all) %>% 
  mutate(panel = case_when(animal_type == "fish" ~ "c) Fish",
                           animal_type == "inverts" ~ "b) Invertebrates",
                           TRUE ~ "a) Inverts + Fish"))

lines = bind_rows(lines_fish,
                  lines_inverts,
                  lines_all) %>% 
  mutate(panel = case_when(animal_type == "fish" ~ "c) Fish",
                           animal_type == "inverts" ~ "b) Invertebrates",
                           TRUE ~ "a) Inverts + Fish"))

ribbons = bind_rows(ribbon_fish,
                    ribbon_inverts,
                    ribbon_all) %>% 
  mutate(panel = case_when(animal_type == "fish" ~ "c) Fish",
                           animal_type == "inverts" ~ "b) Invertebrates",
                           TRUE ~ "a) Inverts + Fish"))


# plot

# three column

(all_isd = ggplot(data = points, aes(x = x, y = y)) + 
    geom_pointrange(aes(ymin = ymin, ymax = ymax),
                    alpha = 0.2, size = 0.01,
                    linewidth = 0.1) + 
    geom_line(data = lines) + 
    geom_ribbon(data = ribbons, aes(ymin = ymin, ymax = ymax), alpha = 0.2) + 
    scale_color_colorblind() + 
    scale_fill_colorblind() + 
    facet_wrap(~panel, ncol = 3) +
    theme_default() + 
    theme(strip.text = element_text(hjust = 0)) +
    guides(color = "none",
           fill = "none") +
    coord_cartesian(ylim = c(-2, -0.4)) +
    labs(y = "\u03bb (ISD exponent)",
         x = "Mean Annual Temperature (\u00b0C)")
)

saveRDS(all_isd, file = "plots/all_isd.rds")
ggview(all_isd, width = 6.5, height = 2)
ggsave(all_isd, file = "plots/all_isd.jpg",dpi = 500,
       width = 6.5, height = 2)