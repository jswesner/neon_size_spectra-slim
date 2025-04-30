library(tidyverse)
library(brms)
library(tidybayes)
library(isdbayes)
# Literature Figure Comparisons -------------------------------------------
theme_set(theme_default())

fit_pareto = readRDS("models/fit_temp_om_gpp.rds")

predictors = readRDS("data/predictors_scaled.rds")

mean_temp = attributes(predictors$mat_s)$`scaled:center`
sd_temp = attributes(predictors$mat_s)$`scaled:scale`
mean_om = attributes(predictors$log_om_s)$`scaled:center`
sd_om = attributes(predictors$log_om_s)$`scaled:scale`
mean_gpp = attributes(predictors$log_gpp_s)$`scaled:center`
sd_gpp = attributes(predictors$log_gpp_s)$`scaled:scale`

# literature comparison

lit <- read_csv("data/temp_summaries_table.csv") %>% 
        filter(Include == "Y") %>% 
  mutate(size_magnitude = round(log10(xmax) - log10(xmin),0),
         size_scale = size_magnitude*scale_km) %>% 
  group_by(Author, organisms, b_diff) %>% 
  mutate(id = cur_group_id())

mod_best <- readRDS("models/fit_temp_om_gpp.rds")
mod_summary = summary(mod_best)
conds = tibble(mat_s = seq(min(fit_pareto$data$mat_s), max(fit_pareto$data$mat_s), length.out = 20)) %>% 
  mutate(log_gpp_s = 0,
         log_om_s = 0,
         no_m2 = 1,
         xmin = 0.0035, 
         xmax = 200000) %>% 
  add_epred_draws(fit_pareto, re_formula = NA) %>% 
  group_by(mat_s) %>% 
  median_qi(.epred)

# plot slopes (scale-independent)
conds_scaled = conds %>% 
  mutate(value = .epred - mean(.epred),
         .lower = .lower - mean(.epred),
         .upper = .upper - mean(.epred),
         Driver = "Temperature",
        x = mat_s,
        id = "Gjoni et al. 2023",
        group = "This Study") %>% 
  mutate(x_raw = (mat_s*sd_temp) + mean_temp)

(lit_plot_unscaled = lit %>% 
  filter(Driver == "Temperature") %>% 
  filter(Author != "Gjoni et al. 2023") %>% 
  mutate(low = 0 - 0.5*direction,
         high = 0 + 0.5*direction,
         group = "Literature Estimates") %>% 
  pivot_longer(cols = c(low, high)) %>% 
  mutate(x = case_when(name == "low" ~ temp_low, TRUE ~ temp_high)) %>% 
  ggplot(aes(x = x, y = value)) + 
  geom_line(aes(group = id), alpha = 0.5) +
  # facet_wrap(~Driver) + 
  geom_line(data = conds_scaled, aes(x = x_raw)) +
  geom_ribbon(data = conds_scaled,
              aes(ymin = .lower,
                  ymax = .upper,
                  x = x_raw),
              alpha = 0.7,
              fill = "#E69F00") +
  geom_text(data = . %>% filter(x == min(x)) %>% distinct(x, Author, value) %>% 
                             group_by(Author) %>%
                             mutate(author_code = cur_group_id()), 
                           aes(label = author_code),
                           size = 3) +
  labs(y = "\u03bb (scaled)",
       x = "Temperature (\u00b0C)"))


ggsave(lit_plot_unscaled, file = "plots/lit_plot_unscaled.jpg", 
       width = 5, height = 5)

