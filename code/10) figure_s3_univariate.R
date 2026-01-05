library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)
library(ggtext)

# Figure S3 ---------------------------------------------------------
#1) Get posterior summaries with conditional_effects()
# univariate plot ---------------------------------------------------------
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")
fit_temp = readRDS("models/fit_temp_year.rds")
predictors = readRDS("data/predictors_scaled.rds")

mean_temp = attributes(predictors$mat_s)$`scaled:center`
sd_temp = attributes(predictors$mat_s)$`scaled:scale`
mean_om = attributes(predictors$log_om_s)$`scaled:center`
sd_om = attributes(predictors$log_om_s)$`scaled:scale`
mean_gpp = attributes(predictors$log_gpp_s)$`scaled:center`
sd_gpp = attributes(predictors$log_gpp_s)$`scaled:scale`

#2) Get individual lambda posterior summaries with add_epred_draws()
sample_dots = fit_temp_om_gpp$data %>% 
  distinct(sample_id, mat_s, log_gpp_s, log_om_s, year, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp)

sample_dots_summary = sample_dots %>% 
  group_by(sample_id, mat_s, year) %>% 
  median_qi(.epred) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) 

site_dots = fit_temp_om_gpp$data %>% 
  distinct(mat_s, log_gpp_s, log_om_s, site_id, xmin, xmax, year) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~ (1|site_id)) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  group_by(mat, site_id, year) %>% 
  median_qi(.epred)

univariate_posts = tibble(mat_s = seq(min(fit_temp$data$mat_s),
                                      max(fit_temp$data$mat_s),
                                      length.out = 30)) %>% 
  mutate(xmin = min(fit_temp$data$xmin),
         xmax = max(fit_temp$data$xmax),
         no_m2 = 1) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  add_epred_draws(fit_temp, re_formula = NA)

#3) Make the marginal plot of temperature and lambda
uni_plot_dots = univariate_posts %>%  
  ggplot(aes(x = mat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6, fill = "black", color = "white",
                  linewidth = 0.2) +
  stat_pointinterval(data = site_dots, aes(y = .epred, group = site_id),
                     size = 0.05, shape = 20,
                     geom = "pointrange") +
  geom_point(data = sample_dots_summary, aes(y = .epred),
             size = 0.05, shape = 1, position = position_jitter(width = 0.2,
                                                                height = 0)) +
  geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
                                        group = site_id),
                  size = 0.1, linewidth = 0.2, shape = 20) +
  theme_default() + 
  labs(y = "\u03bb (ISD exponent)",
       x = "Mean Annual Temperature (\u00b0C)") +
  # guides(color = "none")+ 
  viridis::scale_color_viridis(direction = 1,
                               breaks = seq(27, 1, by = -5)) +
  theme(legend.position = c(0.85, 0.2),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 8)) +
  guides(color = "none")

saveRDS(uni_plot_dots, file = "plots/uni_plot_dots.rds")
ggsave(uni_plot_dots, file = "plots/uni_plot_dots.jpg", width = 5, height = 5, dpi = 400)


as_draws_df(fit_temp) %>% 
  as_tibble() %>% 
  reframe(prob_pos = sum(b_mat_s > 0)/max(nrow(.)))


# by year -----------------------------------------------------------------

univariate_posts_year = tibble(mat_s = seq(min(fit_temp$data$mat_s),
                                      max(fit_temp$data$mat_s),
                                      length.out = 30)) %>% 
  mutate(xmin = min(fit_temp$data$xmin),
         xmax = max(fit_temp$data$xmax),
         no_m2 = 1) %>% 
  mutate(mat = (mat_s*sd_temp) + mean_temp) %>% 
  expand_grid(year = unique(fit_temp$data$year)) %>% 
  add_epred_draws(fit_temp, re_formula = ~ (1 + mat_s|year))


univariate_plot_years = univariate_posts_year %>%  
  ggplot(aes(x = mat, y = .epred)) +
  stat_lineribbon(.width = 0.95, alpha = 0.6, fill = "black", color = "white",
                  linewidth = 0.2) +
  stat_pointinterval(data = site_dots, aes(y = .epred, group = site_id),
                     size = 0.05, shape = 20,
                     geom = "pointrange") +
  geom_point(data = sample_dots_summary, aes(y = .epred),
             size = 0.05, shape = 1, position = position_jitter(width = 0.2,
                                                                height = 0)) +
  geom_pointrange(data = site_dots, aes(y = .epred, ymin = .lower, ymax = .upper,
                                        group = site_id),
                  size = 0.1, linewidth = 0.2, shape = 20) +
  theme_default() + 
  labs(y = "\u03bb (ISD exponent)",
       x = "Mean Annual Temperature (\u00b0C)") +
  # guides(color = "none")+ 
  viridis::scale_color_viridis(direction = 1,
                               breaks = seq(27, 1, by = -5)) +
  theme(legend.position = c(0.85, 0.2),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.3, "cm"),
        legend.title = element_text(size = 8)) +
  guides(color = "none") +
  facet_wrap(~year)
