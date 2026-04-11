library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
library(viridis)
library(patchwork)
theme_set(brms::theme_default())

#1) Load models and data
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")
dat_2022_clauset = as_tibble(fit_temp_om_gpp$data)

# refit with default priors
# fit_temp_om_gpp_wide = brm(dw | vreal(no_m2, xmin, xmax) ~ log_om_s*mat_s*log_gpp_s + (1 |site_id:sample_id) + (1 + log_om_s*mat_s*log_gpp_s|year),
#                            data = dat_2022_clauset,
#                            stanvars = stanvars,     # required for truncated Pareto via isdbayes package
#                            family = paretocounts(), # required for truncated Pareto via isdbayes package
#                            # prior = c(prior(normal(-2, 0.5), class = "Intercept"),  # priors are silenced so that the model uses default priors
#                            #           prior(normal(0, 0.2), class = "b"),
#                            #           prior(exponential(7), class = "sd")),
#                            iter = 1000,
#                            cores = 4,
#                            threads = 12,
#                            chains = 1)
# 
# saveRDS(fit_temp_om_gpp_wide, file = "models/fit_temp_om_gpp_wide.rds")
fit_temp_om_gpp_wide = readRDS(file = "models/fit_temp_om_gpp_wide.rds")

posts_informative = as_draws_df(fit_temp_om_gpp) %>% ungroup() %>% select(.draw, starts_with("b_")) %>% mutate(priors = "informative")
posts_default = as_draws_df(fit_temp_om_gpp) %>% ungroup() %>% select(.draw, starts_with("b_")) %>% mutate(priors = "less informative")

prior_sensitivity_plot = bind_rows(posts_default, posts_informative) %>% 
  pivot_longer(cols = c(-.draw, -priors)) %>% 
  ggplot(aes(y = value, x = name, fill = priors)) +
  geom_violin(width = 0.6) +
  coord_flip() +
  scale_fill_colorblind() +
  labs(fill = "Priors",
       x = "",
       y = "Parameter value")

ggsave(prior_sensitivity_plot, file = "plots/prior_sensitivity_plot.jpg")

post_pred_gm = readRDS(file = "plots/fig_s7_post_pred_gm.rds")

a = (post_pred_gm + labs(subtitle = "a)"))
b = (prior_sensitivity_plot + labs(subtitle = "b)",
                                   x = "Parameter"))

post_pred_priors = a/b
ggsave(post_pred_priors, file = "plots/fig_s7_post_pred_gm_priors.jpg", width = 6.5, height = 6)
