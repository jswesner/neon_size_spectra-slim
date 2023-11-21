library(brms)
library(tidyverse)
library(isdbayes)

#1) load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds")

#2) refit to sample from priors only
# fit_priors = update(fit_temp_om_gpp, sample_prior = "only", iter = 100, chains = 1)

fit_priors = readRDS(file = "models/fit_priors.rds")

#3) Get conditional effects from priors and posts

prior_conds = tibble(mat_s = seq(from = min(fit_priors$data$mat_s), 
                                 to = max(fit_priors$data$mat_s),
                                 length.out = 20)) %>% 
  mutate(xmin = 1, xmax = 1000, no_m2 = 1) %>%  # placeholders. Values don't matter here
  mutate(log_om_s = 0, log_gpp_s = 0) %>% 
  add_epred_draws(fit_priors, re_formula = NA) %>% 
  mutate(model = "a) Prior")
  

post_conds = tibble(mat_s = seq(from = min(fit_priors$data$mat_s), 
                                 to = max(fit_priors$data$mat_s),
                                 length.out = 20)) %>% 
  mutate(xmin = 1, xmax = 1000, no_m2 = 1) %>%  # placeholders. Values don't matter here
  mutate(log_om_s = 0, log_gpp_s = 0) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NA) %>% 
  mutate(model = "b) Posterior")
  
#4) combine prior and posterior conditionals
prior_post = bind_rows(prior_conds, post_conds)

#5) Plot

prior_post_plot = prior_post %>%
  filter(.draw <= 200) %>% 
  ggplot(aes(x = mat_s, y = .epred)) + 
  geom_line(aes(group = .draw), linewidth = 0.1) +
  facet_wrap(~model) +
  labs(y = "\u03bb",
       x = "Mean Annual Temperature (z-score)") +
  theme_default()


ggview::ggview(prior_post_plot, width = 6.5, height = 3)
ggsave(prior_post_plot, width = 6.5, height = 3, dpi = 500,
       file = "plots/ms_plots/figure_s3.jpg")
saveRDS(prior_post_plot, file = "plots/ms_plots/figure_s3.rds")

