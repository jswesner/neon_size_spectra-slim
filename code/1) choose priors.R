library(brms)
library(tidyverse)
library(tidybayes)
library(ggthemes)
library(janitor)

#1) load data (NEON body size data)
dat_all = readRDS("data/derived_data/dat_all.rds")


#2) Simulate priors without using brms (quicker). This allows for easy changes to the prior without having to use brm()
# simulate priors
nsims = 300
int_sd = 0.5   # standard deviation of the intercept
b_sd = 0.1     # standard deviation of the slopes (all betas)
priors = tibble(iter = 1:nsims,
                a = rnorm(nsims, -2, int_sd),
                beta_mat = rnorm(nsims, 0, b_sd),
                beta_om = rnorm(nsims, 0, b_sd),
                beta_gpp = rnorm(nsims, 0, b_sd),
                beta_gpp_om = rnorm(nsims, 0, b_sd),
                beta_gpp_mat = rnorm(nsims, 0, b_sd),
                beta_om_mat = rnorm(nsims, 0, b_sd),
                beta_om_mat_gpp = rnorm(nsims, 0, b_sd))

# labels (For plotting later)
facet_gpp = readRDS(file = "plots/facet_gpp.rds")
facet_om = readRDS(file = "plots/facet_om.rds")


# simulate prior predictive 
prior_sims_byhand = dat_all %>% 
  select(mat_s, site_id) %>% 
  distinct() %>%
  expand_grid(log_om_s = c(-1, 0, 1)) %>%  # add quantiles of gpp and om
  expand_grid(log_gpp_s = c(-1, 0, 1)) %>% 
  expand_grid(priors) %>% 
  mutate(lambda = a + 
           beta_mat*mat_s + 
           beta_gpp*log_gpp_s + 
           beta_om*log_om_s +
           beta_gpp_om*log_gpp_s*log_om_s + 
           beta_gpp_mat*log_gpp_s*mat_s + 
           beta_om_mat*log_om_s*mat_s +
           beta_om_mat_gpp*log_om_s*mat_s*log_gpp_s) %>% 
  mutate(quantile_om = case_when(log_om_s == min(log_om_s) ~ "Low OM",
                                 log_om_s == max(log_om_s) ~ "High OM",
                                 TRUE ~ "Median OM"),
         quantile_gpp = case_when(log_gpp_s == min(log_gpp_s) ~ "Low GPP",
                                  log_gpp_s == max(log_gpp_s) ~ "High GPP",
                                  TRUE ~ "Median GPP")) %>% 
  left_join(facet_gpp) %>% 
  left_join(facet_om) %>% 
  group_by(iter) %>% 
  mutate(r_year = rnorm(1, 0, rexp(1, 7)), # add varying intercepts
         r_site = rnorm(1, 0, rexp(1, 7)),
         r_sample = rnorm(1, 0, rexp(1, 7))) %>% 
  mutate(lambda_varying = lambda + r_year + r_site + r_sample) %>% 
  mutate(quantile_gpp = fct_relevel(quantile_gpp, "Low GPP", "Median GPP"),
         quantile_om = fct_relevel(quantile_om, "Low OM", "Median OM"))

#3) Plot priors
prior_pred = prior_sims_byhand %>% 
  ggplot(aes(x = mat_s, y = lambda_varying, color = -log_om_s)) + 
  geom_line(aes(group = iter), alpha = 0.1) + 
  facet_grid(quantile_gpp ~ quantile_om) +
  labs(y = "\u03bb (ISD exponent)",
       title = "Prior Predictive") +
  scale_color_viridis() +
  guides(color = "none")

ggsave(prior_pred, file = "plots/prior_pred.jpg", width = 6.5, height = 6.5)



