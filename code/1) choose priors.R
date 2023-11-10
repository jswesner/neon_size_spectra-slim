library(brms)
library(tidyverse)
library(tidybayes)
library(ggthemes)
library(janitor)

#1) load data (NEON body size data)
dat_all = readRDS("data/derived_data/dat_all.rds")

#2) fit prior predictive using sample_prior = "only" (note: a quicker way without varying intercepts is below)
prior_sims = brm(dw | vreal(no_m2, xmin, xmax) ~ log_om_s*mat_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id),
    data = dat_all,
    stanvars = stanvars,    # required for truncated Pareto via isdbayes package
    family = paretocounts(),# required for truncated Pareto via isdbayes package
    prior = c(prior(normal(-1.5, 0.2), class = "Intercept"),
              prior(normal(0, 0.1), class = "b"),
              prior(exponential(7), class = "sd")),
    iter = 1000, 
    sample_prior = "only",
    file = "models/prior_sims.rds",          # The model is already saved. This model will load, but not fit. It will only refit if the next line (file_refit) is removed.
    file_refit = "on_change",
    chains = 1)

#3) plot prior regressions

cond_priors = plot(conditional_effects(prior_sims, spaghetti = T, ndraws = 100, ask = F))  # get conditional effects

# plot temp conditional effect
cond_priors$mat_s +
  labs(y = "\u03bb (ISD exponent)") 

# plot om conditional effect
cond_priors$log_om_s +
  labs(y = "\u03bb (ISD exponent)") 

# plot any other two-way conditional effects
cond_priors$`log_om_s:log_gpp_s` +
  labs(y = "\u03bb (ISD exponent)") 


#4) plot prior lambdas
# check sample specific lambdas

prior_lambdas = prior_sims$data %>% 
  select(-dw, -no_m2, -xmax, -xmin) %>%    # select out to make distinct sample grid (without individual size information)
  distinct() %>% 
  mutate(xmin = 0.003, xmax = 200000 , no_m2 = 1) %>% 
  add_epred_draws(prior_sims, re_formula = NULL)


prior_lambdas %>% 
  ggplot(aes(x = mat_s, y = .epred)) + 
  geom_point(shape = 21, alpha = 0.2)


#5) Simulate priors without using brms (quicker but doesn't include varying intercepts). This allows for easy changes to the prior without having to use brm()
# simulate priors
nsims = 100
int_sd = 0.2   # standard deviation of the intercept
b_sd = 0.1     # standard deviation of the slopes (all betas)
priors = tibble(iter = 1:nsims,
                a = rnorm(nsims, -1.5, int_sd),
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
  left_join(facet_om)



#6) Plot priors
prior_sims_byhand %>% 
  ggplot(aes(x = mat_s, y = lambda, color = -log_om_s)) + 
  geom_line(aes(group = iter), alpha = 0.3) + 
  facet_grid(facet_gpp ~ facet_om, labeller = "label_parsed") +
  labs(y = "\u03bb (ISD exponent)",
       title = "Prior Predictive") +
  scale_color_viridis() +
  geom_abline(intercept = -1.3, slope = -0.1, color = "red")  # approximate slope from Pomeranz et al. (2022). Check as a reference.





