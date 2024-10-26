library(tidyverse)
library(ggridges)
library(latex2exp)
library(tidybayes)
library(brms)
library(MCMCglmm)
library(patchwork)
library(isdbayes)


# Make Simulation Plot ----------------------------------------------------

#1) Simulate values of ppmr, te, metab scaling, and subsidies
set.seed(3333333)

ppmr = rlnorm(1000, 12, 2.4) # range between 10^2 and 10^7
# ppmr = rtnorm(1000, 10^5, 50000, lower = 0) # range between 10^2 and 10^7
te = rbeta(1000, 7, 60)
metab = rbeta(1000, 12, 15)
subsidies = rbeta(1000, 14, 25)

#2) Combine simulated values to a tibble
isd_sims = tibble(ppmr = ppmr,
                  te = te,
                  metab = metab,
                  subsidies = subsidies) %>% 
  mutate(`1) Standard Model` = (log10(te)/log10(ppmr)) - 0.75 - 1,
         `3) Standard Model +\nSubsidies` = (log10(te)/log10(ppmr)) - 0.75 + subsidies - 1,
         `2) Standard Model +\nShallow Metabolic Scaling` = (log10(te)/log10(ppmr)) - metab - 1,
         `4) Standard Model +\nShallow Metabolic Scaling +\nSubsidies` = (log10(te)/log10(ppmr)) - metab + subsidies - 1)

isd_sims %>% 
  pivot_longer(cols = c(-ppmr, -te, -metab, -subsidies)) %>% 
  group_by(name) %>% 
  median_qi(value)

#3) Load fitted model and wrangle empirical lambdas
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")

post_sample_lambdas_summary = fit_temp_om_gpp$data %>% 
  distinct(sample_id, site_id, year, log_om_s, log_gpp_s, mat_s) %>% 
  mutate(xmin = 0.003, xmax = 200000, no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL) %>% 
  group_by(sample_id) %>% 
  median_qi(.epred)

#4) Make the simulation plot (panel a) 

simulation_plot = isd_sims %>% 
  pivot_longer(cols = c(-ppmr, -te, -metab, -subsidies)) %>% 
  ggplot(aes(x = value, y = fct_rev(as.factor(name)))) + 
  geom_vline(data = post_sample_lambdas_summary, 
             aes(xintercept = .epred), alpha = 0.1, color = "#1f968bff",
             linewidth = 0.2) +
  stat_binline(bins = 100, scale = 1, draw_baseline = FALSE) + 
  theme_default() +
  labs(x = "\u03bb (ISD exponent)",
       # subtitle = TeX("$\\lambda = \\frac{log_{10}\\alpha}{log_{10}\\beta} - 0.75 - 1$"),
       y = "") + 
  theme(axis.title.y = element_blank()) + 
  annotate(geom = "text", x = -0.9, y = 4.5, label = "Empirical \u03bb's",
           color = "#1f968bff", size = 2.5) +
  coord_cartesian(xlim = c(-2.4, -0.6))

# ggview::# ggview(simulation_plot, width = 6.5, height = 3.9)
saveRDS(simulation_plot, file = "plots/ms_plots/simulation_plot.rds")
ggsave(simulation_plot, width = 6.5, height = 3.9, file = "plots/ms_plots/simulation_plot.jpg", dpi = 500)



# Make Metabolic Scaling Plot ---------------------------------------------

#5) Load metabolic scaling results and summarize posteriors 
brm_metab <- readRDS("models/brm_metab.rds")

cond_posts_metab = brm_metab$data %>% 
  distinct(heat, fish) %>% 
  expand_grid(log_dw_c = seq(min(brm_metab$data$log_dw_c), max(brm_metab$data$log_dw_c), length.out = 20)) %>% 
  add_epred_draws(brm_metab, re_formula = NA) %>% 
  group_by(log_dw_c, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  group_by(log_dw_c) %>% 
  median_qi(.epred)

#6) Plot metabolic scaling regression
metab_plot = cond_posts_metab %>% 
  ggplot(aes(x = log_dw_c, y = .epred))  +
  geom_point(data = brm_metab$data, aes(y = log_resp_c), 
             shape = 1, size = 0.1) + 
  geom_line() + 
  geom_ribbon(aes(ymin = .lower, ymax = .upper), alpha = 0.3) + 
  labs(y = "log10 Metabolic Rate (centered)",
       x = "log10 dry mass (mg centered)") + 
  geom_abline(slope = 0.75, linewidth = 0.5, linetype = "dotted") + 
  theme_default()

saveRDS(metab_plot, file = "plots/metab_plot.rds")

#7) Summarize metabolic scaling slopes
brm_metab$data %>% 
  distinct(heat, fish) %>% 
  expand_grid(log_dw_c = c(0, 1)) %>% 
  add_epred_draws(brm_metab, re_formula = NA) %>% 
  group_by(log_dw_c, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = log_dw_c, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  median_qi(slope)



# Combine All Plots -------------------------------------------------------
#1) Load plots
simulation_plot = readRDS(file = "plots/ms_plots/simulation_plot.rds") + labs(subtitle = "a) Simulation Study") + 
  theme(axis.text.y = element_text(size = 8))
metab_plot = readRDS(file = "plots/metab_plot.rds") + labs(subtitle = "b) Shallow Metabolic Scaling")

#2) Combine plots
sim_metab_plot = simulation_plot + metab_plot

#3) Save plots
# ggview::# ggview(sim_metab_plot, width = 6.5, height = 3, units = "in")
ggsave(sim_metab_plot, file = "plots/ms_plots/sim_metab_plot.jpg",
       width = 6.5, height = 3, units = "in", dpi = 600)
saveRDS(sim_metab_plot, file = "plots/ms_plots/sim_metab_plot.rds")


# Summarize ---------------------------------------------------------------

# metab slopes
brm_metab$data %>% 
  distinct(heat, fish) %>% 
  expand_grid(log_dw_c = c(0, 1)) %>% 
  add_epred_draws(brm_metab, re_formula = NA) %>% 
  group_by(log_dw_c, .draw) %>% 
  reframe(.epred = mean(.epred)) %>% 
  pivot_wider(names_from = log_dw_c, values_from = .epred) %>% 
  mutate(slope = `1` - `0`) %>% 
  median_qi(slope)


# simulation parameters
parameter = c("alpha", "beta", "sigma", "gamma")
first_moment = c(12, 7, 12, 14)
second_moment = c(3, 60, 15, 25)

simulation_table = tibble(parameter = parameter,
       first_moment = first_moment,
       second_moment = second_moment) %>% 
  group_by(parameter) %>% 
  mutate(mean = case_when(parameter == "alpha" ~ mean(rlnorm(100000, first_moment, second_moment)),
                          TRUE ~ mean(rbeta(100000, first_moment, second_moment))),
         median = case_when(parameter == "alpha" ~ median(rlnorm(100000, first_moment, second_moment)),
                          TRUE ~ median(rbeta(100000, first_moment, second_moment))),
         sd = case_when(parameter == "alpha" ~ sd(rlnorm(100000, first_moment, second_moment)),
                          TRUE ~ sd(rbeta(100000, first_moment, second_moment))))

write_csv(simulation_table, file = "tables/simulation_table.csv")

