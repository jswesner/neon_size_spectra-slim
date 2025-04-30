library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(janitor)

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp.rds") # this model is correct. It has cari and oksr non-fish samples. I has the correct fish masses (i.e., 0.2wm = dm) and it has no syca
fit_temp = readRDS("models/fit_temp.rds")

# Simulate lambdas from Reuman equation with lambda = -2 and compare to empirical.
#1) Simulate values of ppmr, te, metab scaling, and subsidies
set.seed(3333333)

meanlog = log(10^4)
sdlog = 2
ppmr = rlnorm(1000, meanlog, sdlog) # range between 10^2 and 10^7

exp(meanlog + (sdlog^2/2))
sqrt((exp(sdlog^2) - 1)*exp(2*meanlog + sdlog^2))

a = 7
b = 60
te = rbeta(1000, a, b)
a/(a+b) #mean
sqrt((a*b)/((a+b)^2*(a + b+ 1))) # sd

metab = rbeta(1000, 9.6, 3.25)
a2 = 9.6
b2 = 3.25
a2/(a2+b2) #mean
sqrt((a2*b2)/((a2+b2)^2*(a2 + b2 + 1))) # sd


# lambda predicted by default functioning
(log10(0.1)/log10(10^4)) - 0.75 - 1

#2) Combine simulated values to a tibble

isd_sims = tibble(ppmr = ppmr,
                  te = te,
                  metab = metab) %>% 
  mutate(lambda_sims = (log10(te)/log10(ppmr)) - metab - 1)

site_dots = fit_temp_om_gpp$data %>% 
  distinct(mat_s, log_gpp_s, log_om_s, year, site_id, xmin, xmax) %>%
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~ (1|site_id) + (1|year))

site_medians = site_dots %>% 
  group_by(site_id) %>% 
  median_qi(.epred)

newsimplot = isd_sims %>% 
  ggplot(aes(x = lambda_sims)) +
  # geom_histogram(bins = 100) +
  stat_dots() +
  geom_vline(data = site_medians, aes(xintercept = .epred), 
             # alpha = 0.1,
             color = "#1f968bff",
             linewidth = 0.2) + 
  theme_default() +
  labs(x = "\u03bb (ISD exponent)",
       # subtitle = TeX("$\\lambda = \\frac{log_{10}\\alpha}{log_{10}\\beta} - 0.75 - 1$"),
       y = "") + 
  theme(axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) + 
  annotate(geom = "text", x = -1.5, y = 0.9, label = "Empirical \u03bb's",
           color = "#1f968bff", size = 2.5)  + 
  annotate(geom = "text", x = -2.33, y = 0.3, label = "Simulated \u03bb's",
           color = "grey50", size = 2.5) +
  coord_cartesian(xlim = c(-2.5, -1.3)) +
  scale_y_continuous(limit = c(0, 1),
                     expand = c(0, 0))

saveRDS(newsimplot, file = "plots/fig_s4_newsimplot.rds")
ggsave(newsimplot, width = 5, height = 5, file = "plots/fig_s4_newsimplot.jpg", dpi = 500)
