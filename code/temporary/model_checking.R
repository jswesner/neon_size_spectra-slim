library(brms)
library(isdbayes)
library(tidyverse)
library(tidybayes)
library(lubridate)


sim_data = tibble(dw = rparetocounts(lambda = -1.5)) %>% 
  mutate(xmin = min(dw),
         xmax = max(dw),
         no_m2 = 1)

fit_test = brm(dw | vreal(no_m2, xmin, xmax) ~ 1,
               data = sim_data, 
               stanvars = stanvars,
               family = paretocounts())

sim_data_2 = sim_data %>% 
  sample_n(100000, weight = no_m2, replace = T)

pp_check(fit_test) +
  scale_x_log10() +
  geom_density(data = sim_data_2, aes(x = dw), color = "red")


pp_data = pp_check(fit_test)


fit_test$data %>% 
  arrange(-dw) %>% 
  mutate(order = 1:nrow(.)) %>% 
  ggplot(aes(x = dw, y = order)) +
  scale_x_log10() +
  scale_y_log10() +
  geom_point()

