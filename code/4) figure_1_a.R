library(isdbayes)
library(tidybayes)
library(tidyverse)
library(brms)
library(viridis)
theme_set(theme_default())

sim_1 = tibble(x = rparetocounts(n = 3000, lambda = -1, xmin = 1, xmax = 1000)) %>% 
  mutate(xmin = 1, xmax = 1000, counts = 1, group = "-1")
sim_2 = tibble(x = rparetocounts(n = 3000, lambda = -1.5, xmin = 1, xmax = 1000)) %>% 
  mutate(xmin = 1, xmax = 1000, counts = 1, group = "-1.5")
sim_3 = tibble(x = rparetocounts(n = 3000, lambda = -2, xmin = 1, xmax = 1000)) %>% 
  mutate(xmin = 1, xmax = 1000, counts = 1, group = "-2")

sim_data = bind_rows(sim_1, sim_2, sim_3) %>% 
  group_by(group) %>% 
  arrange(-x) %>% 
  mutate(order = row_number())

sim_mod = readRDS(file = "models/sim_mod.rds")

# sim_mod = brm(x|vreal(counts, xmin, xmax) ~ group,
#               stanvars = stanvars,
#               family = paretocounts(),
#               chains = 1, iter = 1000,
#               data = sim_data)
# 
# saveRDS(sim_mod, file = "models/sim_mod.rds")

# sim_mod = update(sim_mod, newdata = sim_data)

lineposts = tibble(x = seq(1, 1000, length.out = 50)) %>% 
  mutate(xmin = 1, xmax = 1000, counts = 1) %>% 
  expand_grid(group = unique(sim_data$group)) %>% 
  add_epred_draws(sim_mod) %>% 
  group_by(x, xmin, xmax, group) %>% 
  median_qi(.epred) %>% 
  mutate(prob_yx = (1 - (x^(.epred + 1) - (xmin^(.epred + 1))) / ((xmax)^(.epred + 1) - (xmin^(.epred + 1)))),
         dw_mg = x,
         prob_yx_lower = (1 - (x^(.lower + 1) - (xmin^(.lower + 1))) / ((xmax)^(.lower + 1) - (xmin^(.lower + 1)))),
         prob_yx_upper = (1 - (x^(.upper + 1) - (xmin^(.upper + 1))) / ((xmax)^(.upper + 1) - (xmin^(.upper + 1))))) 

labels = tibble(label = c("\u03bb = -2",
                          "\u03bb = -1.5",
                          "\u03bb = -1"),
                y = c(20, 200, 1500),
                group = c("-2", "-1.5", "-1")) %>% 
  mutate(x = 100)

fig1a = lineposts %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = prob_yx*3000 , group = group)) +
  scale_x_log10() + 
  scale_y_log10() +
  geom_point(data = sim_data, aes(y = order, size = x, color = x), shape = 1) +
  coord_cartesian(ylim = c(1, NA)) +
  labs(x = "mgDM Individual",
       y = "Number of individuals >= x") +
  scale_color_viridis(option = "E") +
  guides(color = "none",
         size = "none") +
  scale_size(range=c(0.1, 2.5)) +
  geom_label(aes(x = x, y = y, label = label),
             data = labels,
             size = 2.5)

# ggview::ggview(fig1a, width = 3.5, height = 3.5)
ggsave(fig1a, file = "plots/fig1a.jpg", width = 3.5, height = 3.5, dpi = 500)


lineposts %>% 
  ggplot(aes(x = x)) +
  geom_line(aes(y = prob_yx*3000 , group = group)) +
  scale_x_log10() + 
  # scale_y_log10() +
  geom_point(data = sim_data, aes(y = order, size = x, color = x), shape = 1) +
  coord_cartesian(ylim = c(1, NA)) +
  labs(x = "mgDM Individual",
       y = "Number of individuals >= x") +
  scale_color_viridis(option = "E") +
  guides(color = "none",
         size = "none") +
  scale_size(range=c(0.1, 2.5)) +
  geom_label(aes(x = x, y = y, label = label),
             data = labels,
             size = 2.5) +
  facet_wrap(~group)

