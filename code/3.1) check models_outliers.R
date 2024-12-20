library(rstan)
library(tidyverse)
library(janitor)
library(tidybayes)
library(brms)
library(ggthemes)
library(isdbayes)
theme_set(brms::theme_default())

# load model
fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_updated12192024.rds")

# load data
sample_size = readRDS("data/sample_size.rds")
dat_clauset_xmins = readRDS(file = "data/derived_data/dat_clauset_xmins.rds")
dat_all = readRDS("data/derived_data/dat_all.rds")

post_medians = fit_temp_om_gpp$data %>% select(-dw, -no_m2) %>% 
  distinct() %>%
  mutate(no_m2 = 1) %>% 
  left_join(sample_size) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = ~ (1|sample_id)) %>% 
  group_by(sample_id, n) %>% 
  median_qi(.epred)


post_medians %>% 
  arrange(.epred) %>% 
  ggplot(aes(x = n, y = .epred)) +
  geom_point()


sample_size_all = fit_temp_om_gpp$data %>% 
  group_by(sample_id) %>% 
  tally(name = "n_fit_temp_om_gpp") %>% 
  left_join(sample_size %>% rename(n_new_data = n)) %>% 
  left_join(dat_all %>% group_by(sample_id) %>% tally(name = "n_original_data"))

dat_clauset_xmins %>% 
  filter(sample_id == 9)


dummy_mod = brm(dw | vreal(no_m2, xmin, xmax) ~ 1,
                stanvars = stanvars,
                family = paretocounts(),
                data = tibble(dw = rparetocounts(),
                              xmin = 1, 
                              xmax = 1000,
                              no_m2 = 1),
                file = "models/dummy_mod.rds")

temp = update(dummy_mod, newdata = dat_clauset_xmins %>% filter(sample_id == 133))
temp2 = update(dummy_mod, newdata = fit_temp_om_gpp$data %>% filter(sample_id == 133))

temp
temp2

dat_1 = fit_temp_om_gpp$data %>% group_by(sample_id) %>% group_split()
dat_2 = dat_clauset_xmins %>% group_by(sample_id) %>% group_split()

mods_1 = list()
mods_2 = list()

for(i in 1:length(dat_1)){
  mods_1[[i]] = update(dummy_mod, iter = 1000, chains = 1, newdata = dat_1[[i]],
                       data2 = list(sample_id = unique(dat_1[[i]]$sample_id),
                                    data = "fit_temp_om_gpp"))
  
  mods_2[[i]] = update(dummy_mod, iter = 1000, chains = 1, newdata = dat_2[[i]],
                       data2 = list(sample_id = unique(dat_1[[i]]$sample_id),
                                    data = "newdata"))
  
}


saveRDS(mods_1, file = "models/old/mods_1.rds")
saveRDS(mods_2, file = "models/old/mods_2.rds")

mods_12 = c(mods_1, mods_2)

posts_12 = list()

for(i in 1:length(mods_12)){
  posts_12[[i]] = as_draws_df(mods_12[[i]]) %>% mutate(sample_id = mods_12[[i]]$data2$sample_id,
                                                  data = mods_12[[i]]$data2$data)
}

bind_rows(posts_12) %>% 
  group_by(sample_id, data) %>% 
  median_qi(b_Intercept) %>% 
  select(sample_id, data, b_Intercept) %>% 
  pivot_wider(names_from = data, values_from = b_Intercept) %>% 
  ggplot(aes(x = fit_temp_om_gpp, y = newdata)) + 
  geom_text(aes(label = sample_id)) +
  geom_abline()


bind_rows(posts_12) %>% 
  group_by(sample_id, data) %>% 
  median_qi(b_Intercept) %>% 
  select(sample_id, data, b_Intercept) %>% 
  pivot_wider(names_from = data, values_from = b_Intercept) %>% 
  left_join(sample_size_all) %>% 
  mutate(diff = newdata - fit_temp_om_gpp) %>% 
  select(-fit_temp_om_gpp, -newdata) %>% 
  pivot_longer(cols = starts_with("n_")) %>% 
  ggplot(aes(x = value, y = diff)) + 
  geom_point() +
  facet_wrap(~name)



