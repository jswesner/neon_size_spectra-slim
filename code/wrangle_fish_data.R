library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based

# load data
dat_fish_all = readRDS("data/derived_data/dat_all.rds")
dat_fish = dat_fish_all %>% 
  group_by(sample_id) %>% 
  sample_n(5000, weight = no_m2, replace = T)

dat_fish_list = dat_fish %>% group_by(sample_id) %>% group_split()

xmin_fish_list = list()

for(i in 1:length(dat_fish_list)){
  powerlaw = conpl$new(dat_fish_list[[i]]$dw)
  xmin_fish_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                          sample_id = unique(dat_fish_list[[i]]$sample_id))
}

xmins_fish_clauset = bind_rows(xmin_fish_list)

dat_fish_clauset_xmins = dat_fish_all %>% left_join(xmins_fish_clauset) %>% 
  group_by(sample_id) %>% 
  filter(dw >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset,
         xmax = max(dw))

saveRDS(dat_fish_clauset_xmins, file = "data/derived_data/dat_fish_clauset_xmins.rds")


# check cutoffs -----------------------------------------------------------

dat_fish %>% left_join(dat_fish_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% 
                         rename(xmin_clauset = xmin)) %>% 
  ggplot(aes(x = dw)) + 
  facet_wrap(~site_id, scales = "free") +
  geom_density(aes(group = sample_id)) +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) +
  NULL

samples = unique(dat_fish$sample_id)
dat_fish %>% left_join(dat__fishclauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% rename(xmin_clauset = xmin)) %>% 
  group_by(sample_id) %>% 
  sample_n(1000, weight = no_m2, replace = T) %>% 
  filter(sample_id == sample(samples, 1)) %>% 
  ggplot(aes(x = dw)) + 
  geom_density() +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) + 
  facet_wrap(~site_id)
