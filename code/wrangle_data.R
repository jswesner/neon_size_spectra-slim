library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based

# load data
dat_all = readRDS("data/derived_data/dat_all.rds")

set.seed(1212124)
dat = dat_all %>% 
  group_by(sample_id) %>% 
  sample_n(5000, weight = no_m2, replace = T)

dat_list = dat %>% group_by(sample_id) %>% group_split()

xmin_list = list()

for(i in 1:length(dat_list)){
  powerlaw = conpl$new(dat_list[[i]]$dw)
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                          sample_id = unique(dat_list[[i]]$sample_id))
}

xmins_clauset = bind_rows(xmin_list)

dat_clauset_xmins = dat_all %>% left_join(xmins_clauset) %>%
  group_by(sample_id) %>%
  filter(dw >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset,
         xmax = max(dw))

saveRDS(dat_clauset_xmins, file = "data/derived_data/dat_clauset_xmins.rds")

dat_clauset_xmins = readRDS(file = "data/derived_data/dat_clauset_xmins.rds")

# check cutoffs -----------------------------------------------------------

dat %>% left_join(dat_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% rename(xmin_clauset = xmin)) %>% 
  ggplot(aes(x = dw)) + 
  facet_wrap(~site_id, scales = "free") +
  geom_density(aes(group = sample_id)) +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) +
  NULL

samples = unique(dat$sample_id)
dat %>% left_join(dat_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% rename(xmin_clauset = xmin)) %>% 
  group_by(sample_id) %>% 
  sample_n(1000, weight = no_m2, replace = T) %>% 
  filter(sample_id == sample(samples, 1)) %>% 
  ggplot(aes(x = dw)) + 
  geom_density() +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) + 
  facet_wrap(~site_id)


# sample size -------------------------------------------------------------

sample_size = dat_clauset_xmins %>% 
  group_by(sample_id) %>% 
  tally() %>% 
  arrange(n)

saveRDS(sample_size, "data/sample_size.rds")

samples = sample_size %>% filter(n<20) %>% pull(sample_id)

dat %>% left_join(dat_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% rename(xmin_clauset = xmin)) %>% 
  filter(sample_id %in% samples) %>% 
  ggplot(aes(x = dw)) + 
  facet_wrap(~site_id, scales = "free") +
  geom_density(aes(group = sample_id)) +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) +
  NULL


# check sensitivity of xmin cutoff
# load data
dat_all = readRDS("data/derived_data/dat_all.rds")

dat = dat_all %>% 
  group_by(sample_id) %>% 
  sample_n(5000, weight = no_m2, replace = T)


dat_list_temp = dat %>% 
  filter(sample_id == 91) %>% 
  group_by(sample_id) %>% group_split()

xmin_list = list()

for(i in 1:length(dat_list_temp)){
  powerlaw = conpl$new(dat_list_temp[[i]]$dw)
  xmin_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                          sample_id = unique(dat_list_temp[[i]]$sample_id))
  print(xmin_list)
}
