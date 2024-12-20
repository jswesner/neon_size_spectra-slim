library(poweRlaw)
library(tidyverse)
# this method follows Clauset et al. 2009 by estimating the minimum
# size for which the data follow a power law using by minimizing the K-S statistic.
# i.e., xmin is model-based

# load data
dat_inverts_all = readRDS("data/derived_data/macro_dw-wrangled.rds")
dat_inverts = dat_inverts_all %>% 
  group_by(sample_id) %>% 
  sample_n(5000, weight = no_m2, replace = T)

dat_inverts_list = dat_inverts %>% group_by(sample_id) %>% group_split()

xmin_inverts_list = list()

for(i in 1:length(dat_inverts_list)){
  powerlaw = conpl$new(dat_inverts_list[[i]]$dw)
  xmin_inverts_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                               sample_id = unique(dat_inverts_list[[i]]$sample_id))
}

xmins_inverts_clauset = bind_rows(xmin_inverts_list)

dat_inverts_clauset_xmins = dat_inverts_all %>% left_join(xmins_inverts_clauset) %>% 
  group_by(sample_id) %>% 
  filter(dw >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset,
         xmax = max(dw))

saveRDS(dat_inverts_clauset_xmins, file = "data/derived_data/dat_inverts_clauset_xmins.rds")

# check cutoffs -----------------------------------------------------------

dat_inverts %>% left_join(dat_inverts_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% 
                         rename(xmin_clauset = xmin)) %>% 
  ggplot(aes(x = dw)) + 
  facet_wrap(~site_id, scales = "free") +
  geom_density(aes(group = sample_id)) +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) +
  NULL

samples = unique(dat_inverts$sample_id)
dat_inverts %>% left_join(dat_inverts_clauset_xmins %>% ungroup %>% distinct(sample_id, xmin) %>% rename(xmin_clauset = xmin)) %>% 
  group_by(sample_id) %>% 
  sample_n(1000, weight = no_m2, replace = T) %>% 
  filter(sample_id == sample(samples, 1)) %>% 
  ggplot(aes(x = dw)) + 
  geom_density() +
  scale_x_log10() +
  geom_vline(aes(xintercept = xmin_clauset)) + 
  facet_wrap(~site_id)
