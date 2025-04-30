library(tidyverse)
library(brms)
library(janitor)
library(tidybayes)

#This code ultimately produces dat_2022_clauset.rds
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
#Re-running the code may not exactly replicate the dat_2022_clauset due to continuous updates by NEON and any algorithmic stochasticity (e.g., in estimate_xmins)


# 1) load fish density posteriors (total number of fish per m2 per collection)
fish_density = readRDS(file = "data/fish_density.rds") %>% as_tibble()
fish = readRDS(file = "data/raw_data/fish.rds") # NEON FISH DATA list

# collection IDs
site_ints = fish_density %>% ungroup %>% distinct(site_int, siteID, domainID, namedLocation, eventID, boutEndDate)

# 2) load fish sizes. Add collection ids
fsh_perFish = fish$fsh_perFish %>%  # sizes
  left_join(site_ints) %>% 
  filter(!is.na(site_int))

# 3) Get proportion of fish sizes per collection id
fish_sizes = fsh_perFish %>%
  mutate(dw_g = 0.2*fishWeight, # convert wet to dry mass. From Brey, T., Müller-Wiegmann, C., Zittier, Z. M., & Hagen, W. (2010). Body composition in aquatic organisms—a global data bank of relationships between mass, elemental composition and energy content. Journal of Sea Research, 64(3), 334-340.
         dw = dw_g*1000) %>% ##### units are in mgDM
  select(-fishWeight) %>% 
  group_by(site_int) %>% 
  add_tally(name = "total") %>% 
  group_by(site_int, siteID, domainID, namedLocation, eventID, boutEndDate, dw, total) %>% 
  tally() %>% 
  mutate(proportion_dw = n/total) 

# 3) merge and calculate size density (number of body sizes per m2 per collection)

fish_dw_sizeonly = fish_density %>% 
  select(site_int, no_fish_per_m2, no_fish_per_m2_lower, no_fish_per_m2_upper) %>% 
  left_join(fish_sizes, by = "site_int") %>% 
  mutate(no_m2 = proportion_dw*no_fish_per_m2) %>% 
  select(site_int, siteID, boutEndDate, dw, no_m2, starts_with("no_fish_per"))

saveRDS(fish_dw_sizeonly, file = "data/fish_dw_sizeonly.rds")

# check
fish_dw_sizeonly %>% 
  filter(siteID == "ARIK") %>% 
  group_by(site_int) %>% 
  sample_n(2000, replace = T, weight = no_m2) %>% 
  group_by(site_int) %>% 
  arrange(-dw) %>% 
  mutate(order = 1:max(row_number())) %>% 
  ggplot(aes(x = dw, y = order)) + 
  geom_point() +
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(~site_int)

max(fish_dw_sizeonly$dw, na.rm = T) # should be ~2e+05
