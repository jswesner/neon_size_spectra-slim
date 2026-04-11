library(tidyverse)

dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")

sites = dat_2022_clauset %>% distinct(site_id) %>% pull()
dates = dat_2022_clauset %>% distinct(collect_date) %>% pull()

# get raw, pre-culled data
fish =  readRDS(file = "data/raw_data/fish.rds")

fish_raw <- fish$fsh_perFish %>% 
  mutate(dates = passStartTime) %>% 
  clean_names() %>% 
  filter(site_id %in% sites) %>% 
  filter(pass_end_time <= "2023-01-01") %>% 
  as_tibble() %>% 
  mutate(year = year(pass_end_time),
         month = month(pass_end_time),
         year_month = paste0(year, "_", month))

macro_dw_raw <- readRDS("data/macro_dw_raw.rds") %>% 
  clean_names() %>% 
  mutate(dates = collect_date) %>% 
  filter(site_id %in% sites) %>% 
  filter(collect_date <= "2023-01-01") %>% 
  as_tibble() %>% 
  mutate(year = year(collect_date),
         month = month(collect_date),
         year_month = paste0(year, "_", month))


# total macros
total_macros = sum(macro_dw_raw$individual_count)

# total fish
total_fish = nrow(fish_raw %>% 
  filter(!is.na(fish_weight)))

# total sizes
total_macros + total_fish

# total invert species
length(unique(macro_dw_raw$scientific_name))

# total fish species
length(unique(fish_raw$taxon_id))

# total collections 
length(unique(dat_2022_clauset$sample_id))

# total streams 
length(unique(dat_2022_clauset$site_id))

# size range
log10(max(dat_2022_clauset$dw)) - log10(min(macro_dw_raw$dw))


# totals after culling ----------------------------------------------------

xmins_clauset = readRDS("data/xmins_clauset.rds") 

# assume that we delete body sizes below an average xmin at a site. We have to assume this b/c
# there is no way to match the derived samples with fish/inverts to the original sample ids (at least not easily)
site_xmins = xmins_clauset %>% 
  group_by(site_id, year) %>%
  reframe(xmin = min(xmin_clauset)) %>% 
  mutate(min_xmin = min(xmin)) 

# total macros
total_macros_afterculling = macro_dw_raw %>% 
  left_join(site_xmins) %>% 
  filter(dw >= min_xmin) %>% ungroup %>% 
  reframe(sum = sum(individual_count))

# total fish
total_fish_afterculling = fish_raw %>% 
                    filter(!is.na(fish_weight)) %>% 
  left_join(site_xmins) %>% 
  filter(fish_weight >= min_xmin) %>% 
  nrow()

# total sizes
total_macros_afterculling$sum + total_fish_afterculling
