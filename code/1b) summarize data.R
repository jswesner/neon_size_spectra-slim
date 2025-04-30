library(tidyverse)

dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")

sites = dat_2022_clauset %>% distinct(site_id) %>% pull()
dates = dat_2022_clauset %>% distinct(collect_date) %>% pull()

# get raw, pre-culled data
fish =  readRDS(file = "data/raw_data/fish.rds")

fish_raw <- fish$fsh_perFish %>% 
  glimpse() %>% 
  mutate(dates = passStartTime) %>% 
  clean_names() %>% 
  filter(site_id %in% sites) %>% 
  filter(pass_end_time <= "2023-01-01") %>% 
  as_tibble()

macro_dw_raw <- readRDS("data/macro_dw_raw.rds") %>% 
  clean_names() %>% 
  mutate(dates = collect_date) %>% 
  filter(site_id %in% sites) %>% 
  filter(collect_date <= "2023-01-01") %>% 
  as_tibble()


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
