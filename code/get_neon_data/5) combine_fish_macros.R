library(tidyverse)
library(brms)
library(janitor)
library(tidybayes)
library(poweRlaw)

#This code ultimately produces dat_2022_clauset.rds
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
#Re-running the code may not exactly replicate the dat_2022_clauset due to continuous updates by NEON and any algorithmic stochasticity (e.g., in estimate_xmins)

# 1) load raw sizes (not culled for xmin or anything)
macro_dw_sizeonly_temp = readRDS(file = "data/macro_dw_sizeonly.rds") %>% mutate(taxon = "macro",
                                                                            year = year(collectDate),
                                                                            month = month(collectDate),
                                                                            macrodate = ymd(as.Date(collectDate))) %>% 
  filter(!is.na(dw)) %>% 
  filter(!is.na(no_m2)) %>% 
  clean_names() %>% 
  group_by(year, month, site_id) %>% 
  mutate(macro_id = cur_group_id(),
         dateplus = collect_date + days(30),
         dateminus = collect_date - days(30)) %>% 
  mutate(collect_date = as.Date(collect_date),
         dateplus = as.Date(dateplus),
         dateminus = as.Date(dateminus))

fish_dw_sizeonly_temp = readRDS(file = "data/fish_dw_sizeonly.rds") %>% 
  rename(collectDate = boutEndDate) %>% 
  mutate(taxon = "fish",
         year = year(collectDate),
         month = month(collectDate),
         fishdate = as.Date(collectDate)) %>% 
  filter(!is.na(dw)) %>%
  filter(!is.na(no_m2)) %>% 
  clean_names() %>% 
  group_by(year, month, site_id) %>% 
  mutate(fish_id = cur_group_id())


# 2) filter to dates within X days of macro collection

fish = fish_dw_sizeonly_temp %>% ungroup %>% select(site_id, collect_date, dw, no_m2, taxon)
macros = macro_dw_sizeonly_temp %>% ungroup %>% select(site_id, collect_date, dw, no_m2, taxon)

data = bind_rows(fish, macros)

# Initialize id column
data$id <- 0

current_id <- 1
for (i in seq_len(nrow(data))) {
  if (data$id[i] == 0) {
    data$id[i] <- current_id
    within_30_days <- which(abs(difftime(data$collect_date[i], data$collect_date, units = "days")) <= 30 & data$id == 0)
    data$id[within_30_days] <- current_id
    current_id <- current_id + 1
  }
}

data2 = data %>% group_by(id, site_id) %>% 
  mutate(sample_id = cur_group_id())

ids_to_keep = data2 %>% 
  distinct(taxon, id, sample_id, site_id) %>% 
  group_by(sample_id, site_id) %>% 
  tally() %>% 
  mutate(n = case_when(site_id %in% c("CARI", "OKSR") ~ 2, TRUE ~ n)) %>% 
  filter(n > 1)

macro_fish_raw = data2 %>% 
  filter(sample_id %in% ids_to_keep$sample_id) %>% 
  group_by(site_id, dw, sample_id) %>% 
  mutate(collect_date == min(collect_date))

saveRDS(macro_fish_raw, file = "data/macro_fish_raw.rds")

dat_2024 = macro_fish_raw %>% 
  filter(sample_id %in% ids_to_keep$sample_id) %>% 
  group_by(site_id, dw, sample_id) %>% 
  mutate(collect_date = min(collect_date)) %>% 
  group_by(site_id, collect_date, dw, sample_id) %>% 
  reframe(no_m2 = mean(no_m2))

saveRDS(dat_2024, file = "data/dat_2024.rds")

dat_2024_macros = macro_dw_sizeonly_temp %>% 
  group_by(site_id, dw, macro_id) %>% 
  mutate(collect_date == min(collect_date)) %>%
  group_by(site_id, collect_date, dw, macro_id) %>% 
  reframe(no_m2 = mean(no_m2))

dat_2024_fish = fish %>% 
  group_by(collect_date, site_id) %>% 
  mutate(sample_id = cur_group_id()) %>% 
  group_by(site_id, dw, sample_id) %>% 
  mutate(collect_date == min(collect_date)) %>% 
  group_by(site_id, collect_date, dw, sample_id) %>% 
  reframe(no_m2 = mean(no_m2)) %>% 
  mutate(taxon = "fish")

saveRDS(dat_2024_macros, file = "data/dat_2024_macros.rds")
saveRDS(dat_2024_fish, file = "data/dat_2024_fish.rds")




