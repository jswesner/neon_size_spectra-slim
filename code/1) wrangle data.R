library(tidyverse)
library(brms)
library(janitor)
library(tidybayes)
library(poweRlaw)

# this code combines macro and fish body size data sets using a fuzzy join for collections within 30 days of each other.
# The original data loaded below was downloaded using the get_neon_body_sizes repo.https://github.com/jswesner/get_neon_body_sizes

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
  distinct(taxon, id, sample_id) %>% 
  group_by(sample_id) %>% 
  tally() %>% 
  filter(n > 1)

macro_fish_raw = data2 %>% 
  filter(sample_id %in% ids_to_keep$sample_id) %>% 
  group_by(site_id, dw, sample_id, no_m2, taxon) %>% 
  mutate(collect_date == min(collect_date))

saveRDS(macro_fish_raw, file = "data/macro_fish_raw.rds")


dat_2024 = macro_fish_raw %>% 
  filter(sample_id %in% ids_to_keep$sample_id) %>% 
  group_by(site_id, dw, sample_id) %>% 
  mutate(collect_date == min(collect_date)) %>% 
  group_by(site_id, collect_date, dw, sample_id) %>% 
  reframe(no_m2 = mean(no_m2))

saveRDS(dat_2024, file = "data/dat_2024.rds")

# don't exclude fishless at CARI and OKSR ---------------------------------

ids_to_keep2 = data2 %>% 
  distinct(taxon, id, sample_id, site_id) %>% 
  group_by(sample_id, site_id) %>% 
  tally() %>% 
  mutate(n = case_when(site_id %in% c("CARI", "OKSR") ~ 2, TRUE ~ n)) %>% 
  filter(n > 1)

dat_2024_cari_oksr = data2 %>% 
  filter(sample_id %in% ids_to_keep2$sample_id) %>% 
  group_by(site_id, dw, sample_id) %>% 
  mutate(collect_date == min(collect_date)) %>%
  group_by(site_id, collect_date, dw, sample_id, taxon) %>% 
  reframe(no_m2 = mean(no_m2))

saveRDS(dat_2024_cari_oksr, file = "data/dat_2024_cari_oksr.rds")




