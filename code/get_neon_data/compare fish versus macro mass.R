library(isdbayes)
library(tidybayes)
library(brms)
library(tidyverse)
library(viridis)
library(ggthemes)
library(ggridges)
library(janitor)
library(ggh4x)
library(viridis)

# remake main figures adding analyses with macros or fish only


# 1) remake interaction plot ----------------------------------------------
mod = readRDS("models/fit_temp_om_gpp_cari_oksr_2022_4chain.rds")
mod_macros = readRDS("models/fit_temp_om_gpp_2022_macros.rds")
mod_fish = readRDS("models/fit_temp_om_gpp_2022_fish.rds")

raw_dat = readRDS(file = "data/dat_2024_cari_oksr.rds") %>% 
  left_join(mod$data %>% distinct(sample_id, xmin)) %>% 
  filter(dw >= xmin)

proportion_by_taxon = raw_dat %>% 
  # filter(sample_id == 5) %>% 
  mutate(dw_m2 = dw*no_m2) %>% 
  group_by(taxon, sample_id, site_id) %>%
  reframe(sum = sum(dw_m2)) %>% 
  group_by(sample_id) %>% 
  mutate(total = sum(sum),
         prop = sum/total)

proportion_by_taxon %>%
  filter(taxon != "fish") %>% 
  ggplot(aes(x = prop)) +
  geom_histogram()

proportion_by_taxon %>%
  filter(taxon != "fish") %>% 
  group_by(site_id) %>% 
  mutate(median = median(prop)) %>% 
  ggplot(aes(y = reorder(site_id, median), x = prop)) +
  geom_point() +
  geom_boxplot(aes(group = site_id, fill = median)) +
  scale_fill_viridis()
