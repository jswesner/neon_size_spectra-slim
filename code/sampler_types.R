library(neonUtilities)
library(tidyverse)
library(lubridate)
library(janitor)


# tally sampler types -----------------------------------------------------

# load raw NEON macro data
macro_dw = readRDS(file = "data/macro_dw_raw.rds") 
stream_sites <- read_csv("data/sites.csv")

# tally sampler types
sampler_types_only = macro_dw %>% 
  filter(siteID %in% unique(stream_sites$siteID)) %>% 
  ungroup %>% 
  distinct(sampleID) %>% 
  mutate(sampleID = str_remove(sampleID, "\\.SS")) %>% 
  separate(sampleID, into = c("site", "date", "samplertype", "samplenumber"), 
           remove = F)

sampler_tallies = sampler_types_only %>% 
  group_by(samplertype) %>% 
  tally() %>% 
  arrange(-n)

macro_dw %>%
  # sample_n(10000) %>% 
  left_join(sampler_types_only %>% distinct(sampleID, samplertype)) %>% 
  filter(!is.na(samplertype)) %>%
  filter(siteID %in% unique(stream_sites$siteID)) %>% 
  group_by(samplertype) %>% 
  mutate(median = median(dw*no_m2)) %>% 
  ggplot(aes(x = reorder(samplertype, median), y = dw*no_m2)) + 
  geom_jitter(width = 0.2, shape = ".") +
  # facet_wrap(~siteID) +
  geom_boxplot(aes(group = samplertype),
               outlier.shape = NA) +
  scale_y_log10() +
  coord_flip()

# re-fit by type ----------------------------------------------------------
# load data that uses only surber samplers (and fish)
dat_surber = readRDS("data/dat_2022_clauset_surber.rds")

dat_type = dat_surber %>% 
  group_by(sample_id) %>% 
  group_split()

mod_samplertypes = list()

dummy_mod <- readRDS("models/dummy_mod.rds")

for(i in 1:length(dat_macros_type)){
  mod_samplertypes[[i]] = update(dummy_mod, newdata = dat_type[[i]],
                                 data2 = dat_type[[i]])
}

saveRDS(mod_samplertypes, file = "models/mod_samplertypes.rds")


post_samplertypes_list = list()

for(i in 1:length(mod_samplertypes)){
  post_samplertypes_list[[i]] = as_draws_df(mod_samplertypes[[i]]) %>% 
    mutate(site_id = unique(mod_samplertypes[[i]]$data2$site_id),
           samplertype = "surber",
           sample_id = unique(mod_samplertypes[[i]]$data2$sample_id))
}

post_samplertypes = bind_rows(post_samplertypes_list) %>% 
  group_by(sample_id, samplertype) %>% 
  median_qi(b_Intercept)

fit_macros = readRDS("models/fit_macros.rds")

post_fits = fit_macros$data %>% 
  distinct(sample_id, .keep_all = T) %>% 
  add_epred_draws(fit_macros) %>% 
  group_by(sample_id) %>% 
  median_qi(.epred)

post_fits %>% 
  left_join(post_samplertypes %>% select(-.lower, -.upper)) %>% 
  ggplot(aes(x = b_Intercept, y = .epred)) + 
  geom_point() 



# refit main model --------------------------------------------------------

fit_temp_om_gpp = readRDS("models/fit_temp_om_gpp_year.rds")

fit_temp_om_gpp_surber = update(fit_temp_om_gpp, newdata = dat_surber,
                                iter = 1000, chains = 1)

saveRDS(fit_temp_om_gpp_surber, file = "models/fit_temp_om_gpp_surber.rds")
