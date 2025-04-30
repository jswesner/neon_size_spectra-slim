library(neonUtilities)
library(tidyverse)
library(lubridate)
library(janitor)

source("code/inverts_dw-functions.R")

# load Length Weight coefficient table (used in part C below)
coeff <- read.csv("data/macro_lw_coeffs.csv")
# neon_token <- source("C:/Users/jfpom/Documents/Wesner/NEON documents/neon_token_source.R")
# stream_sites <- readRDS("data/streams.rds")

neon_token <- "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJqZWZmd2VzbmVyQGdtYWlsLmNvbSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc4NjU1MjkxMiwiaWF0IjoxNjI4ODcyOTEyLCJlbWFpbCI6ImplZmZ3ZXNuZXJAZ21haWwuY29tIn0.VnIZyX8yUCBfQyLOtS2hxr_tB4JW2CBzD46QezlxnIKCc1biYv9BbVZvl72obmKP1uXu4iK_c2pzDmBFW_S9oA"

# 1) Download data ---------------------------------------------------------
macro <- loadByProduct(dpID = "DP1.20120.001",
                       site = "ARIK", # set to "all" 
                       startdate = "2016-01",
                       enddate = NA,
                       check.size = FALSE,
                       token = neon_token,
                       nCores = 4,
                       include.provisional = T)

# 2) Add LW coefficients, estimate dry weights  ------------------------------------

# add length weight coefficients by taxon
MN.lw <- LW_coef(x = macro$inv_taxonomyProcessed,
                 lw_coef = coeff,
                 percent = TRUE)

# questionable measurements ####
# filter out individuals that were "damaged" and measurement was affected
# this is a flag which is added by NEON
MN.no.damage <- MN.lw %>%
  filter(!str_detect(sampleCondition,
                     "measurement")) %>%
  est_dw(fieldData = macro$inv_fieldData)

# 3) filter out NA values in dw
macro_dw <- MN.no.damage %>%
  filter(!is.na(dw), !is.na(no_m2)) 

nrow(MN.no.damage) / nrow(macro_dw)

saveRDS(macro_dw, file = "data/macro_dw_raw.rds")

macro_dw = readRDS(file = "data/macro_dw_raw.rds")

# remove taxonomic information and tally density by size class
macro_dw_sizeonly = macro_dw %>% 
  group_by(siteID, collectDate, dw) %>% 
  reframe(no_m2 = mean(no_m2)) %>% 
  group_by(siteID, collectDate, dw) %>% 
  reframe(no_m2 = mean(no_m2))

saveRDS(macro_dw_sizeonly, file = "data/macro_dw_sizeonly.rds")

# 4) filter to only ranges that follow a power law
dat_inverts = macro_dw_sizeonly %>% 
  clean_names() %>% 
  group_by(site_id, collect_date) %>% 
  sample_n(5000, weight = no_m2, replace = T)

dat_inverts_list = dat_inverts %>% group_by(site_id, collect_date) %>% group_split()

xmin_inverts_list = list()

for(i in 1:length(dat_inverts_list)){
  powerlaw = conpl$new(dat_inverts_list[[i]]$dw)
  xmin_inverts_list[[i]] = tibble(xmin_clauset = estimate_xmin(powerlaw)$xmin,
                                  site_id = unique(dat_inverts_list[[i]]$site_id),
                                  collect_date = unique(dat_inverts_list[[i]]$collect_date))
}

xmins_inverts_clauset = bind_rows(xmin_inverts_list)

saveRDS(xmins_inverts_clauset, file = "data/xmins_inverts_clauset.rds")

dat_inverts_clauset_xmins = macro_dw_sizeonly %>% 
  clean_names() %>% 
  left_join(xmins_inverts_clauset) %>% 
  group_by(site_id, collect_date) %>% 
  filter(dw >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset,
         xmax = max(dw))

saveRDS(dat_inverts_clauset_xmins, file = "data/dat_inverts_clauset_xmins.rds")
