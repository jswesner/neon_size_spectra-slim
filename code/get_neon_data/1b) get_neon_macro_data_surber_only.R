library(neonUtilities)
library(tidyverse)
library(lubridate)
library(janitor)

#This code ultimately produces dat_2022_clauset.rds
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
#Re-running the code may not exactly replicate the dat_2022_clauset due to continuous updates by NEON and any algorithmic stochasticity (e.g., in estimate_xmins)


source("code/get_neon_data/inverts_dw-functions.R")

# load Length Weight coefficient table (used in part C below)
coeff <- read.csv("data/macro_lw_coeffs.csv")
# stream_sites <- readRDS("data/streams.rds")

neon_token = readRDS("data/neon_token.rds") #!!!Users will need their own token.

# 1) Download data ---------------------------------------------------------
# assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
# 
# macro <- loadByProduct(dpID = "DP1.20120.001",
#                        # site = "PRLA", 
#                        startdate = "2016-01",
#                        enddate = NA,
#                        check.size = FALSE,
#                        token = neon_token,
#                        nCores = 4,
#                        include.provisional = T)
# 
# saveRDS(macro, file = "data/raw_data/macro.rds")

macro = readRDS(file = "data/raw_data/macro.rds")

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

stream_sites <- read_csv("data/sites.csv")

macro_dw = readRDS(file = "data/macro_dw_raw.rds") 

# remove taxonomic information and tally density by size class
macro_dw_sizeonly_surber = macro_dw %>% 
  mutate(sampleID = str_remove(sampleID, "\\.SS")) %>% 
  separate(sampleID, into = c("site", "date", "samplertype", "samplenumber"), 
           remove = F) %>% 
  filter(samplertype == "SURBER") %>% 
  group_by(siteID, collectDate, dw) %>% 
  reframe(no_m2 = mean(no_m2))

saveRDS(macro_dw_sizeonly_surber, file = "data/macro_dw_sizeonly_surber.rds")


