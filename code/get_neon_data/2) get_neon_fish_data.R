library(neonUtilities)
library(tidyverse)
library(lubridate)
library(janitor)

#This code ultimately produces dat_2022_clauset.rds
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
#Re-running the code may not exactly replicate the dat_2022_clauset due to continuous updates by NEON and any algorithmic stochasticity (e.g., in estimate_xmins)

neon_token <- "eyJ0eXAiOiJKV1QiLCJhbGciOiJFUzI1NiJ9.eyJhdWQiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnL2FwaS92MC8iLCJzdWIiOiJqZWZmd2VzbmVyQGdtYWlsLmNvbSIsInNjb3BlIjoicmF0ZTpwdWJsaWMiLCJpc3MiOiJodHRwczovL2RhdGEubmVvbnNjaWVuY2Uub3JnLyIsImV4cCI6MTc4NjU1MjkxMiwiaWF0IjoxNjI4ODcyOTEyLCJlbWFpbCI6ImplZmZ3ZXNuZXJAZ21haWwuY29tIn0.VnIZyX8yUCBfQyLOtS2hxr_tB4JW2CBzD46QezlxnIKCc1biYv9BbVZvl72obmKP1uXu4iK_c2pzDmBFW_S9oA"

# 1) Download data ---------------------------------------------------------
# assign("has_internet_via_proxy", TRUE, environment(curl::has_internet))
# 
# fish <- loadByProduct(dpID = "DP1.20107.001",
#                        # site = "PRLA", 
#                        startdate = "2016-01",
#                        enddate = NA,
#                        check.size = FALSE,
#                        token = neon_token,
#                        nCores = 4,
#                        include.provisional = T)
# 
# saveRDS(fish, file = "data/raw_data/fish.rds")

fish = readRDS(file = "data/raw_data/fish.rds")


fsh_bulkCount = fish$fsh_bulkCount
fsh_fieldData = fish$fsh_fieldData
fsh_perPass = fish$fsh_perPass
fsh_perFish = fish$fsh_perFish
mean_wetted_width = readRDS(file = "data/mean_wetted_width.rds")


# Derive a L1 version of the passTaxaCount table
bulk_for_passTaxaCount <- fsh_bulkCount %>% select(eventID, boutEndDate, passNumber, taxonID, namedLocation, bulkFishCount) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  summarize(bulk_sum = sum(bulkFishCount))

pass_taxa_count_table <- fsh_perFish %>% select(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation, taxonID) %>%
  tally() %>%
  rename(fish_sum = n) %>%
  left_join(., bulk_for_passTaxaCount, by = c("eventID", "boutEndDate", "passNumber", "taxonID", "namedLocation")) |>
  mutate(bulk_sum = ifelse(is.na(bulk_sum), 0, bulk_sum),
         taxa_total = bulk_sum + fish_sum) |>
  ungroup()|>
  select(-bulk_sum, -fish_sum)


# Derive a L1 version of the fishPerPassSum in the perPass table
bulk_perPass_fishSums <- fsh_bulkCount %>% select(eventID, boutEndDate, passNumber, namedLocation, bulkFishCount) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation) %>%
  summarize(bulk_sum = sum(bulkFishCount))

perPassFishSums <- fsh_perFish %>% select(eventID, boutEndDate, passNumber, namedLocation) %>%
  group_by(eventID, boutEndDate, passNumber, namedLocation) %>%
  tally() %>%
  rename(fish_sum = n) %>%
  left_join(., bulk_perPass_fishSums, by = c("eventID", "boutEndDate", "passNumber", "namedLocation")) %>%
  mutate(bulk_sum = ifelse(is.na(bulk_sum), 0, bulk_sum),
         fish_total = bulk_sum + fish_sum) |>
  ungroup()|>
  select(-bulk_sum, -fish_sum) %>%
  left_join(fsh_perPass, ., by = c("eventID", "boutEndDate", "passNumber", "namedLocation")) 

# Make wide and filter false zeros
true_zeros = fsh_perPass %>%
  select(siteID, eventID, boutEndDate, passNumber, namedLocation, targetTaxaPresent) %>%
  distinct()  %>% # removes duplicates. JSW confirmed that these were true duplicates on 2023-03-01
  filter(targetTaxaPresent == "N")

# 8) make wide format. Replace 0's using information in fsh_perPass$target_taxa_present
three_pass_data_wide = perPassFishSums %>% 
  group_by(siteID, domainID, namedLocation, eventID, boutEndDate, passNumber) %>%
  reframe(fish_total = sum(fish_total, na.rm = T)) %>% 
  bind_rows(true_zeros) %>% 
  filter(passNumber <= 3) %>%
  pivot_wider(names_from = passNumber, values_from = fish_total) %>%
  replace_na(list(`1` = 0, # replace NA with zeros (assumes zero fish if there were no values entered)
                  `2` = 0,
                  `3` = 0)) %>%
  mutate(`1` = case_when(`1` == 0 & is.na(targetTaxaPresent) ~ 1e9,   # create silly number to filter out false zeros
                         TRUE ~ `1`),
         `2` = case_when(`2` == 0 & is.na(targetTaxaPresent) ~ 1e9,
                         TRUE ~ `2`),
         `3` = case_when(`3` == 0 & is.na(targetTaxaPresent) ~ 1e9,
                         TRUE ~ `3`)) %>%
  filter(`1` < 1e9) %>% # filter out false zeros
  filter(`2` < 1e9) %>%
  filter(`3` < 1e9) %>%
  mutate(site_int = as.factor(row_number()), #sample identifier
         increased = case_when(`3` > `1` ~ "no depletion",
                               TRUE ~ "depletion")) %>% 
  filter(!is.na(boutEndDate)) %>% 
  left_join(fsh_fieldData %>% glimpse() %>% distinct(siteID, namedLocation, eventID, 
                                                     measuredReachLength) %>% 
              filter(!is.na(measuredReachLength))) %>% # two collectons in MAYF had different reach lengths. The next pipe removes them.
  group_by(site_int) %>% 
  add_tally() %>% 
  filter(n == 1) %>% 
  left_join(mean_wetted_width %>% rename(siteID = site_id)) %>% 
  mutate(area_m2 = mean_wetted_width_m*measuredReachLength)

# write_csv(three_pass_data_wide, file = "data/three_pass_data_wide.csv")


# Derive a L1 version of the fishfieldDataSum in the fieldData table
bulk_fieldData_fishSums <- fsh_bulkCount %>% select(eventID, boutEndDate, namedLocation, bulkFishCount) %>%
  group_by(eventID, boutEndDate, namedLocation) %>%
  summarize(bulk_sum = sum(bulkFishCount))

fieldDataFishSums <- fsh_perFish %>% select(eventID, boutEndDate, namedLocation) %>%
  group_by(eventID, boutEndDate, namedLocation) %>%
  tally() %>%
  rename(fish_sum = n) %>%
  left_join(., bulk_fieldData_fishSums, by = c("eventID", "boutEndDate", "namedLocation")) %>%
  mutate(bulk_sum = ifelse(is.na(bulk_sum), 0, bulk_sum),
         fish_total = bulk_sum + fish_sum) |>
  ungroup()|>
  select(-bulk_sum, -fish_sum) %>%
  left_join(fsh_fieldData, ., by = c("eventID", "boutEndDate", "namedLocation"))

