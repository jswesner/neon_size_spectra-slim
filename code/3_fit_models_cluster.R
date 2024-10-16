library(brms)
library(isdbayes)
library(tidyverse)

# This is the code used to fit models on USD's Lawrence Cluster on 10/4/2024

# Instead of re-fitting, you can load the fitted models below
fit_temp_om_gpp = readRDS("models/temporary/fit_temp_om_gpp.rds")

# Re-fit models -----------------------------------------------------------
#1) load data (NEON body size data)
# dat_all = readRDS("data/derived_data/dat_all.rds") %>% 
#  group_by(sample_id) %>% 
#  mutate(xmin = min(dw))

dat_sumnorm = readRDS("data/derived_data/dat_clauset_xmins.rds") 

# # 1) temp*gpp*om
# fit_temp_om_gpp_newxmin_sumnorm = update(fit_temp_om_gpp,
#                                          newdata = dat_sumnorm,
#                                          chains = 4,
#                                          iter = 2000)
# 
# saveRDS(fit_temp_om_gpp_newxmin_sumnorm, file = "models/fit_temp_om_gpp_newxmin_sumnorm_clauset.rds")
# 
# # 2) temp*gpp
# fit_temp_gpp_newxmin_sumnorm = update(fit_temp_om_gpp,
#                                      formula = ~ mat_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
#                                      newdata = dat_sumnorm,
#                                      chains = 4,
#                                      iter = 2000)

# saveRDS(fit_temp_gpp_newxmin_sumnorm, file = "models/fit_temp_gpp_newxmin_sumnorm_clauset.rds")
# 
# 3) om*gpp
#fit_om_gpp_newxmin_sumnorm = update(fit_temp_om_gpp,
#                                    formula = ~ log_om_s*log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
#                                   newdata = dat_sumnorm,
#                                    chains = 4,
#                                    iter = 2000)

#saveRDS(fit_om_gpp_newxmin_sumnorm, file = "models/fit_om_gpp_newxmin_sumnorm_clauset.rds")

# # 4) om*temp
fit_om_temp_newxmin_sumnorm = update(fit_temp_om_gpp,
                                     formula = ~ log_om_s*mat_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
                                     newdata = dat_sumnorm,
                                     chains = 4,
                                     iter = 2000, 
                                     cores = 4)

saveRDS(fit_om_temp_newxmin_sumnorm, file = "models/fit_om_temp_newxmin_sumnorm_clauset.rds")
# # 5) om
fit_om_newxmin_sumnorm = update(fit_temp_om_gpp,
                                formula = ~ log_om_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
                                newdata = dat_sumnorm,
                                chains = 4,
                                iter = 2000, 
                                cores = 4)

saveRDS(fit_om_newxmin_sumnorm, file = "models/fit_om_newxmin_sumnorm_clauset.rds")

# 6) gpp
# fit_gpp_newxmin_sumnorm = update(fit_temp_om_gpp,
#                                 formula = ~ log_gpp_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
#                                  newdata = dat_sumnorm,
#                                  chains = 4,
#                                  iter = 2000)

# saveRDS(fit_gpp_newxmin_sumnorm, file = "models/fit_gpp_newxmin_sumnorm_clauset.rds")

# 7) temp
#fit_temp_newxmin_sumnorm = update(fit_temp_om_gpp,
#                                  formula = ~ mat_s + (1 | sample_id) + (1 | year) + (1 | site_id) ,
#                                  newdata = dat_sumnorm,
#                                  chains = 4,
#                                  iter = 2000)

#saveRDS(fit_temp_newxmin_sumnorm, file = "models/fit_temp_newxmin_sumnorm_clauset.rds")
