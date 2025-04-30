library(tidyverse)


# load data ---------------------------------------------------------------
# individual body sizes
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
# The body size data is wrangled using six scripts in the folder get_neon_data
# Those scripts take several hours to run, mostly due to large raw downloads of NEON data and fitting fish density models
# The end product is loaded above. It is the data used in the isd models for this manuscript.

# predictors (backtransformed. These equations were verified most recently by jsw on 04/27/2025. They recapture the raw values of temp, om, and gpp
# that were first calculated ~2022)
predictors_scaled = readRDS("data/predictors.rds") %>% 
  mutate(temp_deg_c = mat_s*attr(predictors$mat_s, 'scaled:scale') + attr(predictors$mat_s, 'scaled:center'),
         om = exp(log_om_s*attr(predictors$log_om_s, 'scaled:scale') + attr(predictors$log_om_s, 'scaled:center')),
         gpp = exp(log_gpp_s*attr(predictors$log_gpp_s, 'scaled:scale') + attr(predictors$log_gpp_s, 'scaled:center')))


saveRDS(predictors_scaled, file = "data/predictors_scaled.rds")
