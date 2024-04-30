library(tidyverse)
library(sizeSpectra)
# source("scripts/custom_functions.R")
library(lubridate)
library(isdbayes)


d = tibble(dw = rparetocounts(n = 10000, lambda = -1.5, xmin = 0.0001, xmax = 1000000))
dat_all <- readRDS("C:/Users/Jeff.Wesner/OneDrive - The University of South Dakota/USD/Github Projects/neon_size_spectra-slim/data/derived_data/dat_all.rds")


d_list = dat_all %>% group_by(sample_id) %>% group_split()

d_resampled = NULL

for(i in 1:length(d_list)){
  d = d_list[[i]] %>% 
    sample_n(1000, weight = no_m2, replace = T)
  
  breaks = 2^seq(floor(range(log2(d$dw)))[1],
                 ceiling(range(log2(d$dw)))[2])
  
  # bin values using hist()
  binned_hist = hist(d$dw, 
                     breaks = breaks,
                     include.lowest = TRUE, plot = FALSE)              
  
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]                              
  # total bin width = right edge - left edge
  bin_width = breaks_offset - breaks_orig
  count = binned_hist$counts
  log_mids = log10(binned_hist$mids)
  biomass = count * 10**log_mids
  nbiomass = log10(biomass / bin_width)

  dataout = data.frame(
    count = count,
    log_count = log10(count),
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / bin_width),
    # original midpoint of bin log10 transformed
    log_mids = log_mids,
    bin_width = bin_width,
    biomass = biomass,
    nbiomass = nbiomass)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  
  
  breaks2 <- 2^(floor(log2(min(range(d$dw)))):
                  ceiling(log2(max(range(d$dw)))) )
  mid_bin_2 <- log10(breaks2[floor(length(breaks2)/2)]) 
  
  dataout$log_mids_center <- dataout$log_mids - mid_bin_2
  dataout$id = i
  
  d_resampled[[i]] = dataout
  
}


NAS_lms = NULL

for(i in 1:length(d_resampled)){
  NAS_lms[[i]] <- lm(log_count_corrected~log_mids_center,
                       data = d_resampled[[i]])
}

NAS_coefs = NULL

for(i in 1:length(NAS_lms)){
  NAS_coefs[[i]] = tibble(slope_mean = summary(NAS_lms[[i]])$coefficients[2],
                          slope_se = summary(NAS_lms[[i]])$coefficients[4],
                          id = i)
}

id_info = dat_all %>% distinct(sample_id, site_id, xmin, xmax, log_gpp_s, log_om_s, mat_s) %>% 
  ungroup %>% 
  mutate(id = 1:nrow(.))

NAS_lambdas = bind_rows(NAS_coefs) %>% 
  left_join(id_info)

NAS_lambdas %>% 
  ggplot(aes(x = mat_s, y = slope_mean, ymin = slope_mean - slope_se,
             ymax = slope_mean + slope_se)) + 
  geom_pointrange() +
  geom_smooth(method = lm)



