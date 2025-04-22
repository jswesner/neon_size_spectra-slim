library(tidyverse)
library(brms)
library(janitor)
library(tidybayes)
library(poweRlaw)

# load predictors. these are standardized after removing SYCA
predictors = readRDS(file = "data/predictors.rds")
dat_2022_macros <- readRDS("data/dat_2024_macros.rds") %>% 
  left_join(predictors) %>% 
  filter(collect_date <= "2022-12-31") %>% 
  filter(!is.na(log_om_s),
         !is.na(log_gpp_s)) %>% 
  mutate(sample_id = macro_id)

# dat_2024 = readRDS(file = "data/dat_2024_fish.rds") %>% 
#   mutate(no_m2 = case_when(taxon == "fish" ~ no_m2*1000,
#                            T ~ no_m2)) # do this for fish to place densities closer to the reach (i.e., the scale of capture)
#                             #. other wise densities are too low to be reliably estimated


get_clauset_xmin = function(data) {
  dat_list = data %>% 
    group_by(sample_id) %>% 
    sample_n(5000, weight = no_m2, replace = TRUE) %>% 
    group_by(sample_id) %>% 
    group_split()
  
  xmin_list = list()
  
  set.seed(202002)
  for(i in seq_along(dat_list)) {
    powerlaw = conpl$new(dat_list[[i]]$dw)
    xmin_value = estimate_xmin(powerlaw)$xmin
    xmin_list[[i]] = tibble(
      xmin_clauset = xmin_value,
      sample_id = unique(dat_list[[i]]$sample_id)
    )
  }
  
  bind_rows(xmin_list)  # Properly combine results
}

cull_xmins = function(data, min_data){
  data %>% 
    left_join(min_data) %>% 
    group_by(sample_id) %>% 
    filter(dw >= xmin_clauset) %>%
    mutate(xmin = xmin_clauset,
           xmax = max(dw)) %>% 
    mutate(year = as.integer(year(collect_date))) %>% 
    left_join(predictors) %>% 
    filter(!is.na(log_om_s)) %>% 
    filter(!is.na(log_gpp_s)) %>% 
    filter(!is.na(mat_s))
  
}

xmins_macros = get_clauset_xmin(data = dat_2022_macros)
dat_2022_macros_clauset = cull_xmins(data = dat_2022_macros, min_data = xmins_macros)

saveRDS(dat_2022_macros_clauset)

dat_2022_clauset_xmins = dat_2024 %>% 
  left_join(xmins_clauset) %>% 
  group_by(sample_id) %>% 
  filter(dw >= xmin_clauset) %>%
  mutate(xmin = xmin_clauset,
         xmax = max(dw)) %>% 
  mutate(year = as.integer(year(collect_date))) %>% 
  left_join(predictors) %>% 
  filter(!is.na(log_om_s)) %>% 
  filter(!is.na(log_gpp_s)) %>% 
  filter(!is.na(mat_s))


# check for low sample sizes
dat_2024_clauset_xmins %>% 
  group_by(sample_id, xmin) %>% 
  tally() %>% 
  arrange(n)

resample = dat_2024_clauset_xmins %>% 
  group_by(sample_id) %>% 
  add_tally() %>% 
  sample_n(1000, weights = no_m2, replace = T) %>% 
  group_by(sample_id, n, xmin, site_id) %>% 
  reframe(gm = exp(mean(log(dw), na.rm = T)),
          sd = sd(dw, na.rm = T)) 

# no apparent funnel for either gm, sd, or xmin, so keep all data
resample %>% 
  ggplot(aes(x = n, y = gm)) +
  # ggplot(aes(x = n, y = sd)) +
  # ggplot(aes(x = n, y = xmin)) +
  geom_point() +
  scale_y_log10() +
  facet_wrap(~site_id)

# save with different year cutoffs
saveRDS(dat_2024_clauset_xmins, file = "data/dat_2024_clauset.rds")
saveRDS(dat_2024_clauset_xmins %>% filter(year <= 2022), 
        file = "data/dat_2022_clauset.rds")
saveRDS(dat_2024_clauset_xmins %>% filter(year >= 2022), 
        file = "data/dat_20232024_clauset.rds")




