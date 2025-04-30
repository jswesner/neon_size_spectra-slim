library(tidyverse)
library(ubms)
library(brms)
library(janitor)
library(tidybayes)

#This code ultimately produces dat_2022_clauset.rds
dat_2022_clauset = readRDS(file = "data/dat_2022_clauset.rds")
#Re-running the code may not exactly replicate the dat_2022_clauset due to continuous updates by NEON and any algorithmic stochasticity (e.g., in estimate_xmins)


# Fits a multinomial Poisson depletion model to estimate fish density per collection
# The result is combined in the next script with fish sizes to generate fish size density so that fish can be combined with macroinvertebrates.
streamsites=c("HOPB", "LEWI", "POSE", "CUPE",
              "GUIL", "KING", "MCDI", "LECO",
              "WALK", "MAYF", "ARIK", "BLUE",
              "PRIN", "BLDE", "COMO", "WLOU",
              "SYCA", "REDB", "MART", "MCRA",
              "BIGC", "TECR", "OKSR", "CARI")

# load data (filter out zeros...no fish collected also won't have body sizes)
three_pass_data_wide = read_csv("data/three_pass_data_wide.csv") %>% 
  mutate(total_fish = `1` + `2` + `3`) %>% 
  filter(total_fish > 0) %>% 
  filter(siteID %in% streamsites)

# fit multinomial poisson (three pass depletion model) -----------------------------------------------

# put passes in a matrix (for ubms)
three_pass_matrix = three_pass_data_wide %>% 
  select(`1`,`2`,`3`) %>% 
  as.matrix()

# assign covariates (for ubms)
three_pass_frame <- unmarkedFrameMPois(three_pass_matrix,
                                       siteCovs=as.data.frame(three_pass_data_wide %>% select(site_int)),
                                       type = "removal")
saveRDS(three_pass_frame, file = "data/three_pass_frame.rds")

# fit model
# this estimates population size and capture efficiency for each reach_id (called site_int here)
three_pass_frame = readRDS(file = "data/three_pass_frame.rds")

# three_pass_model = stan_multinomPois(formula = ~1 + (1|site_int) ~ 1 + (1|site_int),
#                                      data = three_pass_frame,
#                                      chains = 4, iter = 2000)
# 
# saveRDS(three_pass_model, file = "models/three_pass_model.rds")

three_pass_model = readRDS(file = "models/three_pass_model.rds")

three_pass_population = as_draws_df(three_pass_model@stanfit) %>% 
  select(contains("_state")) %>% 
  pivot_longer(cols = contains("site_int")) %>% 
  clean_names() %>% 
  mutate(value = beta_state_intercept + value) %>% 
  # bind_rows(sample_1) %>% 
  select(name, value) %>% 
  mutate(site_int = as.factor(parse_number(name))) %>% 
  group_by(site_int) %>% # group and summarize
  median_qi(pop_threepass = exp(value)) %>% # summarize on the probability scale (via link function)
  select(-.width, -.point, -.interval) %>% 
  rename(.lower_threepass = .lower,
         .upper_threepass = .upper) %>% #get original group names
  left_join(three_pass_data_wide %>% ungroup %>% 
              mutate(site_int = as.factor(site_int))) %>% 
  mutate(no_fish_per_m2 = pop_threepass/area_m2,
         no_fish_per_m2_lower = .lower_threepass/area_m2,
         no_fish_per_m2_upper = .upper_threepass/area_m2,
         raw_total_per_m2 = total_fish/area_m2)

three_pass_population %>% 
  ggplot(aes(x = no_fish_per_m2)) + 
  geom_histogram() +
  scale_x_log10()

three_pass_population %>% 
  ggplot(aes(x = raw_total_per_m2, y = no_fish_per_m2)) +
  geom_point() +
  geom_linerange(aes(ymin = no_fish_per_m2_lower, ymax = no_fish_per_m2_upper)) +
  scale_x_log10() +
  scale_y_log10() +
  facet_wrap(~siteID) +
  geom_abline()

# post_fish_density = as_draws_df(three_pass_model@stanfit) %>% 
#   select(contains("_state")) %>% 
#   mutate(.draw = 1:nrow(.)) %>% 
#   pivot_longer(cols = contains("site_int")) %>% 
#   clean_names() %>% 
#   mutate(value = exp(beta_state_intercept + value)) %>%
#   select(name, value, draw) %>% 
#   mutate(site_int = as.factor(parse_number(name))) %>% #get original group names
#   left_join(three_pass_data_wide %>% ungroup %>% 
#               distinct(site_int, reach_id, mean_wetted_width_m, measured_reach_length, collect_date) %>% 
#               mutate(site_int = as.factor(site_int))) %>% 
#   mutate(site_id = str_sub(reach_id, 1, 4)) %>% 
#   group_by(site_id) %>% 
#   mutate(median = median(value)) %>% 
#   group_by(site_int) %>% 
#   mutate(fill_color = median(value),
#          no_fish_per_m2 = value/(mean_wetted_width_m*measured_reach_length))
# 
# 
# fish_density = post_fish_density %>% 
#   group_by(site_int, collect_date, reach_id, site_id) %>%
#   median_qi(no_fish_per_m2) 

saveRDS(three_pass_population, file = "data/fish_density.rds")
