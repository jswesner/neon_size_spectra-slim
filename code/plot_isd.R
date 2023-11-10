plot_isd_posts
get_isd_posts



get_isd_posts(fit_temp_om_gpp, group = c("log_om_s", "mat_s", "log_gpp_s", "site_id", "sample_id", "year"), countname = "no_m2",
              re_formula = NULL)


posts = fit_temp_om_gpp$data %>% distinct(log_om_s, mat_s, log_gpp_s, site_id, sample_id, year) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, 
                  ndraws = 100)


df = fit_temp_om_gpp$data %>% 
  # filter(site_id == "ARIK") %>% 
  group_by(sample_id) %>% 
  add_tally() %>%
  mutate(xmin = min(dw),
         xmax = max(dw)) %>% 
  distinct(sample_id, n, xmin, xmax)


# Create an empty data frame to store the results
result_df <- data.frame()

# Iterate over each row in the original data frame
for (i in 1:nrow(df)) {
  # Generate a log sequence between xmin and xmax for each sample_id
  dw <- 10^seq(log10(df$xmin[i]), log10(df$xmax[i]), length.out = df$n[i])
  
  # Append the results to the new data frame
  result_df <- rbind(result_df, data.frame(sample_id = rep(df$sample_id[i], df$n[i]), dw))
}

# Rename the columns in the result data frame
colnames(result_df) <- c("sample_id", "dw")

# Print the result
as_tibble(result_df) %>% 
  ggplot(aes(x = dw, y = dw)) + 
  geom_point()

# get posts
posts = fit_temp_om_gpp$data %>% 
  group_by(sample_id) %>% 
  mutate(xmin = min(dw), xmax = max(dw)) %>% 
  distinct(log_om_s, log_gpp_s, mat_s, year, site_id, sample_id, xmin, xmax) %>% 
  mutate(no_m2 = 1) %>% 
  add_epred_draws(fit_temp_om_gpp, re_formula = NULL, ndraws = 10)



post_isd = result_df %>%
  left_join(df %>% distinct(sample_id, n)) %>% 
  left_join(posts, relationship = "many-to-many") %>% 
  mutate(prob_yx = (1 - (dw^(.epred + 1) - (xmin^(.epred + 
                                                   1)))/((xmax)^(.epred + 1) - (xmin^(.epred + 1)))), 
         n_yx = prob_yx * n)



post_isd %>% 
  ggplot(aes(x = dw, y = n_yx)) +
  geom_line(aes(group = interaction(sample_id, .draw), 
                color = site_id)) +
  scale_y_log10() +
  scale_x_log10() +
  coord_cartesian(ylim = c(1, NA)) +
  facet_wrap(~site_id) +
  guides(site_id = "none")


isd_data = fit_temp_om_gpp$data %>%
  filter(sample_id == 6) %>%
  mutate(log_weight = log10(no_m2)) %>%
  group_by(log_weight) %>%
  sample_n(size = 1, replace = TRUE) %>%
  ungroup() %>%
  select(-log_weight) %>% 
  arrange(-dw) %>% 
  mutate(n_yx = 1:max(row_number()),
         prob_yx = n_yx/max(n_yx))

post_isd %>% 
  filter(sample_id == 6) %>% 
  ggplot(aes(x = dw, y = prob_yx)) +
  geom_line(aes(group = interaction(sample_id, .draw), 
                color = site_id)) +
  scale_y_log10() +
  scale_x_log10() +
  # coord_cartesian(ylim = c(1, NA)) +
  facet_wrap(~site_id) +
  guides(site_id = "none") +
  geom_point(data = isd_data)


            