library(tidyverse)
library(viridis)

# simulate parameters
nsims = 1
temp_maineffect = rnorm(nsims, -0.2, 0.) # temp has a generally negative main effect on lambda
resource_maineffect = rnorm(nsims, 0.2, 0.) # resources have average zero main effect on lambda
temp_resource_interaction = rnorm(nsims, 0.05, 0.0) #interaction between temp and resources 

# simulate data
temp = seq(-1, 1, length.out = 20)
resources = seq(-1, 1, length.out = 20)

# simulate lambdas
lambdas = tibble(a1 = temp_maineffect,
       a2 = resource_maineffect,
       a3 = temp_resource_interaction,
       intercept = -1.5) %>% 
  expand_grid(temp = temp) %>% 
  expand_grid(resources = resources) %>% 
  mutate(.epred = intercept + a1*temp + a2*resources + a3*temp*resources)

lambdas %>% 
  ggplot(aes(y = temp, x =resources)) + 
  geom_tile(aes(fill = .epred)) +
  scale_color_viridis() +
  scale_fill_viridis_c(direction = 1, na.value="white") +
  brms::theme_default() +
  labs(fill = "\u03bb",
       x = "Resources",
       y = "Temperature")


