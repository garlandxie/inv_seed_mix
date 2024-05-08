# libraries -----
library(dplyr)    # for manipulating data
library(here)     # for creating relative file paths
library(ggplot2)  # for visualizing data 

# import ----

inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))

# exploratory data analysis -----

## check for sample size of density and richness treatments ----
# to determine type of sum of squares (I, II, or III) required for ANOVA 
inv_bm %>%
  group_by(Richness_ID, Density_ID) %>%
  summarize(
    n_avg_bm = length(na.omit(mean_inv_bm_mg)),
    n_tot_bm = length(na.omit(tot_inv_bm_mg)),
    n_inv_ab = length(na.omit(num_invaders))
  )

## check for outliers -----
# reason: outliers increase the estimate of sample variance
# decreasing the calculated F statistics for ANOVA's
# lowering the chance of rejecting the null hypothesis 
inv_bm %>%
  mutate(ID = paste(Richness_ID, Density_ID, Rep)) %>%
  ggplot(aes(x = mean_inv_bm_mg, y = ID)) + 
  geom_point() +
  labs(
    x = "Mean Invader Biomass (mg)",
    y = "ID"
    ) + 
  theme(axis.text.y = element_blank())
  
inv_bm %>%
  mutate(ID = paste(Richness_ID, Density_ID, Rep)) %>%
  ggplot(aes(x = tot_inv_bm_mg, y = ID)) + 
  geom_point() +
  labs(
    x = "Total Invader Biomass (mg)",
    y = "ID"
  ) + 
  theme(axis.text.y = element_blank())

inv_bm %>%
  mutate(ID = paste(Richness_ID, Density_ID, Rep)) %>%
  ggplot(aes(x = num_invaders, y = ID)) + 
  geom_point() +
  labs(
    x = "Invader Abundance",
    y = "ID"
  ) + 
  theme(axis.text.y = element_blank())



