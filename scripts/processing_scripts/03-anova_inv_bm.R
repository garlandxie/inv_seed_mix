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

# check box-plots ----
# reason: visualize summary statistics to get a sense of median and spread
# across seeding density and seeding richness treatments
(plot_tot_inv_bm <- inv_bm %>%
   mutate(
     Density_ID = case_when(
       Density_ID == "-"  ~ "20 C. arvense seeds only",
       Density_ID == "D1" ~ "48 native + 20 C. arvense seeds",
       Density_ID == "D2" ~ "200 native + 20 C. arvense seeds", 
       Density_ID == "D3" ~ "400 native + 20 C. arvense seeds") %>%
       factor(levels = c(
         "20 C. arvense seeds only", 
         "48 native + 20 C. arvense seeds", 
         "200 native + 20 C. arvense seeds", 
         "400 native + 20 C. arvense seeds")
       )
   ) %>%
   
   mutate(
     Richness_ID = case_when(
       Density_ID  == "20 invasive seeds"  ~ "Invasive-Only", 
       Richness_ID == "M1" ~ "Native Monoculture", 
       Richness_ID == "M2" ~ "2-species Mixture", 
       Richness_ID == "M4" ~ "4-species Mixture", 
       TRUE ~ Richness_ID
     ), 
     
     Richness_ID = factor(Richness_ID, levels = c(
       "Invasive-Only", 
       "Native Monoculture", 
       "2-species Mixture", 
       "4-species Mixture")
     )
   ) %>%
   ggplot(aes(x = Richness_ID, y = tot_inv_bm_mg)) +
   geom_boxplot(
     aes(fill = Density_ID), 
     position = "dodge2"
   ) + 
   labs(
     fill = "Sown Seed Density", 
     x = NULL, 
     y = "Community-level biomass (mg) of C. arvense"
   ) + 
   theme_bw()
)

(plot_avg_inv_bm <- inv_bm %>%
    mutate(
      Density_ID = case_when(
        Density_ID == "-"  ~ "20 C. arvense seeds only",
        Density_ID == "D1" ~ "48 native + 20 C. arvense seeds",
        Density_ID == "D2" ~ "200 native + 20 C. arvense seeds", 
        Density_ID == "D3" ~ "400 native + 20 C. arvense seeds") %>%
        factor(levels = c(
          "20 C. arvense seeds only", 
          "48 native + 20 C. arvense seeds", 
          "200 native + 20 C. arvense seeds", 
          "400 native + 20 C. arvense seeds")
        )
    ) %>%
    
    mutate(
      Richness_ID = case_when(
        Density_ID  == "20 invasive seeds"  ~ "Invasive-Only", 
        Richness_ID == "M1" ~ "Native Monoculture", 
        Richness_ID == "M2" ~ "2-species Mixture", 
        Richness_ID == "M4" ~ "4-species Mixture", 
        TRUE ~ Richness_ID
      ), 
      
      Richness_ID = factor(Richness_ID, levels = c(
        "Invasive-Only", 
        "Native Monoculture", 
        "2-species Mixture", 
        "4-species Mixture")
      )
    ) %>%
    ggplot(aes(x = Richness_ID, y = mean_inv_bm_mg)) +
    geom_boxplot(
      aes(fill = Density_ID), 
      position = "dodge2"
    ) + 
    labs(
      fill = "Sown Seed Density", 
      x = NULL, 
      y = "Average biomass (mg) of C. arvense"
    ) + 
    theme_bw()
)

(plot_ab_inv_bm <- inv_bm %>%
    mutate(
      Density_ID = case_when(
        Density_ID == "-"  ~ "20 C. arvense seeds only",
        Density_ID == "D1" ~ "48 native + 20 C. arvense seeds",
        Density_ID == "D2" ~ "200 native + 20 C. arvense seeds", 
        Density_ID == "D3" ~ "400 native + 20 C. arvense seeds") %>%
        factor(levels = c(
          "20 C. arvense seeds only", 
          "48 native + 20 C. arvense seeds", 
          "200 native + 20 C. arvense seeds", 
          "400 native + 20 C. arvense seeds")
        )
    ) %>%
    
    mutate(
      Richness_ID = case_when(
        Density_ID  == "20 invasive seeds"  ~ "Invasive-Only", 
        Richness_ID == "M1" ~ "Native Monoculture", 
        Richness_ID == "M2" ~ "2-species Mixture", 
        Richness_ID == "M4" ~ "4-species Mixture", 
        TRUE ~ Richness_ID
      ), 
      
      Richness_ID = factor(Richness_ID, levels = c(
        "Invasive-Only", 
        "Native Monoculture", 
        "2-species Mixture", 
        "4-species Mixture")
      )
    ) %>%
    ggplot(aes(x = Richness_ID, y = num_invaders)) +
    geom_boxplot(
      aes(fill = Density_ID), 
      position = "dodge2"
    ) + 
    labs(
      fill = "Sown Seed Density", 
      x = NULL, 
      y = "Average Invader Abundance"
    ) + 
    theme_bw()
)


