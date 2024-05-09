# libraries -----
library(dplyr)    # for manipulating data
library(here)     # for creating relative file paths
library(ggplot2)  # for visualizing data 

# import ----

inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))

# exploratory data analysis -----

## outliers -----
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

## box-plots ----
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

# statistical analyses ----

# change the values for invader monocultures
# so that you can compare native monocultures + mixtures for the ANOVAs
inv_bm_tidy <- inv_bm %>%
  
  mutate(
    Density_ID = case_when(
      Density_ID == "-" ~ "D0", 
      TRUE ~ Density_ID),
    
    Richness_ID = case_when(
      Richness_ID == "M1" & Density_ID == "D0" ~ "Inv", 
      TRUE ~ Richness_ID)
    ) %>%

  mutate(
    Richness_ID = factor(
      Richness_ID, 
      levels = c("Inv", "M1", "M2", "M4")),
    
    Density_ID = factor(
      Density_ID, 
      levels = c("D0", "D1", "D2", "D3")
    )
  ) %>%
  
  rename(
    Density = Density_ID, 
    Richness = Richness_ID
  )

## model fit ----

lm_tot_inv_bm1 <- lm(
  tot_inv_bm_mg ~ Richness*Density, 
  data = inv_bm_tidy)

# check which regression coefficients are due to singularities
# check for highly collinear terms 
summary(lm_tot_inv_bm1)
alias(lm_tot_inv_bm1) # look for complete aliasing

# drop invader monocultures to avoid aliased coefficients
inv_bm_minus_d0 <- inv_bm_tidy %>%
  mutate(Density = as.character(Density)) %>%
  dplyr::filter(Density %in% c("D1", "D2", "D3")) 

lm_inv_bm_minus_d0 <- lm(
  tot_inv_bm_mg ~ Richness*Density, 
  data = inv_bm_minus_d0
  )

# check assumptions for two-way ANOVA
plot(lm_inv_bm_minus_d0, 2) # normality of residuals
plot(lm_inv_bm_minus_d0, 1) # homogeneity of variance
plot(lm_inv_bm_minus_d0, 4) # influential outliers

# refit model to stabilize variance
# using a log (x+1) transformation
log_bm_minus_d0 <- mutate(inv_bm_minus_d0, 
  log_tot_inv_bm = log(tot_inv_bm_mg + 1)
  ) 

lm_log_inv_bm_minus_d0 <- lm(
  log_tot_inv_bm ~ Richness*Density, 
  data = log_bm_minus_d0 
)

# double-check assumptions for two-way ANOVA
plot(lm_log_inv_bm_minus_d0, 2) # normality of residuals
plot(lm_log_inv_bm_minus_d0, 1) # homogeneity of variance
plot(lm_log_inv_bm_minus_d0, 4) # influential outliers

## two-way ANOVA ----

# run a type II ANOVA to account for unbalanced sample designs
aov_t2_inv_bm <- car::Anova(lm_log_inv_bm_minus_d0, type = 2)

## pairwise comparisons ----

# calculate estimated marginal means
# to account for unbalanced sample design
ref_grid_inv_bm <- emmeans::ref_grid(lm_log_inv_bm_minus_d0)
ref_grid_inv_bm@grid

# visualize potential interactions before doing any statistical comparisons
emmeans::emmip(lm_log_inv_bm_minus_d0, Richness ~ Density)

# run pairwise comparisons
EMM <- emmeans::emmeans(lm_log_inv_bm_minus_d0, ~Richness*Density)
pairs_dense <- test(pairs(EMM, by = "Density"))
pairs_rich<- test(pairs(EMM, by = "Richness"))

# clean datasets for placing into the supp info
pairs_dense_df <- pairs_dense %>%
  as.data.frame() %>%
  rename(treatment = Density)

pairs_rich_df <- pairs_rich %>%
  as.data.frame() %>%
  rename(treatment = Richness)

table_s2 <- pairs_dense_df %>%
  rbind(pairs_rich_df) %>%
  mutate(
    estimate = round(estimate, digits = 2),
    SE = round(SE, digits = 2), 
    t.ratio = round(t.ratio, digits = 2),
    p.value = round(p.value, digits = 2)
    )

# save to disk ----

write.csv(
  table_s2, 
  file = here("output", "data_appendix_output", "table_s2_pairwise_inv_bm.csv"),
  row.names = FALSE
)

