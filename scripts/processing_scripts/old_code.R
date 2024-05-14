# dump some code for exploratory data analysis (might use it later)

height_final %>%
  dplyr::filter(Density_ID %in% c("D1", "D2", "D3")) %>%
  ggplot(aes(x = Richness_ID, y = rgr_height, fill = Spp)) + 
  geom_boxplot() + 
  facet_wrap(~Density_ID) + 
  labs(
    x = "Seeding Richness", 
    y = "Relative Growth Rate (Height)") + 
  theme_bw()

# clean: resident community ----------------------------------------------------

## resident community biomass per tray -----------------------------------------
biomass_ab_native <- bm_raw %>%
  dplyr::filter(Biomass_Type == "AB", Species_ID != "CIAR") %>%
  mutate(Biomass_g = as.numeric(Biomass_g)) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(tot_native_bm_g = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(status = "N")

# species identity - monocultures ---------------------------------------------
biomass_mono_native <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB", 
    Species_ID != "CIAR", 
    Richness_ID == "M1") %>%
  mutate(Biomass_g = as.numeric(Biomass_g)) %>%
  arrange(Richness_ID, Density_ID, Rep) 

# plot: resident community -----------------------------------------------------

## resident community biomass --------------------------------------------------
(ab_biomass_dens_rich <- biomass_ab_native %>%
   
   mutate(
     Richness_ID = case_when(
       Richness_ID == "M1" ~ "1", 
       Richness_ID == "M2" ~ "2", 
       Richness_ID == "M4" ~ "4", 
       TRUE ~ Richness_ID) %>% factor(levels = c("1", "2", "4")),
     
     Density_ID = case_when(
       Density_ID == "D1" ~ "100 seeds", 
       Density_ID == "D2" ~ "200 seeds", 
       Density_ID == "D3" ~ "400 seeds", 
       TRUE ~ Density_ID)
     
   ) %>%
   ggplot(aes(x = Richness_ID, y = tot_native_bm_g)) +
   geom_boxplot(aes(fill = Density_ID), position = "dodge2") +
   labs(
     fill = "Sown Seed Density",
     x = "Sown Species Richness", 
     y = "Resident Aboveground Community Biomass (g)") + 
   theme_bw()
)

## species identity - monocultures -------------------------------------
ab_biomass_spp_mono <- biomass_mono_native %>%
  mutate(
    
    Species_ID = factor(
      Species_ID, 
      levels = c("HEHE", "RUHI", "ANGE", "MOFI", "OEBI")
    ), 
    
    Density_ID = case_when(
      Density_ID == "D1" ~ "48 seeds", 
      Density_ID == "D2" ~ "200 seeds", 
      Density_ID == "D3" ~ "400 seeds", 
      TRUE ~ Density_ID
    )
  )%>%
  ggplot(aes(x = Species_ID, y = Biomass_g)) + 
  geom_boxplot(aes(fill = Density_ID), position = "dodge2") + 
  labs(
    fill = "Sown Seed Density",
    x = "Native Species within Seed Mix",
    y = "Aboveground Community Biomass (g)"
  ) + 
  theme_bw()


# data exploration: invasive species -------------------------------------------

## relationship between mean invader and community-level invader biomass -------
biomass_ab_invasive %>%
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
  ggplot(aes(
    x = mean_inv_bm, 
    y = tot_inv_bm, 
    fill = Density_ID, 
    col = Density_ID)
  ) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~Richness_ID) + 
  labs(
    x = "Mean invader Biomass (mg)",
    y = "Community-level invader Biomass (mg)"
  ) + 
  theme_bw()

## invader survival -----------------------------------------------------------
biomass_ab_invasive %>%
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
      ),
    
    Richness_ID = case_when(
      Density_ID  == "20 invasive seeds"  ~ "Invasive-Only", 
      Richness_ID == "M1" ~ "Native Monoculture", 
      Richness_ID == "M2" ~ "2-species Mixture", 
      Richness_ID == "M4" ~ "4-species Mixture", 
      TRUE ~ Richness_ID
    ) %>%
      factor(levels = c(
        "Native Monoculture", "2-species Mixture", "4-species Mixture")
      )
  ) %>%
  ggplot(aes(x = Richness_ID, y = num_invaders)) + 
  geom_boxplot(aes(fill = Density_ID), position = "dodge2") +
  ylim(0,15) + 
  labs(y = "Invader survival") + 
  theme_bw()

## community-level invader biomass ---------------------------------------------
(plot_tot_inv_bm <- biomass_ab_invasive %>%
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
   )%>%
   
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
   ggplot(aes(x = Richness_ID, y = tot_inv_bm)) +
   geom_boxplot(
     aes(fill = Density_ID), 
     position = "dodge2"
   ) + 
   labs(
     fill = "Sown Seed Density", 
     x = NULL, 
     y = "Community-level biomass (g) of C. arvense"
   ) + 
   theme_bw()
)

## average invader biomass -----------------------------------------------------
(plot_avg_inv_bm <- biomass_ab_invasive %>%
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
   )%>%
   
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
   ggplot(aes(x = Richness_ID, y = tot_inv_bm)) +
   geom_boxplot(
     aes(fill = Density_ID), 
     position = "dodge2"
   ) + 
   labs(
     fill = "Sown Seed Density", 
     x = NULL, 
     y = "Average biomass (g) of C. arvense"
   ) + 
   theme_bw()
)

# save to disk -----------------------------------------------------------------

ggsave(
  plot = ab_biomass_dens_rich, 
  file = here("output", "results", "native_bm_dens_rich.png"), 
  device = "png", 
  height = 4, 
  width = 7
)

ggsave(
  plot = ab_biomass_spp_mono, 
  file = here("output", "results", "ab_biomass_spp_mono.png"), 
  device = "png", 
  height = 4, 
  width = 7
)

ggsave(
  plot = biomass_ab_invasive_plot, 
  file = here("output", "results", "biomass_ab_invasive.png"), 
  device = "png", 
  height = 4, 
  width = 7
)

# libraries -----
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data 
library(ggplot2)    # for visualizing data  

# import ----
c_germ <- read.csv(here("data", "input_data", "cumulative_germination.csv"))
bm_res <- read.csv(here("data", "intermediate_data", "resident_biomass.csv"))

# data cleaning ----

## cumulative germination ----
c_germ_tidy <- c_germ %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ, names_from = spp) %>%
  ungroup() %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep, 
    
    cum_germ_ciar = CIAR, 
    cum_germ_oebi = OEBI, 
    cum_germ_ruhi = RUHI, 
    cum_germ_ange = ANGE, 
    cum_germ_hehe = HEHE, 
    cum_germ_mofi = MOFI) %>%
  mutate(across(cum_germ_ciar:cum_germ_mofi, ~replace(., is.na(.), 0))) %>%
  mutate(
    cum_germ_res = 
      cum_germ_oebi + 
      cum_germ_ruhi + 
      cum_germ_ange + 
      cum_germ_hehe + 
      cum_germ_mofi
  ) 

## cumulative germination (as a percentage) ----

# invader-only
c_germ_perc_ciar <- c_germ %>%
  filter(spp == "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ_perc, names_from = spp) %>%
  ungroup() %>%
  rename(cum_germ_perc_ciar = CIAR)

c_germ_perc_res <- c_germ %>% 
  filter(spp != "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  summarize(cum_germ_res = sum(cum_germ, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    sown_seeds = case_when(
      density_id == "D1" ~ 48,
      density_id == "D2" ~ 200, 
      density_id == "D3" ~ 400
    ), 
    cum_germ_perc_res = cum_germ_res/sown_seeds,
    cum_germ_perc_res = round(cum_germ_perc_res, digits = 2)
  ) 

# percentage for all native species within a tray 
c_germ_perc <- c_germ_perc_ciar %>% 
  inner_join(
    c_germ_perc_res, 
    by = c("week", "richness_id", "density_id", "rep")
  ) %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep,
    sown_seeds = sown_seeds.y, 
    cum_germ_perc_ciar, 
    cum_germ_perc_res
  )

# for data visualizing
c_germ_perc_tidy <- c_germ %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ, names_from = spp) %>%
  ungroup() %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep, 
    
    cum_germ_perc_ciar = CIAR, 
    cum_germ_perc_oebi = OEBI, 
    cum_germ_perc_ruhi = RUHI, 
    cum_germ_perc_ange = ANGE, 
    cum_germ_perc_hehe = HEHE, 
    cum_germ_perc_mofi = MOFI) %>%
  
  mutate(across(
    cum_germ_perc_ciar:cum_germ_perc_mofi, 
    ~replace(., is.na(.), 0))
  ) %>%
  
  mutate(
    cum_germ_perc_res = 
      cum_germ_perc_oebi + 
      cum_germ_perc_ruhi + 
      cum_germ_perc_ange + 
      cum_germ_perc_hehe + 
      cum_germ_perc_mofi
  ) 

# data visualization ----

## cumulative germination ----
c_germ_tidy %>%
  group_by(week, richness_id, density_id) %>%
  summarize(
    mean_cum_germ_ciar = mean(cum_germ_ciar),
    mean_cum_germ_oebi = mean(cum_germ_oebi),
    mean_cum_germ_ruhi = mean(cum_germ_ruhi), 
    mean_cum_germ_ange = mean(cum_germ_ange), 
    mean_cum_germ_hehe = mean(cum_germ_hehe), 
    mean_cum_germ_mofi = mean(cum_germ_mofi)
  ) %>%
  ungroup() %>%
  dplyr::filter(
    richness_id %in% c("M1", "M2", "M4"), 
    density_id %in% c("D1", "D2", "D3")
  ) %>%
  ggplot() +
  geom_line(aes(x = week, y = mean_cum_germ_ciar), col = "purple") +
  geom_line(aes(x = week, y = mean_cum_germ_oebi), col = "orange") + 
  geom_line(aes(x = week, y = mean_cum_germ_ruhi), col = "red") + 
  geom_line(aes(x = week, y = mean_cum_germ_ange), col = "blue") + 
  geom_line(aes(x = week, y = mean_cum_germ_mofi), col = "green") + 
  geom_line(aes(x = week, y = mean_cum_germ_hehe), col = "grey") +
  labs(
    x = "Week of Observation", 
    y = "Cumulative germination") + 
  facet_wrap(richness_id~density_id) + 
  theme_bw()

## cumulative germination (as a percentage) -----

c_germ_perc_tidy %>%
  group_by(week, density_id, richness_id) %>%
  summarize(
    mean_cum_perc_germ_ciar = mean(cum_germ_perc_ciar),
    mean_cum_perc_germ_oebi = mean(cum_germ_perc_oebi),
    mean_cum_perc_germ_ruhi = mean(cum_germ_perc_ruhi), 
    mean_cum_perc_germ_ange = mean(cum_germ_perc_ange), 
    mean_cum_perc_germ_hehe = mean(cum_germ_perc_hehe), 
    mean_cum_perc_germ_mofi = mean(cum_germ_perc_mofi)
  ) %>%
  ungroup() %>%
  dplyr::filter(
    density_id %in% c("D1", "D2", "D3") & 
      richness_id %in% c("M1", "M2", "M4")
  ) %>%
  ggplot() +
  geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
  geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "grey") + 
  labs(x = "Week of Observation", y = "Cumulative germination") + 
  facet_wrap(density_id~richness_id) + 
  theme_bw() 

# some analyses ----
c_germ_w12 <- dplyr::filter(c_germ_perc, week == 12) 
glm_perc_ciar <- glm(
  cum_germ_perc_ciar ~ richness_id*density_id, 
  family = binomial(link = "logit"), 
  data = c_germ_w12
)

sum(residuals(glm_perc_ciar, "pearson")^2) / glm_perc_ciar$df.residual

glm_perc_ciar2 <- glm(
  cum_germ_perc_ciar ~ richness_id*density_id, 
  family = quasibinomial(link = "logit"), 
  data = c_germ_w12
)

car::Anova(glm_perc_ciar2, type = "II")

# resident community biomass
lm_bm_res <- lm(res_comm_biomass_mg ~ Richness_ID*Density_ID, data = bm_res)
car::Anova(lm_bm_res, type = "II")  
emmeans::emmeans(lm_bm_res, ~Richness_ID) %>% pairs()

height_final <- janitor::clean_names(height_final)

df_clean <- bm_res %>%
  janitor::clean_names() %>%
  inner_join(c_germ_w12, by = c("richness_id", "density_id", "rep")) %>%
  inner_join(height_final, by = c("richness_id", "density_id", "rep"))

df_clean %>%
  ggplot(aes(x = res_comm_biomass_mg, y= cum_germ_perc_ciar)) + 
  geom_point() + 
  geom_smooth(method = "lm")

glm_res_ciar <- glm(cum_germ_perc_ciar ~ res_comm_biomass_mg, 
                    family = quasibinomial(link = "logit"), 
                    data = df_clean
)

summary(glm_res_ciar)

df_clean %>%
  ggplot(aes(x = cum_germ_perc_res, y = res_comm_biomass_mg)) + 
  geom_point(aes(col = rep)) + 
  geom_smooth(method = "lm")

res_germ_bm <- lm(res_comm_biomass_mg ~ cum_germ_perc_res, data = df_clean)
summary(res_germ_bm)
plot(res_germ_bm)

df_clean %>%
  ggplot(aes(x = cum_germ_perc_res, y = rgr_height)) + 
  geom_point() + 
  geom_smooth(method = "lm")

lm_rgr_perc_res <- lm(rgr_height ~ cum_germ_perc_res, data = df_clean) 
summary(lm_rgr_perc_res)
plot(lm_rgr_perc_res)

df_clean %>%
  ggplot(aes(x = rgr_height, y = res_comm_biomass_mg)) + 
  geom_point() + 
  geom_smooth(method = "lm") 

lm_rgr_bm_res <- lm(res_comm_biomass_mg ~ rgr_height, data = df_clean)
summary(lm_rgr_bm_res)


# save to disk -----

write.csv(
  c_germ_perc, 
  file = here("data", "intermediate_data", "c_germ_clean.csv")
)


# libraries -----
library(dplyr)    # for manipulating data
library(here)     # for creating relative file paths
library(ggplot2)  # for visualizing data 

# import ----

inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))

# exploratory data analysis -----

## test correlation between total invader biomass and mean invader biomass
# total invader biomass: sum of biomass from all invaders within a tray
# mean invader biomass: average biomass per tray across all individual invaders
inv_bm %>%
  ggplot(aes(
    x = mean_inv_bm_mg, 
    y = tot_inv_bm_mg, 
    fill = Density_ID, 
    col = Density_ID)
  ) + 
  geom_point(alpha = 0.5) +
  geom_smooth(method = "lm", se = FALSE) + 
  facet_wrap(~Richness_ID) + 
  labs(
    x = "Mean invader Biomass (mg)",
    y = "Community-level invader Biomass (mg)"
  ) + 
  theme_bw()

bm_cor <- inv_bm[complete.cases(inv_bm),]
cor(x = bm_cor$mean_inv_bm_mg,y = bm_cor$tot_inv_bm_mg)

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

germ_perc_ciar <- read.csv(here("data", "intermediate_data", "c_germ_perc_ciar.csv"))
inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))


germ_perc_inv <- germ_perc %>%
  dplyr::filter(week == 12) %>%
  dplyr::select(richness_id, density_id, rep, cum_germ_perc_ciar)

germ_perc_res <- germ_perc %>%
  dplyr::filter(week == 12) %>%
  dplyr::select(richness_id, density_id, rep, cum_germ_perc_res)

df_clean <- inv_bm %>%
  janitor::clean_names() %>% 
  inner_join(germ_perc_inv, by = c("richness_id", "density_id", "rep")) 

df_clean %>%
  filter(
    density_id %in% c("D1", "D2", "D3") & richness_id %in% c("M1", "M2", "M4")) %>%
  ggplot(aes(x = cum_germ_perc_ciar, y = tot_inv_bm_mg)) + 
  geom_point(aes(col = rep)) + 
  geom_smooth(method = "lm") 

plot(lm(log(tot_inv_bm_mg + 0.1) ~ cum_germ_perc_ciar, data = df_clean))
