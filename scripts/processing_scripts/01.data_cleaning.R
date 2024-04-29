# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(ggplot2)    # for visualizing plots
library(tidyr)      # for handling missing values
library(tidyr)      # for replacing missing values
library(lubridate)  # for dealing with dates

# import -----------------------------------------------------------------------

bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))
height_raw <- read.csv(here("data", "input_data", "height.csv"))

# clean: height ----------------------------------------------------------------

height_tidy <- height_raw %>%
  mutate(
    
    # plant height should be a numeric value
    # ensure that NA's are explicit in this dataset 
    # for calculating summary statistics
    Height_mm = dplyr::na_if(Height_mm, "-"),
    Height_mm = as.numeric(Height_mm),
  
    # group into three distinct observations
    # hopefully this makes it a bit easier to code downstream
    Session = case_when(
      Date_ymd == "2023-04-17" ~ "early", 
      Date_ymd == "2023-04-18" ~ "early",
      Date_ymd == "2023-05-16" ~ "mid", 
      Date_ymd == "2023-05-17" ~ "mid", 
      Date_ymd == "2023-05-18" ~ "mid", 
      Date_ymd == "2023-06-22" ~ "late", 
      Date_ymd == "2023-06-23" ~ "late", 
      TRUE ~ Date_ymd
    )
  )

# calculate summary statistics
height_summ1 <- height_tidy %>%
  group_by(Session, Tray_ID, Richness_ID, Density_ID, Rep, Spp) %>%
  summarize(
    
    # calculate average plant height (and its associated standard deviation)
    mean_height = mean(Height_mm, na.rm = TRUE) %>% round(digits = 2),
    sd_height = sd(Height_mm, na.rm = TRUE) %>% round(digits = 2),
    
    # replace missing value with zeros for standard deviations
    sd_height = tidyr::replace_na(sd_height, 0),
    
    # calculate sample size 
    # number of sampled individual plants per species per tray
    sample_size = n()
    
  ) %>%
  ungroup()
  
# this is not an elegant solution, but it works!
height_summ_earl <- height_summ1 %>%
  dplyr::filter(Session == "early") %>%
  tidyr::pivot_wider(names_from = Session, values_from = mean_height) %>%
  dplyr::select(Richness_ID, Density_ID, Rep, Spp, early_height = early)

height_summ_mid <- height_summ1 %>%
  dplyr::filter(Session == "mid") %>%
  tidyr::pivot_wider(names_from = Session, values_from = mean_height) %>%
  dplyr::select(Richness_ID, Density_ID, Rep, Spp, mid_height = mid)

height_summ_late <- height_summ1 %>%
  dplyr::filter(Session == "late") %>%
  tidyr::pivot_wider(names_from = Session, values_from = mean_height) %>%
  dplyr::select(Richness_ID, Density_ID, Rep, Spp, late_height = late)

height_summ <- height_summ_earl %>%
  inner_join(height_summ_mid, by = c("Richness_ID", "Density_ID", "Rep", "Spp")) %>%
  inner_join(height_summ_late, by = c("Richness_ID", "Density_ID", "Rep", "Spp"))

time_diff_earl_late <- as.numeric(ymd("2023-06-22") - ymd("2023-04-18"))
  
height_final <- height_summ %>%
  mutate(
    log_early_height = log(early_height + 1),
    log_late_height = log(late_height), 
    log_diff = log_late_height - log_early_height, 
    rgr_height = log_diff/time_diff_earl_late,
    rgr_height =round(rgr_height, digits = 3)
  ) 

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

# clean: invasive species ------------------------------------------------------

# total biomass for Cirsium arvense
biomass_ab_invasive <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB" &
    Species_ID == "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g), 
    Biomass_mg = Biomass_g * 1000
    ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(
    mean_inv_bm = mean(Biomass_mg, na.rm = TRUE) %>% round(digits = 2), 
    tot_inv_bm = sum(Biomass_mg, na.rm = TRUE),
    num_invaders = n()
    ) %>%
  ungroup() %>%
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

