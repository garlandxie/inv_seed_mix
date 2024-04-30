# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(ggplot2)    # for visualizing plots
library(tidyr)      # for handling missing values
library(tidyr)      # for replacing missing values
library(lubridate)  # for dealing with dates
library(stringr)    # for dealing with string characters

# import -----------------------------------------------------------------------

bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))
height_raw <- read.csv(here("data", "input_data", "height.csv"))
leaf_area_raw <- read.csv(here("data", "input_data", "leaf_area.csv"))
leaf_bm_raw <- read.csv(here("data", "input_data", "leaf_biomass.csv"))

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
# get the mean height of early, mid, and late sampling periods 
# as a separate data frame; can use inner join to merge into a single dataset
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
  inner_join(
    height_summ_mid, 
    by = c("Richness_ID", "Density_ID", "Rep", "Spp")
    ) %>%
  inner_join(
    height_summ_late, 
    by = c("Richness_ID", "Density_ID", "Rep", "Spp")
    )

# calculate time difference (t2 - t1) between late and early sampling
# for the relative growth rate (log(r2) - log(r1))/t2 - t1
time_diff_earl_late <- as.numeric(ymd("2023-06-22") - ymd("2023-04-18"))
  
# calculate relative growth rate of height (rgr) as a column variable 
height_final <- height_summ %>%
  mutate(
    log_early_height = log(early_height + 1),
    log_late_height = log(late_height), 
    log_diff = log_late_height - log_early_height, 
    rgr_height = log_diff/time_diff_earl_late,
    rgr_height = round(rgr_height, digits = 3)
  ) 

# clean: specific leaf area ----------------------------------------------------

leaf_area_tidy <- leaf_area_raw %>%
  mutate(
    Richness_ID = stringr::str_split_i(tray_id, pattern = "-", 1),
    Density_ID = stringr::str_split_i(tray_id, pattern = "-", 2),
    Rep = stringr::str_split_i(tray_id, pattern = "-", 3))

sla <- leaf_area_tidy %>%
  dplyr::inner_join(
    leaf_bm_raw, by = c(
      "Richness_ID", "Density_ID", "Rep", "species", "ind", "leaf_rep")
    )


