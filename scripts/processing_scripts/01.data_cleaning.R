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



