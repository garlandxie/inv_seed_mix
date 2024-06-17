# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(ggplot2)    # for visualizing plots
library(tidyr)      # for handling missing values
library(tidyr)      # for replacing missing values
library(lubridate)  # for dealing with dates
library(stringr)    # for dealing with string characters
library(GerminaR)   # for summarizing germination indices

# import -----------------------------------------------------------------------

bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))
height_raw <- read.csv(here("data", "input_data", "height.csv"))
leaf_area_raw <- read.csv(here("data", "input_data", "leaf_area.csv"))
leaf_bm_raw <- read.csv(here("data", "input_data", "leaf_biomass.csv"))
c_germ_raw <- read.csv(here("data", "input_data", "cumulative_germination.csv"))

# data cleaning ----------------------------------------------------------------

## invasive biomass ------------------------------------------------------------

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
    mean_inv_bm_mg = mean(Biomass_mg, na.rm = TRUE) %>% round(digits = 2), 
    tot_inv_bm_mg = sum(Biomass_mg, na.rm = TRUE),
    num_invaders = n()
  ) %>%
  ungroup() %>%
  arrange(Richness_ID, Density_ID, Rep) 

# replace a zero with a missing value 
# M4-D2-5 had a missing aboveground sample, so there should be NA's for 
# average invader biomass and total invader biomass

# also: there is a root sample for one invader (C1) so number of invaders 
# for the M4-D2-5 tray can still be 1 

biomass_ab_invasive$mean_inv_bm_mg[which(
  biomass_ab_invasive$Richness_ID == "M4" &
  biomass_ab_invasive$Density_ID == "D2" &
  biomass_ab_invasive$Rep == "5")] <- NA

biomass_ab_invasive$tot_inv_bm_mg[which(
  biomass_ab_invasive$Richness_ID == "M4" &
  biomass_ab_invasive$Density_ID == "D2" &
  biomass_ab_invasive$Rep == "5")] <- NA

## height ----------------------------------------------------------------------

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


## specific leaf area ----------------------------------------------------------

leaf_area_tidy <- leaf_area_raw %>%
  mutate(
    Richness_ID = stringr::str_split_i(tray_id, pattern = "-", 1),
    Density_ID = stringr::str_split_i(tray_id, pattern = "-", 2),
    Rep = stringr::str_split_i(tray_id, pattern = "-", 3)
    )

sla <- leaf_area_tidy %>%
  
  # perform inner join to get the leaf biomass and area per individual specimen
  dplyr::inner_join(
    leaf_bm_raw, 
    by = c("Richness_ID", "Density_ID", "Rep", "species", "ind", "leaf_rep")
    ) %>%
  
  # select appropriate columns required to do the SLA calculations
    dplyr::select(
      Richness_ID, 
      Density_ID, 
      Rep, 
      species, 
      ind, 
      leaf_rep,
      leaf_area_cm2, 
      biomass_mg
      ) %>%
  
    # calculate specific leaf area (leaf_area_cm2/biomass_mg)
    # per individual specimen
    mutate(sla = leaf_area_cm2/biomass_mg, 
           sla = round(sla, digits = 2)
           ) 

## resident community biomass --------------------------------------------------

biomass_ab_res <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB" &
      Species_ID != "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g), 
    Biomass_mg = Biomass_g * 1000
  ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(res_comm_biomass_mg = sum(Biomass_mg, na.rm = TRUE)) %>%
  ungroup()

## cumulative percentage germination -------------------------------------------

# invader-only
c_germ_perc_ciar <- c_germ_raw %>%
  filter(spp == "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ_perc, names_from = spp) %>%
  ungroup() %>%
  rename(cum_germ_perc_ciar = CIAR)

c_germ_perc_res <- c_germ_raw %>% 
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

# for now, obtain the cumulative percentage germination 
# at the end of the experiment 
c_germ_w12 <- dplyr::filter(c_germ_perc, week == 12) 

## germination speed ------------------------------------------------------------

# calculate the number of seeds germinated per week for each samples species
c_gs_w1 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 1 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W1" = "1")

c_gs_w2 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 2 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W2" = "2")

c_gs_w3 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 3 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W3" = "3")

c_gs_w4 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 4 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ)  %>%
  ungroup() %>%
  rename("W4" = "4")

c_gs_w5 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 5 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W5" = "5")

c_gs_w6 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 6 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ)  %>%
  ungroup() %>%
  rename("W6" = "6")

c_gs_w7 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 7 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W7" = "7")

c_gs_w8 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 8 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W8" = "8")
  
c_gs_w9 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 9 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W9" = "9")

c_gs_w10 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 10 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W10" = "10")

c_gs_w11 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 11 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W11" = "11")

c_gs_w12 <- c_germ_raw %>%
  dplyr::select(week, richness_id, density_id, rep, spp, sown_seeds, cum_germ) %>%
  filter(week == 12 & richness_id != "-") %>%
  group_by(richness_id, density_id, rep, spp) %>%
  tidyr::pivot_wider(names_from = week, values_from = cum_germ) %>%
  ungroup() %>%
  rename("W12" = "12")

multi_key_id <- c("richness_id", "density_id", "rep", "spp", "sown_seeds")

c_gs <- c_gs_w1 %>%
  inner_join(c_gs_w2, by = multi_key_id) %>%
  inner_join(c_gs_w3, by = multi_key_id) %>%
  inner_join(c_gs_w4, by = multi_key_id) %>%
  inner_join(c_gs_w5, by = multi_key_id) %>%
  inner_join(c_gs_w6, by = multi_key_id) %>%
  inner_join(c_gs_w7, by = multi_key_id) %>%
  inner_join(c_gs_w8, by = multi_key_id) %>%
  inner_join(c_gs_w9, by = multi_key_id) %>%
  inner_join(c_gs_w10, by = multi_key_id) %>%
  inner_join(c_gs_w11, by = multi_key_id) %>%
  inner_join(c_gs_w12, by = multi_key_id) 

c_gs2 <- c_gs %>%
  mutate(
    W2_diff = W2 - W1, 
    W3_diff = W3 - W2, 
    W4_diff = W4 - W3, 
    W5_diff = W5 - W4, 
    W6_diff = W6 - W5, 
    W7_diff = W7 - W6, 
    W8_diff = W8 - W7, 
    W9_diff = W9 -W8, 
    W10_diff = W10 - W9, 
    W11_diff = W11 - W10, 
    W12_diff = W12 - W11
    ) %>%
  select(
    richness_id, 
    density_id, 
    rep, 
    spp, 
    sown_seeds, 
    W1, 
    W2_diff, 
    W3_diff, 
    W4_diff, 
    W5_diff, 
    W6_diff, 
    W7_diff, 
    W8_diff, 
    W9_diff, 
    W10_diff, 
    W11_diff, 
    W12_diff
  ) 

# calculate germination speed per sampled species
c_gs2 <- c_gs2 %>%
    mutate(
    mgt = ger_MGT(evalName = "W", data = .),
    gsp = ger_GSP(evalName = "W", data = .)   
    ) 

# calculate germination speed and mean germination time per tray
# for the resident community 
c_gs_res <- c_gs2 %>%
  dplyr::filter(!(spp %in% c("CIAR"))) %>%
  group_by(richness_id, density_id, rep) %>%
  summarize(
    mean_gsp = mean(gsp, na.rm = TRUE),
    mean_mgt = mean(mgt, na.rm = TRUE)
  ) %>%
  ungroup() 

# create a dataset with germination speed, mean germination time, 
# and percent germination for each tray
c_germ <- left_join(
  c_germ_w12, 
  c_gs_res, 
  by = c("richness_id", "density_id", "rep")
)

# some exploratory analyses ----------------------------------------------------

# may move this to another R script
c_germ %>%
  ggplot(aes(x = density_id, y = mean_gsp)) + 
  geom_boxplot()

summary(lm(mean_gsp ~ density_id*richness_id, 
            data = c_germ))

c_germ %>%
  ggplot(aes(x = mean_gsp, y = cum_germ_perc_res)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE)

summary(glm(cum_germ_perc_res ~ mean_gsp, 
            weights = sown_seeds,
            family = binomial(link = "logit"), 
            data = c_germ))

# save to disk -----------------------------------------------------------------

## invader biomass -----------------------s--------------------------------------

write.csv(
  biomass_ab_invasive, 
  file = here("data", "intermediate_data", "invader_biomass.csv"),
  row.names = FALSE
)

## height ----------------------------------------------------------------------

write.csv(
  height_final, 
  file = here("data", "intermediate_data", "rgr_height.csv"),
  row.names = FALSE
)

## specific leaf area ----------------------------------------------------------

write.csv(
  sla,
  file = here("data", "intermediate_data", "specific_leaf_area.csv"),
  row.names = FALSE
)

## resident community biomass --------------------------------------------------

write.csv(
  biomass_ab_res, 
  file = here("data", "intermediate_data", "resident_biomass.csv"),
  row.names = FALSE
)

## cumulative percentage germination -------------------------------------------

write.csv(
  c_germ_w12, 
  file = here("data", "intermediate_data", "c_germ_week.csv")
)
