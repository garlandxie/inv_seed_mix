# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(ggplot2)    # for visualizing plots

# import -----------------------------------------------------------------------
bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))

# aboveground biomass ----------------------------------------------------------

## invasive species ------------------------------------------------------------
bm_abg_inv <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB" &
      Species_ID == "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(
    inv_abg_bm_g = sum(Biomass_g, na.rm = TRUE),
    num_invaders = n()
  ) %>%
  ungroup() %>%
  arrange(Richness_ID, Density_ID, Rep) 

# replace a zero with a missing value 
# M4-D2-5 had a missing aboveground sample, so there should be NA's for 
# average invader biomass and total invader biomass

# also: there is a root sample for one invader (C1) so number of invaders 
# for the M4-D2-5 tray can still be 1 

bm_abg_inv$inv_abg_bm_g[which(
  bm_abg_inv$Richness_ID == "M4" &
    bm_abg_inv$Density_ID == "D2" &
    bm_abg_inv$Rep == "5")] <- NA

## native species --------------------------------------------------------------
bm_abg_res <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB" &
      Species_ID != "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g) 
  ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(
    res_abg_bm_g = sum(Biomass_g, na.rm = TRUE),
    realized_sr = length(unique(Species_ID))) %>%
  ungroup()

# below-ground biomass ---------------------------------------------------------

## invasive species ------------------------------------------------------------
bm_roots_inv <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "BM" &
      Species_ID == "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(inv_roots_bm_g = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(Richness_ID, Density_ID, Rep) 

## native species --------------------------------------------------------------
bm_roots_res <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "BM" &
      Species_ID != "CIAR"
  ) %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g) 
  ) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(
    res_root_bm_g = sum(Biomass_g, na.rm = TRUE)
    ) %>%
  ungroup()

# data joins -------------------------------------------------------------------

roots_tidy <- dplyr::inner_join(
  bm_roots_res,
  bm_roots_inv, 
  by = c("Richness_ID", "Density_ID", "Rep")
  )

abg_tidy <- dplyr::inner_join(
  bm_abg_res,
  bm_abg_inv, 
  by = c("Richness_ID", "Density_ID", "Rep")
)

bm_tidy <- dplyr::inner_join(
  roots_tidy, 
  abg_tidy, 
  by = c("Richness_ID", "Density_ID", "Rep")
)

# save to disk -----------------------------------------------------------------

write.csv(
  x = bm_tidy, 
  file = here("data", "intermediate_data", "biomass.csv")
)
