# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(FD)         # for calculating community-weighted means values
library(tidyr)      # for handling missing values

# import -----------------------------------------------------------------------

bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))
leaf_area_raw <- read.csv(here("data", "input_data", "leaf_area.csv"))
leaf_bm_raw <- read.csv(here("data", "input_data", "leaf_biomass.csv"))

# data clean -------------------------------------------------------------------

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

# get SLA for each species
# across three plant density treatments
sla_summary1 <- sla %>%
  group_by(Richness_ID, Density_ID, Rep, species, ind) %>%
  summarize(mean_sla = mean(sla, na.rm = TRUE)) %>%
  ungroup() %>% 
  group_by(Density_ID, species) %>%
  summarize(mean_sla = mean(mean_sla, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mean_sla = round(mean_sla, digits = 2)) %>%
  arrange(species)

# from a previous analysis, SLA for each species 
# across three plant density treatments were very similar
# so, just group the density treatments to get an average SLA value
# for each species 
sla_summary2 <- sla %>%
  filter(species != "CIAR") %>%
  group_by(Richness_ID, Density_ID, Rep, species, ind) %>%
  summarize(mean_sla = mean(sla, na.rm = TRUE)) %>%
  ungroup() %>% 
  group_by(species) %>%
  summarize(mean_sla = mean(mean_sla, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(mean_sla = round(mean_sla, digits = 2)) %>%
  tibble::column_to_rownames(var = "species")

sla_ciar <- sla %>%
  filter(species == "CIAR") %>%
  group_by(species) %>%
  summarize(sla = mean(sla, na.rm = TRUE)) 

## relative abundance ----------------------------------------------------------

# calculate relative abundance of each resident species for 
# each community 
biomass_ab_res_rel <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB" &
      Species_ID != "CIAR"
  ) %>%
  mutate(
    spp_biomass_g = na_if(Biomass_g, "-"),
    spp_biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Richness_ID, Density_ID, Rep, Species_ID) %>%
  summarize(
    spp_biomass_g = sum(spp_biomass_g, na.rm = TRUE),
    spp_biomass_g = round(spp_biomass_g, digits = 2)) %>%
  
  # fix a typo
  mutate(Species_ID = case_when(
    Species_ID == "OEBI " ~ "OEBI", 
    TRUE ~ Species_ID)
  ) %>%
  ungroup() %>%
  
  # create a community data matrix 
  # with abundance as cells and tray ID's as rownames
  group_by(Richness_ID, Density_ID, Rep) %>%
  tidyr::pivot_wider(names_from = Species_ID, values_from = spp_biomass_g) %>%
  mutate(
    across(ANGE:RUHI, ~replace_na(., 0)),
    tray = paste(Richness_ID, Density_ID, Rep, sep = "-")
  ) %>%
  ungroup() %>%
  dplyr::select(-c(Richness_ID, Density_ID, Rep)) %>%
  tibble::column_to_rownames(var = "tray") %>%
  dplyr::select(ANGE, HEHE, MOFI, OEBI, RUHI)

# remove zero sum abundances
biomass_ab_res_rel2 <- biomass_ab_res_rel[rowSums(biomass_ab_res_rel)>0,]

## weighted community dissimilarity --------------------------------------------

# calculate community-weighted mean SLA values
# for the resident community
cwm_sla_res <- FD::dbFD(
  x = sla_summary2, a = as.matrix(biomass_ab_res_rel2)
)

# calculate weighted community dissimilarity 
cwm_res <- cwm_sla_res$CWM %>%
  tibble::rownames_to_column(var = "tray_id") %>%
  mutate(
    Richness_ID = sapply(strsplit(tray_id, split = "-"), "[[", 1),
    Density_ID = sapply(strsplit(tray_id, split = "-"), "[[", 2),
    Rep = sapply(strsplit(tray_id, split = "-"), "[[", 3)
  ) %>%
  rename(cwm_res_sla = mean_sla) %>%
  mutate(wds_sla = abs(cwm_res_sla - sla_ciar$sla)) %>%
  dplyr::select(Richness_ID, Density_ID, Rep, cwm_res_sla, wds_sla)

# save to disk -----------------------------------------------------------------

write.csv(
  cwm_res, 
  file = here("data", "intermediate_data", "cwm_res.csv"),
  row.names = FALSE
)