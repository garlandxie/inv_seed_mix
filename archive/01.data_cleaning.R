# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(ggplot2)    # for visualizing plots
library(tidyr)      # for handling missing values
library(tidyr)      # for replacing missing values
library(lubridate)  # for dealing with dates
library(stringr)    # for dealing with string characters
library(GerminaR)   # for summarizing germination indices
library(FD)         # for calculating community-weighted means values

# import -----------------------------------------------------------------------

bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))
height_raw <- read.csv(here("data", "input_data", "height.csv"))
leaf_area_raw <- read.csv(here("data", "input_data", "leaf_area.csv"))
leaf_bm_raw <- read.csv(here("data", "input_data", "leaf_biomass.csv"))
c_germ_raw <- read.csv(here("data", "input_data", "cumulative_germination.csv"))

# data cleaning ----------------------------------------------------------------

## community-weighted mean SLA -------------------------------------------------

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

# calculate community-weighted mean SLA values
# for the resident community

cwm_sla_res <- FD::dbFD(
  x = sla_summary2, a = as.matrix(biomass_ab_res_rel2)
)

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

# save to disk -----------------------------------------------------------------

## invader biomass -------------------------------------------------------------

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

## community-weighted mean SLA -------------------------------------------------

write.csv(
  cwm_res, 
  file = here("data", "intermediate_data", "cwm_res.csv"),
  row.names = FALSE
)

## cumulative percentage germination -------------------------------------------

write.csv(
  c_germ, 
  file = here("data", "intermediate_data", "c_germ_week.csv")
)
