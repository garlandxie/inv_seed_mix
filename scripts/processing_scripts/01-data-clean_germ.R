# libraries --------------------------------------------------------------------
library(dplyr)      # for manipulating data
library(here)       # for creating relative file paths
library(GerminaR)   # for calculating germination traits

# import -----------------------------------------------------------------------
c_germ_raw <- read.csv(here("data", "input_data", "cumulative_germination.csv"))

# data clean -------------------------------------------------------------------

## cumulative percentage germination -------------------------------------------

### invasive species ------------------------------------------------------------
c_germ_perc_ciar <- c_germ_raw %>%
  filter(spp == "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ_perc, names_from = spp) %>%
  ungroup() %>%
  rename(cum_germ_perc_ciar = CIAR)

### native species --------------------------------------------------------------
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

## mean germination rate -------------------------------------------------------

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
  mutate(c_gs2, mgt = ger_MGT(evalName = "W", data = .),
  ) 

# calculate germination speed and mean germination time per tray
# for the resident community 
c_gs_res <- c_gs2 %>%
  dplyr::filter(!(spp %in% c("CIAR"))) %>%
  group_by(richness_id, density_id, rep) %>%
  summarize(
    mean_mgt = mean(mgt, na.rm = TRUE),
    mean_mgr = 1/mean_mgt
  ) %>%
  ungroup() 

# data joins -------------------------------------------------------------------

# create a dataset with mean germination rate + percent germination 
c_germ <- left_join(
  c_germ_w12, 
  c_gs_res, 
  by = c("richness_id", "density_id", "rep")
)

# save to disk -----------------------------------------------------------------

write.csv(
  c_germ, 
  file = here("data", "intermediate_data", "c_germ.csv")
)