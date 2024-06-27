# libraries --------------------------------------------------------------------
library(here)      # for creating relative file-paths
library(dplyr)     # for manipulating data 
library(janitor)   # for cleaning column names

# import -----------------------------------------------------------------------

rgr_height <- read.csv(here("data", "intermediate_data", "rgr_height.csv"))
cum_germ_perc <- read.csv(here("data", "intermediate_data", "c_germ.csv"))
biomass <- read.csv(here("data", "intermediate_data", "biomass.csv"))
cwm_res <- read.csv(here("data", "intermediate_data", "cwm_res.csv"))

# perform inner joins ----------------------------------------------------------

multi_key_id <- c("richness_id", "density_id", "rep")

sem_df <- biomass %>%
  janitor::clean_names() %>%
  full_join(clean_names(rgr_height), by = multi_key_id) %>%
  full_join(clean_names(cum_germ_perc), by = multi_key_id) %>%
  full_join(clean_names(biomass), by = multi_key_id) %>%
  full_join(clean_names(cwm_res), by = multi_key_id) %>%
  dplyr::select(
    
    # experimental treatments
    richness_id, 
    density_id, 
    rep, 
    sown_seeds_res = sown_seeds,
    
    # invader 
    inv_roots_bm_g = inv_roots_bm_g.x,
    inv_abg_bm_g = inv_abg_bm_g.x,
    cum_germ_perc_ciar,
    mean_rgr_height_ciar,
    
    # resident community 
    res_root_bm_g = res_root_bm_g.x,
    res_abg_bm_g = res_abg_bm_g.x, 
    realized_sr = realized_sr.x, 
    mean_rgr_height_res, 
    cum_germ_perc_res, 
    mean_mgr, 
    wds_sla
  ) %>%
  dplyr::filter(density_id != "-") %>%
  mutate(
    richness_id = factor(richness_id, levels = c("M1", "M2", "M4")),
    density_id = factor(density_id, levels = c("D1", "D2", "D3")),
    )

# save to disk -----------------------------------------------------------------

write.csv(
  x = sem_df, 
  file = here("data", "analysis_data", "sem.csv")
)