# libraries --------------------------------------------------------------------
library(dplyr)          # for manipulating data
library(here)           # for creating relative file paths
library(piecewiseSEM)   # for calculating piecewise structural equation models

# import -----------------------------------------------------------------------

sem_df <- read.csv(here("data", "analysis_data", "sem.csv"))

# statistical analyses: resident community -------------------------------------

## % germination (resident) ----------------------------------------------------

# model fit
glm_res_germ <- glm(
  cum_germ_perc_res ~  mean_mgr, 
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df)

summary(glm_res_germ)

## germination speed (resident) <- richness + density --------------------------

# model fit
glm_res_mgr <- glm(
  mean_mgr ~ density_id + richness_id, 
  family = Gamma(link = "log"), 
  data = sem_df)

## RGR height (resident) <- richness + density ---------------------------------

# model fit 
lm_res_height <- lm(
  mean_rgr_height_res ~ richness_id + density_id, 
  data = sem_df)
 
summary(lm_res_height) 

## abg biomass (resident) <- % germination (resident) + height (resident) ------

# model fit
lm_res_abg_bm <- lm(
  res_abg_bm_g ~ cum_germ_perc_res + mean_rgr_height_res + realized_sr, 
  data = sem_df)

summary(lm_res_bm)

## species richness (resident) <- % germination (resident) ---------------------

lm_res_sr <- lm(realized_sr ~ cum_germ_perc_res + richness_id, data = sem_df)

summary(lm_sr)

## root biomass (resident) <- % germination ------------------------------------

lm_res_root_bm <- lm(
  res_root_bm_g ~ cum_germ_perc_res, 
  data = sem_df
)

summary(lm_res_root_bm)

## CWM SLA (resident) <- community biomass (resident) --------------------------

# model fit
glm_res_sla <- lm(
  wds_sla ~ res_abg_bm_g,
  data = sem_df
)

summary(glm_res_sla)

# statistical analyses: invader =-----------------------------------------------

## % germination (invader) <- community biomass (resident) ---------------------

# model fit
glm_inv_germ <- glm(
  cum_germ_perc_ciar ~ res_abg_bm_g + realized_sr, 
  family = binomial(link = "logit"),
  data = sem_df
)

summary(glm_inv_germ)

## RGR height (invader) <- community biomass (resident) ------------------------

# model fit
lm_inv_height <- lm(
  mean_rgr_height_ciar ~ 
    res_abg_bm_g + realized_sr, 
  data = sem_df
)

summary(lm_inv_height)

## abg biomass (invader) <- % germination (invader) + height (invader) + SLA ---

# model fit
lm_inv_bm <- lm(
  inv_abg_bm_g ~ 
    cum_germ_perc_ciar + 
    mean_rgr_height_ciar + 
    wds_sla,
  data = sem_df
) 

summary(lm_inv_bm)

# run piecewise SEM 1-----------------------------------------------------------

sem_germ <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ,
  glm_res_mgr,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height,
  lm_inv_bm
)

summary_psem <- summary(sem_germ, conserve = TRUE) 
dSep1 <- sem_germ %>%
  piecewiseSEM::dSep(conserve = TRUE) %>%
  data.frame() %>%
  mutate(run = 1) 

# revised models: resident community -------------------------------------------

# add missing pathways from d-sep tests sequentially 

## 1.RGR height (invader) <- .~. + height (resident) ---------------------------

lm_inv_height2 <- update(lm_inv_height, . ~ . + mean_rgr_height_res)

## 2. run piecewise SEM 2-----------------------------------------------------------

rev_sem_germ2 <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ,
  glm_res_mgr,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height2,
  lm_inv_bm
)

rev_summary_psem2 <- summary(rev_sem_germ2, conserve = TRUE) 
dSep2 <- rev_sem_germ2 %>%
  piecewiseSEM::dSep(conserve = TRUE) %>%
  data.frame() %>%
  mutate(run = 2)

## 3. CWM SLA (resident) <- .~. germinability (resident) ------------------------

glm_res_sla2 <- update(glm_res_sla, .~. + cum_germ_perc_res)

## 4. run piecewise SEM 3--------------------------------------------------------

rev_sem_germ3 <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ,
  glm_res_mgr,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla2,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height2,
  lm_inv_bm
)

rev_summary_psem3 <- summary(rev_sem_germ3, conserve = TRUE) 
dSep3 <- rev_sem_germ3 %>%
  piecewiseSEM::dSep(conserve = TRUE) %>%
  data.frame() %>%
  mutate(run = 3)

## 5. Germination rate (resident) <- .~. height (resident) ---------------------

glm_res_mgr2 <- update(glm_res_mgr, .~. + mean_rgr_height_res)

## 6. run piecewise SEM 4 -------------------------------------------------------

rev_sem_germ4 <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ,
  glm_res_mgr2,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla2,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height2,
  lm_inv_bm
)

rev_summary_psem4 <- summary(rev_sem_germ4, conserve = TRUE) 
dSep4 <- rev_sem_germ4 %>%
  piecewiseSEM::dSep(conserve = TRUE) %>%
  data.frame() %>%
  mutate(run = 4)

## 7. CWM SLA (resident) <- .~. species richness (resident) --------------------

glm_res_sla3 <- update(glm_res_sla2, .~. + realized_sr)

## 8. run piecewise SEM 5 ------------------------------------------------------

rev_sem_germ5 <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ,
  glm_res_mgr2,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla3,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height2,
  lm_inv_bm
)

rev_summary_psem5 <- summary(rev_sem_germ5, conserve = TRUE) 
dSep5 <- rev_sem_germ5 %>%
  piecewiseSEM::dSep(conserve = TRUE) %>%
  data.frame() %>%
  mutate(run = 5)

## 9. Percent germination <- .~. height (resident) -----------------------------

glm_res_germ2 <- update(glm_res_germ, .~. + mean_rgr_height_res)

## 10. final piecewise SEM -----------------------------------------------------

rev_sem_germ6 <- piecewiseSEM::psem(
  
  # resident community 
  glm_res_germ2,
  glm_res_mgr2,
  lm_res_height,
  lm_res_abg_bm,
  glm_res_sla3,
  lm_res_sr,
  
  # invader
  glm_inv_germ,
  lm_inv_height2,
  lm_inv_bm
)

rev_summary_psem6 <- summary(rev_sem_germ6, conserve = TRUE) 
stand_coefs <- piecewiseSEM::coefs(rev_sem_germ6, standardize = "range", 
                                   standardize.type = "Menard.OE")

# save to disk -----------------------------------------------------------------

# clean up the d-sep test a bit before saving it to disk
d_sep_df <- dSep1 %>%
  rbind(dSep2) %>%
  rbind(dSep3) %>%
  rbind(dSep4) %>%
  rbind(dSep5) %>%
  janitor::clean_names() %>%
  dplyr::filter(p_value < 0.05)  %>%
  mutate(crit_value = round(crit_value, digits = 2)) %>%
  select(
    run, 
    independ_claim, 
    df, 
    crit_value, 
    var_6)
  
write.csv(
  x = d_sep_df, 
  file = here("output", "data_appendix_output", "table_s2_dsep.csv"), 
  row.names = FALSE
)
