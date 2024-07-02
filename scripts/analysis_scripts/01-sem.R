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

## root biomass (invader) <- root biomass (resident) ---------------------------

lm_inv_roots_bm <- lm(inv_roots_bm_g ~ res_root_bm_g, 
   data = sem_df)

summary(lm_inv_roots_bm)

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
piecewiseSEM::dSep(sem_germ, conserve = TRUE)

# revised models: resident community -------------------------------------------

# add missing pathways from d-sep tests sequentially 

## 1.RGR height (invader) <- .~. + height (resident) ---------------------------

lm_inv_height2 <- update(lm_inv_height, . ~ . + mean_rgr_height_res)

## 2.CWM SLA (resident) <- .~. germinability (resident) ------------------------

glm_res_sla2 <- update(glm_res_sla, .~. + cum_germ_perc_res)

## 3. Germination rate <- .~. height (resident) --------------------------------

glm_res_mgr2 <- update(glm_res_mgr, .~. + mean_rgr_height_res)

## 4. CWM SLA (resident) <- .~. species richness (resident) --------------------

glm_res_sla3 <- update(glm_res_sla2, .~. + realized_sr)

## 5. Percent germination <- .~. height (resident) -----------------------------

glm_res_germ2 <- update(glm_res_germ, .~. + mean_rgr_height_res)

# run piecewise SEM 2-----------------------------------------------------------

rev_sem_germ <- piecewiseSEM::psem(
  
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

rev_summary_psem <- summary(rev_sem_germ, conserve = TRUE) 
piecewiseSEM::dSep(rev_sem_germ, conserve = TRUE)

View(rev_summary_psem$coefficients)
