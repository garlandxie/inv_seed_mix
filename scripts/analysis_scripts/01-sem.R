# libraries --------------------------------------------------------------------
library(dplyr)          # for manipulating data
library(here)           # for creating relative file paths
library(piecewiseSEM)   # for calculating piecewise structural equation models
library(ggplot2)        # for obtaining partial residuals
library(patchwork)

# import -----------------------------------------------------------------------

sem_df <- read.csv(here("data", "analysis_data", "sem.csv"))

# statistical analyses: resident community -------------------------------------

## % germination (resident) ----------------------------------------------------

# model fit
glm_res_germ <- glm(
  cum_germ_perc_res ~ mean_mgr,
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df
  )
 
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

summary(lm_res_abg_bm)

## species richness (resident) <- % germination (resident) ---------------------

lm_res_sr <- lm(realized_sr ~ cum_germ_perc_res + richness_id, data = sem_df)

summary(lm_res_sr)

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

# double-check collinearity
car::vif(glm_res_sla3)

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

# data visualization -----------------------------------------------------------

## get partial residuals -------------------------------------------------------

# for wds sla 
visreg_wds_res <- visreg::visreg(
  glm_res_sla3, "res_abg_bm_g", 
  type = "conditional") 

visreg_wds_germ <- visreg::visreg(
  glm_res_sla3, "cum_germ_perc_res", 
  type = "conditional") 

visreg_wds_sr <- visreg::visreg(
  glm_res_sla3, "realized_sr", 
  type = "conditional") 

# for invader biomass
visreg_sla <- visreg::visreg(
  lm_inv_bm, "wds_sla", 
  type= "conditional") 

## plot wds sla ----------------------------------------------------------------

# partial R2
part2_res_sla <- sensemakr::partial_r2(glm_res_sla3)

# p-values for each predictor variable
sum_res_sla <- summary(glm_res_sla3)

p_bm <- sum_res_sla$coefficients["res_abg_bm_g", "Pr(>|t|)"] 
p_bm <- ifelse(p_bm < 0.01, "p<0.01")

p_germ <- sum_res_sla$coefficients["cum_germ_perc_res", "Pr(>|t|)"] 
p_germ <- ifelse(p_germ < 0.001, "p<0.001")

p_sr <- sum_res_sla$coefficients["realized_sr", "Pr(>|t|)"]
p_sr <- ifelse(p_sr < 0.001, "p<0.001")

(plot_wds_abg <- 
  ggplot() + 
   
  # partial residuals
  geom_point(
    aes(x = res_abg_bm_g, y = visregRes), 
    alpha = 0.5, 
    data = visreg_wds_res$res
    ) +
   
   # 95% confidence interval
   geom_line(
     aes(x = res_abg_bm_g, y = visregFit), 
     col = "blue", 
     linewidth = 0.75, 
     data = visreg_wds_res$fit
   ) +
   geom_ribbon(
     aes(x = res_abg_bm_g, y = visregFit,
         ymin = visregLwr, ymax = visregUpr), 
     col = "grey", 
     alpha = 0.1, 
     data = visreg_wds_res$fit
   ) + 
   
  labs(
    x = "Aboveground community biomass (g)",
    y = expression(
      paste("Partial residuals (",
            "| CWM SLA - ", 
            italic("Cirsium arvense"), " |)"
            )
      )
    ) + 
  annotate("text", x = 1,  y = 0.25, label = "A)") + 
  annotate(
    "text", 
    x = 12, 
    y = 0.23, 
    label = paste(
      "partial R2:", 
      round(part2_res_sla["res_abg_bm_g"], digits = 2)
      )
    ) +
  annotate(
    "text", 
    x = 11, 
    y = 0.215, 
    label = p_bm
  ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25)) + 
  xlim(0, 15) + 
  ylim(-0.05, 0.25) + 
  theme_bw()
)

(plot_wds_germ <- 
    ggplot() + 
    
    # partial residuals 
    geom_point(
      aes(x = cum_germ_perc_res, y = visregRes), 
      alpha = 0.5,
      data = visreg_wds_germ$res
    ) + 
    
    # 95% confidence interval
    geom_line(
      aes(x = cum_germ_perc_res, y = visregFit), 
      col = "blue", 
      linewidth = 0.75, 
      data = visreg_wds_germ$fit
    ) +
    geom_ribbon(
      aes(x = cum_germ_perc_res, y = visregFit,
          ymin = visregLwr, ymax = visregUpr), 
      col = "grey", 
      alpha = 0.1, 
      data = visreg_wds_germ$fit
    ) + 
    
  # add annotations
  annotate("text", x = 0.1, y = 0.25, label = "B)") + 
  annotate(
    "text", 
    x = 0.80, 
    y = 0.23, 
    label = paste(
      "partial R2:", 
      round(part2_res_sla["cum_germ_perc_res"], digits = 2)
      )
    ) +
    annotate(
      "text", 
      x = 0.72, 
      y = 0.215, 
      label = p_bm
    ) + 
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25)) + 
  xlim(0, 1) + 
  ylim(-0.05, 0.25) +
  labs(
    x = "Germinability of resident community (%)",
    y = NULL) + 
  theme_bw() + 
  theme(axis.text.y = element_blank())
)

(plot_wds_sr <- 
  ggplot() + 
  
  # partial residuals
  geom_point(
    aes(x = realized_sr, y = visregRes), 
    alpha = 0.5, 
    data = visreg_wds_sr$res
  ) + 
  
  # 95% confidence interval
  geom_line(
    aes(x = realized_sr, y = visregFit), 
    col = "blue", 
    linewidth = 0.65, 
    data = visreg_wds_sr$fit
    ) +
  geom_ribbon(
    aes(x = realized_sr, y = visregFit,
        ymin = visregLwr, ymax = visregUpr), 
    col = "grey", 
    alpha = 0.1, 
    data = visreg_wds_sr$fit
    ) +
  annotate("text", x = 1.2, y = 0.25, label = "C)") + 
  annotate(
    "text", 
    x = 3.4, 
    y = 0.23, 
    label = paste(
      "partial R2:", 
      round(part2_res_sla["realized_sr"], digits = 2)
      )
    ) +
  annotate(
    "text", 
    x = 3.2, 
    y = 0.215, 
    label = p_sr
    ) +
  scale_y_continuous(breaks = c(0, 0.05, 0.10, 0.15, 0.20, 0.25)) + 
  ylim(-0.05, 0.25) + 
  labs(
    x = "Realized species richness",
    y = NULL) + 
  theme_bw() + 
  theme(axis.text.y = element_blank())
)

(plot_wds <- plot_wds_abg + plot_wds_germ + plot_wds_sr)

## plot invader biomass---------------------------------------------------------

# get info for annotation

# partial R2
part2_inv_bm <-sensemakr::partial_r2(lm_inv_bm)
part2_inv_bm["wds_sla"]

# p-value of the full linear model
f_inv_bm <- summary(lm_inv_bm)$fstatistic
p_inv_bm <- pf(
   f_inv_bm["value"], 
   f_inv_bm["numdf"],
   f_inv_bm["dendf"], 
   lower.tail = FALSE
   ) 
p_inv_bm <- ifelse(p_inv_bm < 0.001, "<0.001")

(plot_sla_v_inv <- 
  ggplot() +
    
  # partial residuals
  geom_point(
    aes(x = wds_sla, y = visregRes), 
    data = visreg_sla$res) + 
 
  # 95% confidence interval
    geom_line(
    aes(x = wds_sla, y = visregFit), 
    col = "blue", 
    linewidth = 0.75, 
    data = visreg_sla$fit
  ) +
  geom_ribbon(
    aes(x = wds_sla, y = visregFit,
        ymin = visregLwr, ymax = visregUpr), 
    col = "grey", 
    alpha = 0.1, 
    data = visreg_sla$fit
  ) + 
  annotate(
     "text", 
     x = 0.20, 
     y = 4.9, 
     label = paste("partial R2:", round(part2_inv_bm["wds_sla"], digits = 2))
   ) + 
  annotate(
     "text", 
     x = 0.185, 
     y = 4.6, 
     label = paste0(
       "F(", f_inv_bm["numdf"], ",", f_inv_bm["dendf"], ")", " = ",
       round(f_inv_bm["value"], digits = 2), ", p", p_inv_bm)
   ) + 
  labs(
    x = expression(paste("| CWM SLA - ", italic("Cirsium arvense"), " |")),
    y = "Partial residuals (invader biomass)",
    ) +  
  theme_bw() 
)

## plot early plant stage ------------------------------------------------------

visreg_cum_germ <- visreg::visreg(
  glm_res_germ2, "mean_mgr", 
  type = "conditional"
)

visreg_res_height <- visreg::visreg(
  glm_res_germ2, "mean_rgr_height_res", 
  type = "conditional"
)

(plot_mgr_germ <- ggplot() +
  geom_rug(aes(x = mean_mgr), data = visreg_cum_germ$res) + 
  geom_line(
    aes(x = mean_mgr, y = exp(visregFit)), 
    col = "blue", 
    data = visreg_cum_germ$fit
    ) + 
    geom_ribbon(
      aes(x = mean_mgr, y = exp(visregFit),
          ymin = exp(visregLwr), ymax = exp(visregUpr)
          ), 
      col = "grey", 
      alpha = 0.1, 
      data = visreg_cum_germ$fit
    ) + 
   ylim(0.3, 0.7) + 
   labs(
     x = "Mean Germination Rate", 
     y = "Odds of higher germinability") + 
   theme_bw()
)

visreg_mgr_height <- visreg::visreg(
  glm_res_mgr2, "mean_rgr_height_res", 
  type = "conditional"
)

(plot_mgr <- ggplot() + 
    geom_point(
      aes(x = mean_rgr_height_res, y = visregRes), 
      data = visreg_mgr_height$res
      )  + 
    geom_line(
      aes(x = mean_rgr_height_res, y = visregFit), 
      col = "blue", 
      data = visreg_mgr_height$fit
    ) + 
    geom_ribbon(
      aes(x = mean_rgr_height_res, y = visregFit,
          ymin = visregLwr, ymax = visregUpr
      ), 
      col = "grey", 
      alpha = 0.1, 
      data = visreg_mgr_height$fit
    ) + 
  labs(
    x = "RGR height (resident)",
    y = "Mean germination rate") + 
  theme_bw()
)

(plot_height <- ggplot() + 
    geom_rug(aes(x = mean_rgr_height_res), data = visreg_res_height$res) + 
    geom_line(
      aes(x = mean_rgr_height_res, y = exp(visregFit)), 
      col = "blue", 
      data = visreg_res_height$fit
    ) + 
    geom_ribbon(
      aes(x = mean_rgr_height_res, y = exp(visregFit),
          ymin = exp(visregLwr), ymax = exp(visregUpr)
      ), 
      col = "grey", 
      alpha = 0.1, 
      data = visreg_res_height$fit
    ) + 
    labs(
      x = "RGR height (resident community)",
      y = "Odds of higher germinability") + 
    theme_bw()
)

res_height_emm <- emmeans::emmeans(lm_res_height, "density_id")
multcomp::cld(res_height_emm)

(plot_rich_dens <- sem_df %>% 
    ggplot(aes(x = density_id, y = mean_rgr_height_res)) + 
    geom_boxplot() + 
    geom_jitter(width = 0.1, alpha = 0.1) + 
    labs(
      x = "Seeding density", 
      y = "RGR height (resident)"
      ) + 
    theme_bw() 
)

(plot_early <- plot_mgr_germ + plot_mgr + plot_height + plot_rich_dens)

# save to disk ----------------------------------------------------------------

ggsave(
  plot = plot_sla_v_inv, 
  filename = here("output", "results", "fig_sla_v_inv.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 5
)

ggsave(
  plot = plot_wds, 
  filename = here("output", "results", "fig_wds_sla.png"),
  device = "png", 
  units = "in", 
  height = 4, 
  width = 10
)
