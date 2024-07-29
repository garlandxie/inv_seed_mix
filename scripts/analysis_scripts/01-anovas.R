# libraries --------------------------------------------------------------------
library(dplyr)          # for manipulating data
library(here)           # for creating relative file paths
library(ggplot2)        # for visualizing data
library(emmeans)        # for doing pairwise comparisons
library(stringr)        # for manipulating string characters
library(patchwork)      # for creating multi-panel figures
library(mgcv)           # for test non-linear relationships 

# import -----------------------------------------------------------------------

sem_df <- read.csv(here("data", "analysis_data", "sem.csv"))
bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))

# clean data -------------------------------------------------------------------

mono_ab <- sem_df %>%
  dplyr::filter(richness_id == "M1") %>%
  mutate(
    
    # codes for species identity
    # O = Oenothera biennis, 
    # M = Monarda fistulosa, 
    # A = Andropogon gerardii, 
    # H = Heliopsis helianthoides, 
    # R = Rudbeckia hirta
    spp_identity = substr(rep, start = 1, stop = 1),
    rep_spp = substr(rep, start = 2, stop = 2),
    
    # to stabilize variance in residual plots
    log_inv_bm = log(inv_abg_bm_g),
    
    # convert into factors
    density_id = factor(density_id),
    spp_identity = factor(spp_identity)
  ) %>%
  select(
    density_id, 
    spp_identity, 
    rep_spp, 
    sown_seeds_res, 
    log_inv_bm, 
    inv_abg_bm_g, 
    cum_germ_perc_ciar,
    cum_germ_perc_res,
    mean_mgr,
    res_abg_bm_g,
    res_root_bm_g,
    inv_roots_bm_g) 

mono_ab_inv_bm <- sem_df %>%
  dplyr::filter(richness_id == "M1") %>%
  mutate(
    
    # codes for species identity
    # O = Oenothera biennis, 
    # M = Monarda fistulosa, 
    # A = Andropogon gerardii, 
    # H = Heliopsis helianthoides, 
    # R = Rudbeckia hirta
    spp_identity = substr(rep, start = 1, stop = 1),
    rep_spp = substr(rep, start = 2, stop = 2),
    
    # convert into factors
    density_id = factor(density_id),
    spp_identity = factor(spp_identity)
  ) %>%
  dplyr::select(
    density_id, spp_identity, rep_spp, inv_abg_bm_g
  )

bm_tidy <- bm_raw %>%
  dplyr::filter(Density_ID == "-" & Biomass_Type == "AB") %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Density_ID, Rep) %>%
  summarize(biomass = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Density_ID = case_when(
      Density_ID == "-" ~ "D0", 
      TRUE ~ Density_ID
    ), 
    Rep = substr(Rep, start = 2, stop = 2),
    spp_identity = "None"
    ) %>%
  dplyr::select(
    density_id = Density_ID, 
    spp_identity, 
    rep_spp = Rep, 
    inv_abg_bm_g = biomass
  )

inv_bm <- rbind(bm_tidy, mono_ab_inv_bm)

lm_inv_bm_int <- lm(inv_abg_bm_g ~ spp_identity*density_id, data = inv_bm)

aov_inv_bm <- aov(lm_inv_bm)

emm_inv_bm <- emmeans(lm_inv_bm, ~"spp_identity")

# statistical analyses: ANOVAs -------------------------------------------------

## native: percent seedling emergence ------------------------------------------

lm_mono_perc_res <- glm(cum_germ_perc_res ~ spp_identity*density_id, 
                        family = binomial(link = "logit"), 
                        weight = sown_seeds_res, 
                        data = mono_ab)
plot(lm_mono_perc_res)

aov_mono_perc_res <- aov(lm_mono_perc_res)
summary(aov_mono_perc_res)

emm_mono_perc_res <- emmeans(lm_mono_perc_res, spec = "spp_identity")
pairs_mono_perc_res <- pairs(emm_mono_perc_res)

# get backtransformed emmeans
emm_mono_perc_res_df <- data.frame(emm_mono_perc_res)
emm_mono_perc_res_df <- emm_mono_perc_res_df %>%
  mutate(
    back_emmean = boot::inv.logit(emmean), 
    back_lcl    = boot::inv.logit(asymp.LCL),
    upper_lcl   = boot::inv.logit(asymp.UCL)
  )

cld_perc_res <- multcomp::cld(emm_mono_perc_res)
cld_perc_res_df <- cld_perc_res %>%
  data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    group = as.character(group),
    group = stringr::str_replace(group, pattern = "1", replace = "A"),
    group = stringr::str_replace(group, pattern = "2", replace = "B"),
    group = stringr::str_replace(group, pattern = "3", replace = "C"),
    group = stringr::str_replace(group, pattern = "4", replace = "D"),
    group = stringr::str_replace(group, pattern = "5", replace = "E"),
    group = stringr::str_trim(group), 
    spp_identity = forcats::fct_reorder(spp_identity, emmean, .desc = FALSE)
  ) 

## native: mean germination rate -----------------------------------------------

lm_mono_mgr <- lm(mean_mgr ~ spp_identity*density_id, data = mono_ab)
plot(lm_mono_mgr)

aov_mono_mgr <- aov(lm_mono_mgr)
summary(aov_mono_mgr)

## native: aboveground biomass -------------------------------------------------

lm_mono_res_ab <- lm(res_abg_bm_g ~ spp_identity*density_id, data = mono_ab)
plot(lm_mono_res_ab)

aov_mono_res_ab <- aov(lm_mono_res_ab)
summary(aov_mono_res_ab)

emm_mono_res_ab <- emmeans(lm_mono_res_ab, ~spp_identity*density_id)
multcomp::cld(emm_mono_res_ab)

## invader: percent seedling emergence -----------------------------------------

lm_mono_perc <- lm(cum_germ_perc_ciar ~ spp_identity*density_id, data = mono_ab)
plot(lm_mono_perc)

aov_mono_ab <- aov(lm_mono_perc)
summary(aov_mono_ab)

## monocultures: invader biomass -----------------------------------------------

# model fit
lm_mono_ab <- lm(log_inv_bm ~ spp_identity*density_id, data = mono_ab) 
plot(lm_mono_ab)

# conduct global significance for interaction of predictor variables 
aov_mono_ab <- aov(lm_mono_ab)

# run pairwise comparisons
emm_spp_id <- emmeans(lm_mono_ab, specs = "spp_identity")
pairs_spp_id <- pairs(emm_spp_id)

emm_density_id <- emmeans(lm_mono_ab, specs = "density_id")   
pairs_density_id <- pairs(emm_density_id)

# obtain compact letter displays
cld_spp_id <- multcomp::cld(emm_spp_id)
cld_spp_id_df <- cld_spp_id %>%
  data.frame() %>%
  janitor::clean_names() %>%
  mutate(
    back_emm = exp(emmean),
    back_low_cl = exp(lower_cl), 
    back_upper_cl = exp(upper_cl), 
    group = as.character(group), 
    group = stringr::str_replace(group, pattern = "1", replace = "A"),
    group = stringr::str_replace(group, pattern = "2", replace = "B"),
    group = stringr::str_replace(group, pattern = "3", replace = "C"),
    group = stringr::str_replace(group, pattern = "4", replace = "D"),
    group = stringr::str_trim(group), 
    spp_identity = forcats::fct_reorder(spp_identity, back_emm, .desc = FALSE)
    ) 

## invader: root biomass -------------------------------------------------------

# relationship between aboveground and root biomass of resident community
lm_res_ab_bm <- lm(res_root_bm_g ~ res_abg_bm_g, data = sem_df)
plot(lm_res_ab_bm)

# relationship between aboveground biomass of resident community
# and invader root biomass 
gam_inv_res_bm <- mgcv::gam(inv_roots_bm_g ~ s(res_abg_bm_g), data = sem_df)
lm_inv_res_bm <- lm(inv_roots_bm_g ~ res_abg_bm_g, data = sem_df)

AIC(gam_inv_res_bm, lm_inv_res_bm)

par(mfrow = c(2,2))
mgcv::gam.check(gam_inv_res_bm)

# visualize data ---------------------------------------------------------------

## ANOVA: native germinability -------------------------------------------------
(plot_perc_res <- cld_perc_res_df %>%
   ggplot(aes(x = spp_identity, y = emmean)) + 
   geom_point() + 
   geom_pointrange(aes(x = spp_identity, y = emmean, ymin = asymp_lcl , ymax = asymp_ucl)) + 
   geom_text(aes(label = group), size = 3, nudge_x = 0.17, nudge_y = 0.01) + 
   scale_x_discrete(labels = c(
     expression(italic("Oenothera biennis")), 
     expression(italic("Rudbeckia hirta")), 
     expression(italic("Monardo fistulosa")), 
     expression(italic("Andropogon gerardii")), 
     expression(italic("Heliopsis helianthoides"))
     )
   ) + 
   coord_flip() + 
   labs(
     x = "Species Identity (in monocultures)", 
     y = expression(paste("Logit(Germinability)"))
   ) + 
   theme_bw()
)

## ANOVA: invader biomass ------------------------------------------------------
(plot_spp_id <- cld_spp_id_df %>%
  ggplot(aes(x = spp_identity, y = emmean)) + 
  geom_point() + 
  geom_pointrange(aes(x = spp_identity, y = emmean, ymin = lower_cl, ymax = upper_cl)) + 
  geom_text(aes(label = group), size = 3, nudge_x = 0.17, nudge_y = 0.1) + 
  scale_x_discrete(labels = c(
    "Oenothera biennis", 
     "Rudbeckia hirta", 
     "Monardo fistulosa", 
     "Andropogon gerardii", 
     "Heliopsis helianthoides")
    ) + 
  coord_flip() + 
  labs(
    x = NULL, 
    y = expression(paste("Log(Biomass of ", italic("C. arvense"), " [grams])"))
    ) + 
  theme_bw() + 
  theme(axis.text.y=element_blank())
)

## ANOVA: multi-panel figure ---------------------------------------------------

(plot_spp_identity <- plot_perc_res + plot_spp_id)

## root biomass ----------------------------------------------------------------

# get summary statistics from lm model object
rsq_res_ab_bm <- as.character(round(summary(lm_res_ab_bm)$adj.r.squared, digits = 2))
p_res_ab_bm <- summary(lm_res_ab_bm)$coefficients["res_abg_bm_g", "Pr(>|t|)"]
p_res_ab_bm <- ifelse(p_res_ab_bm<0.001, "p<0.001")

(plot_res_ab_vs_res_bm <- sem_df %>%
  ggplot(aes(x = res_abg_bm_g, y = res_root_bm_g)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  geom_text(label = "A)", x = 0.1, y = 17.5) + 
  geom_text(
    label = paste("~adj-R^{2} ==",  rsq_res_ab_bm), 
    x = 11.2, 
    y = 16, 
    parse = TRUE
  ) + 
  geom_text(label = paste("OLS:", p_res_ab_bm), x = 11.5, y = 17) + 
  labs(
    x = "Aboveground biomass of native community (in grams)", 
    y = "Root biomass of native community (in grams)"
    ) + 
  theme_bw()
)

# get summary statistics from gam model object
rsq_gam_inv_bm <- as.character(round(summary(gam_inv_res_bm)$r.sq, digits = 2))
edf_gam_inv_bm <- as.character(round(summary(gam_inv_res_bm)$edf, digits = 2))
p_gam_inv_bm <- summary(gam_inv_res_bm)$s.table[, "p-value"]
p_gam_inv_bm <- ifelse(p_gam_inv_bm < 0.001, "p<0.001")  
  
(plot_res_bm_vs_inv_root <- sem_df %>%
  ggplot(aes(x = res_abg_bm_g, y = inv_roots_bm_g)) + 
  geom_point() + 
  geom_smooth(method = "gam", se = FALSE) +
  geom_text(label = "B)", x = 0.2, y= 7.1) + 
  geom_text(label = paste("GAM:", "edf = ", edf_gam_inv_bm, ",", p_gam_inv_bm), x = 10, y = 7) + 
  geom_text(label = paste("adj-R^{2} == ", rsq_gam_inv_bm), x = 10, y = 6.5, parse = TRUE) + 
  labs(
    x = "Aboveground biomass of native community (in grams)",
    y = expression(paste("Root biomass of ", italic("Cirsium arvense")," (in grams)"))
    ) + 
  theme_bw()
)

(plot_inv_roots <- plot_res_ab_vs_res_bm + plot_res_bm_vs_inv_root)

## save to disk ----------------------------------------------------------------

ggsave(
  plot = plot_spp_identity, 
  filename = here("output", "results", "mono_ab.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 8
)

ggsave(
  plot = plot_inv_roots, 
  filename = here("output", "results", "inv_roots.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 9
)



  
  


  
