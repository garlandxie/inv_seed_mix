# libraries --------------------------------------------------------------------
library(dplyr)          # for manipulating data
library(here)           # for creating relative file paths
library(ggplot2)        # for visualizing data
library(emmeans)        # for doing pairwise comparisons
library(stringr)        # for manipulating string characters
library(patchwork)      # for creating multipanel figures

# import -----------------------------------------------------------------------

sem_df <- read.csv(here("data", "analysis_data", "sem.csv"))

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

## invader: biomass ------------------------------------------------------------

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

## statistical analysis: regression --------------------------------------------

lm_res_ab_bm <- lm(res_root_bm_g ~ res_abg_bm_g, data = sem_df)
plot(lm_res_ab_bm)

# relationship between aboveground biomass of resident community
# and invader root biomass 
gam_inv_res_bm <- mgcv::gam(inv_roots_bm_g ~ s(res_abg_bm_g), data = sem_df)
lm_inv_res_bm <- lm(inv_roots_bm_g ~ res_abg_bm_g, data = sem_df)

AIC(gam_inv_res_bm, lm_inv_res_bm)

gam.check(gam_inv_res_bm)

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

plot_res_ab_vs_res_bm <- sem_df %>%
  ggplot(aes(x = res_abg_bm_g, y = res_root_bm_g)) + 
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) + 
  labs(
    x = "Aboveground biomass of native community (in grams)", 
    y = "Root biomass of native community (in grams)"
    ) + 
  theme_bw()

plot_res_bm_vs_inv_root <- sem_df %>%
  ggplot(aes(x = res_abg_bm_g, y = inv_roots_bm_g)) + 
  geom_point() + 
  geom_smooth(method = "gam", se = FALSE) +
  labs(
    x = "Aboveground biomass of native community (in grams)",
    y = "Invader root biomass (in grams)") + 
  theme_bw()

plot_res_bm <- plot_res_ab_vs_res_bm + plot_res_bm_vs_inv_root

## save to disk ----------------------------------------------------------------

ggsave(
  plot = plot_spp_identity, 
  filename = here("output", "results", "mono_ab.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 8
)


  
  


  
