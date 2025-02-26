# libraries --------------------------------------------------------------------
library(dplyr)          # for manipulating data
library(here)           # for creating relative file paths
library(ggplot2)        # for visualizing data
library(stringr)        # for manipulating string characters
library(patchwork)      # for creating multi-panel figures
library(mgcv)           # for test non-linear relationships 
library(broom)          # for tidying lm objects

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

mix_inv_bm <- sem_df %>%
  dplyr::filter(richness_id %in% c("M2", "M4")) %>%
  select(density_id, spp_identity = richness_id, rep_spp = rep, inv_abg_bm_g)


inv_bm <- rbind(bm_tidy, mono_ab_inv_bm, mix_inv_bm)

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

## invader: aboveground biomass ------------------------------------------------

# clean data 
inv_bm_d0_d1 <- inv_bm %>%
  dplyr::filter(density_id %in% c("D0", "D1")) %>%
  mutate(spp_identity = factor(
    spp_identity, levels = c("None","O", "A", "R", "M", "H", "M2", "M4"))
    )

inv_bm_d0_d2 <- inv_bm %>%
  dplyr::filter(density_id %in% c("D0", "D2")) %>%
  mutate(spp_identity = factor(
    spp_identity, levels = c("None","O", "A", "R", "M", "H", "M2", "M4"))
  )

inv_bm_d0_d3 <- inv_bm %>%
  dplyr::filter(density_id %in% c("D0", "D3")) %>%
  mutate(spp_identity = factor(
    spp_identity, levels = c("None","O", "A", "R", "M", "H", "M2", "M4"))
  )

# model fit
lm_inv_bm_d0_d1 <- lm(inv_abg_bm_g ~ spp_identity, data = inv_bm_d0_d1) 
lm_inv_bm_d0_d2 <- lm(log(inv_abg_bm_g) ~ spp_identity, data = inv_bm_d0_d2) 
lm_inv_bm_d0_d3 <- lm(log(inv_abg_bm_g) ~ spp_identity, data = inv_bm_d0_d3) 

# sanity checks
plot(lm_inv_bm_d0_d1)
plot(lm_inv_bm_d0_d2)
plot(lm_inv_bm_d0_d3)

# run two-way ANOVA
aov_inv_bm_d0_d1 <- summary(aov(lm_inv_bm_d0_d1))
aov_inv_bm_d0_d2 <- summary(aov(lm_inv_bm_d0_d2))
aov_inv_bm_d0_d3 <- summary(aov(lm_inv_bm_d0_d3))

# visualize the data 
broom_inv_bm_d0_d1 <- lm_inv_bm_d0_d1 %>%
  broom::tidy() %>%
  mutate(density_id = "D1")

broom_inv_bm_d0_d2 <- lm_inv_bm_d0_d2 %>%
  broom::tidy() %>%
  mutate(density_id = "D2")

broom_inv_bm_d0_d3 <- lm_inv_bm_d0_d3 %>%
  broom::tidy() %>%
  mutate(density_id = "D3")

ci_inv_bm_d0_d1 <- lm_inv_bm_d0_d1 %>%
  confint() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(density_id = "D1")

ci_inv_bm_d0_d2 <- lm_inv_bm_d0_d2 %>%
  confint() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(density_id = "D2")

ci_inv_bm_d0_d3 <- lm_inv_bm_d0_d3 %>%
  confint() %>%
  as.data.frame() %>%
  tibble::rownames_to_column(var = "term") %>%
  mutate(density_id = "D3")

ci_inv_bm <- rbind(
  ci_inv_bm_d0_d1, 
  ci_inv_bm_d0_d2, 
  ci_inv_bm_d0_d3
)

broom_inv_bm <- rbind(
  broom_inv_bm_d0_d1, 
  broom_inv_bm_d0_d2, 
  broom_inv_bm_d0_d3
  )

tidy_inv_bm <- inner_join(broom_inv_bm, ci_inv_bm, by = c("term", "density_id"))

emm_inv_bm_d0_d1 <- emmeans(lm_inv_bm_d0_d1, ~"spp_identity")
contr_inv_bm_d0_d1 <- contrast(emm_inv_bm_d0_d1, method = "trt.vs.ctrl")

emm_inv_bm_d0_d2 <- emmeans(lm_inv_bm_d0_d2, ~"spp_identity")
contr_inv_bm_d0_d2 <- contrast(emm_inv_bm_d0_d2, method = "trt.vs.ctrl")
 
emm_inv_bm_d0_d3 <- emmeans(lm_inv_bm_d0_d3, ~"spp_identity")
contr_inv_bm_d0_d3 <- contrast(emm_inv_bm_d0_d3, method = "trt.vs.ctrl")
  
## mixtures: invader biomass ---------------------------------------------------

# clean data 
mixtures <- sem_df %>%
  dplyr::filter(richness_id %in% c("M2", "M4")) %>%
  select(density_id, spp_identity = richness_id, rep_spp = rep, inv_abg_bm_g) %>%
  rbind(bm_tidy) %>%
  mutate(spp_identity = factor(spp_identity, levels = c("None", "M2", "M4")))
    
# model fit    
inv_bm_mix_d0_d1 <- dplyr::filter(mixtures, density_id %in% c("D0", "D1"))
inv_bm_mix_d0_d2 <- dplyr::filter(mixtures, density_id %in% c("D0", "D2"))
inv_bm_mix_d0_d3 <- dplyr::filter(mixtures, density_id %in% c("D0", "D3"))

# sanity checks
lm_inv_mix_d0_d1 <- lm(inv_abg_bm_g ~ spp_identity, data = inv_bm_mix_d0_d1) 
lm_inv_mix_d0_d2 <- lm(inv_abg_bm_g ~ spp_identity, data = inv_bm_mix_d0_d2) 
lm_inv_mix_d0_d3 <- lm(inv_abg_bm_g ~ spp_identity, data = inv_bm_mix_d0_d3) 

## invader: root biomass -------------------------------------------------------

# relationship between aboveground and root biomass of resident community
lm_res_ab_bm <- lm(res_root_bm_g ~ res_abg_bm_g, data = sem_df)
plot(lm_res_ab_bm)

# relationship between aboveground biomass of resident community
# and invader root biomass 
gam_inv_res_bm <- mgcv::gam(inv_roots_bm_g ~ s(res_abg_bm_g), data = sem_df)
lm_inv_res_bm <- lm(inv_roots_bm_g ~ res_abg_bm_g, data = sem_df)

AIC(gam_inv_res_bm, lm_inv_res_bm)
summary(gam_inv_res_bm)
coef(gam_inv_res_bm)

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

# summary statistics -----------------------------------------------------------

summary_inv_bm <- inv_bm %>%
  group_by(density_id, spp_identity) %>%
  summarize(
    mean_inv_bm_g = mean(inv_abg_bm_g, na.rm = TRUE) %>% round(digits = 2),
    sd_inv_bm_g = sd(inv_abg_bm_g, na.rm = TRUE) %>% round(digits = 2)
    ) %>%
  ungroup()

perc_diff_inv_bm <- summary_inv_bm %>%
  dplyr::filter(density_id %in% c("D0", "D3")) %>%
  select(spp_identity,mean_inv_bm_g) %>%
  tidyr::pivot_wider(values_from = mean_inv_bm_g, names_from = spp_identity)

percent_diff <- function(old_value, new_value) {
  ((new_value - old_value) / mean(old_value, new_value)) * 100
}

percent_diff(perc_diff_inv_bm$None, perc_diff_inv_bm$A)

# save to disk -----------------------------------------------------------------

## contrasts -------------------------------------------------------------------

# prep the data frame
contr_inv_bm_d0_d1 <- as.data.frame(contr_inv_bm_d0_d1)
contr_inv_bm_d0_d1$density_id <- "D1"
contr_inv_bm_d0_d1$estimate <- round(contr_inv_bm_d0_d1$estimate, digits = 2)

contr_inv_bm_d0_d2 <- as.data.frame(contr_inv_bm_d0_d2)
contr_inv_bm_d0_d2$density_id <- "D2"
contr_inv_bm_d0_d2$estimate <- round(contr_inv_bm_d0_d2$estimate, digits = 2)

contr_inv_bm_d0_d3 <- as.data.frame(contr_inv_bm_d0_d3)
contr_inv_bm_d0_d3$density_id <- "D3"
contr_inv_bm_d0_d3$estimate <- round(contr_inv_bm_d0_d3$estimate, digits = 2)

contr_inv_bm <- rbind(
  contr_inv_bm_d0_d1, 
  contr_inv_bm_d0_d2, 
  contr_inv_bm_d0_d3
  )

write.csv(
  x = contr_inv_bm, 
  file = here("output", "data_appendix_output", "table_contrasts.csv"),
  row.names = FALSE
)

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



  
  


  
