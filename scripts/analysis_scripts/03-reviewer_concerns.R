# libraries --------------------------------------------------------------------
library(dplyr)       # for manipulating data
library(ggplot2)     # for visualizing data
library(here)        # for creating relative file-paths
library(emmeans)     # for doing pairwise comparisons 
library(multcomp)    # for getting compact letter displays
library(patchwork)   # for visualizing multi-figure panels
library(nlme)        # for running generalized least squared models
library(car)         # for running type II ANOVAs

# import data ------------------------------------------------------------------
bm_raw <- read.csv(here("data", "input_data", "biomass.csv"))

# clean data -------------------------------------------------------------------

## species identity within monocultures ----------------------------------------

mono_ab <- bm_raw %>%
  dplyr::filter(
    Biomass_Type == "AB", 
    Rep %in% c(
      "H1", "H2", "H3", "H4", "H5",
      "A1", "A2", "A3", "A4", "A5",
      "M1", "M2", "M3", "M4", "M5", 
      "O1", "O2", "O3", "O4", "O5",
      "R1", "R2", "R3", "R4", "R5"
      ),
    Species_ID == "CIAR",
    Biomass_g != "-"
  ) %>%
  mutate(Biomass_g = as.numeric(Biomass_g)) %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(biomass = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(spp_identity = substring(Rep, first = 1, last = 1))

## aboveground biomass ---------------------------------------------------------

ab_bm_tidy <- bm_raw %>%
  dplyr::filter(Density_ID != "-" & Biomass_Type == "AB") %>%
  dplyr::filter(Species_ID == "CIAR") %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Density_ID, Richness_ID, Rep) %>%
  summarize(biomass = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Richness_ID = factor(Richness_ID, levels = c("M1", "M2", "M4")),
    Density_ID = factor(Density_ID, levels = c("D1", "D2", "D3"))
  )

## belowground biomass ---------------------------------------------------------

bg_bm_tidy <- bm_raw %>%
  dplyr::filter(Density_ID != "-" & Biomass_Type == "BM") %>%
  dplyr::filter(Species_ID == "CIAR") %>%
  mutate(
    Biomass_g = na_if(Biomass_g, "-"),
    Biomass_g = as.numeric(Biomass_g)
  ) %>%
  group_by(Density_ID, Richness_ID, Rep) %>%
  summarize(biomass = sum(Biomass_g, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    Richness_ID = factor(Richness_ID, levels = c("M1", "M2", "M4")),
    Density_ID = factor(Density_ID, levels = c("D1", "D2", "D3"))
  )

# ANOVAs -----------------------------------------------------------------------

## species identity within monocultures ----------------------------------------

# apply generalized least-squares to account for heteroscedastic errors
# include interactions between species identity and seeding density
lm_id_dens <- nlme::gls(
  biomass ~ spp_identity*Density_ID, 
  weights = varPower(),
  data = mono_ab
)

# double-check residual distributions
plot(lm_id_dens) 

# fit type II Anova
aov_id_dens <- car::Anova(lm_id_dens, type = "II")

# get estimated marginal means
emm_id_by_dens  <- emmeans(lm_id_dens, pairwise ~ spp_identity | Density_ID)
emm_dens_by_id  <- emmeans(lm_id_dens, pairwise ~ Density_ID | spp_identity)

# double-check interaction plots
emmip(emm_id_dens, spp_identity ~ Density_ID)
emmip(emm_id_dens, Density_ID ~ spp_identity)

# get compact letter displays
cld_id_by_dens <- cld(emm_id_by_dens)
cld_dens_by_id <- cld(emm_dens_by_id)

## aboveground biomass ---------------------------------------------------------

lm_abg_bm_rich_den <- lm(biomass ~ Density_ID*Richness_ID, data = ab_bm_tidy)
plot(lm_abg_bm_rich_den)
car::ncvTest(lm_abg_bm_rich_den)

gls_abg_bm_rich_den <- gls(biomass ~ Density_ID*Richness_ID, weights = varPower(), data = ab_bm_tidy)
aov_abg_bm_rich_den <- car::Anova(gls_abg_bm_rich_den, type = "II")
emmip(lm_abg_bm_rich_den, Richness_ID ~ Density_ID) 

emm_abg_bm_dens <- emmeans(gls_abg_bm_rich_den, specs = "Density_ID")
pairs_abg_bm_dens <- pairs(emm_abg_bm_dens)
cld_abg_bm_dens <- cld(emm_abg_bm_dens)

emm_abg_bm_rich <- emmeans(gls_abg_bm_rich_den, specs = "Richness_ID")
pairs_abg_bm_rich <- pairs(emm_abg_bm_rich)
cld_abg_bm_rich <- cld(emm_abg_bm_rich)

## belowground biomass ---------------------------------------------------------
lm_bg_bm_rich_den <- lm(biomass ~ Richness_ID*Density_ID, data = bg_bm_tidy)
plot(lm_bg_bm_rich_den)

aov_bg_bm_rich_den <- car::Anova(lm_bg_bm_rich_den, type = "II")
emm_bg_bm_dens <- emmeans(lm_bg_bm_rich_den, specs = "Density_ID")
cld_bg_bm_dens <- cld(emm_bg_bm_dens)

# effect sizes -----------------------------------------------------------------

# for writing the results section
calc_percent_change <- function(V1, V2) {
  
  perc_diff <- V2 - V1
  abs_V2 <- abs(V1)
  perc_change <- (perc_diff/abs_V2)*100
  return(perc_change)
  
}

abg_bm_dens_D1 <- cld_abg_bm_dens[3,2]
abg_bm_dens_D2 <- cld_abg_bm_dens[2,2]
abg_bm_dens_D3 <- cld_abg_bm_dens[1,2]

abg_bm_rich_M1 <- cld_abg_bm_rich[3,2]
abg_bm_rich_D2 <- cld_abg_bm_rich[2,2]
abg_bm_rich_M4 <- cld_abg_bm_rich[1,2]

calc_percent_change(abg_bm_dens_D1, abg_bm_dens_D2)
calc_percent_change(abg_bm_dens_D2, abg_bm_dens_D3)
calc_percent_change(abg_bm_rich_M1, abg_bm_rich_M4)

# visualize --------------------------------------------------------------------

## species identity within monocultures ----------------------------------------

density_names <- as_labeller(
  c(
  `D1` = "48 native seeds", 
  `D2` = "200 native seeds",
  `D3` = "400 native seeds"
  )
)

(mono_spp_id <- cld_id_by_dens %>%
  as.data.frame() %>%
  mutate(
    
    # help make the error bars in the figure
    upper_se = emmean + SE,
    lower_se = emmean - SE,
    
    # create more appealing compact letter displays
    cld_abc = case_when(
      .group == " 1 "  ~ "A",
      .group == " 12"  ~ "AB",
      .group == "  2"  ~ "B",
      .group == " 1  " ~ "A",
      .group == " 12 " ~ "AB",
      .group == "  23" ~ "BC", 
      .group == "   3" ~ "C",
      TRUE ~ .group
    ),
    
  spp_identity = factor(
    spp_identity, 
    levels = c("A", "H", "M", "R", "O")
    )
  ) %>%
  ggplot(aes(x = spp_identity, y = emmean)) + 
  labs(
    x = "Species Identity (as monocultures)",
    y = "Aboveground biomass of Cirsium arvense (in grams)") + 
   geom_col() + 
   geom_errorbar(aes(ymax = upper_se, ymin = lower_se), width =.3) + 
   scale_x_discrete(
     labels = c(
       "A. gerardii", 
       "H. helianthoides", 
       "M. fistulosa", 
       "R. hirta",
       "O. biennis")
   ) + 
   geom_text(aes(label = cld_abc), nudge_x = 0.31, nudge_y = 0.22) + 
   facet_wrap(~Density_ID, labeller = density_names) + 
   theme_bw() +
   theme(axis.text.x = element_text(angle = 45, hjust = 1))
)

## aboveground biomass ---------------------------------------------------------

(plot_ab_bm_by_dens <- cld_abg_bm_dens %>%
  mutate(
    upper_se = emmean + SE, 
    lower_se = emmean - SE,
    cld_letters = case_when(
      .group == " 1 " ~ "A", 
      .group == "  2" ~ "B", 
      TRUE ~ .group)
  ) %>%
  ggplot(aes(x = Density_ID, y = emmean)) + 
  geom_col() + 
  geom_errorbar(aes(ymax = upper_se, ymin = lower_se), width =.3) + 
  geom_text(aes(label = cld_letters), nudge_x = 0.25, nudge_y = 0.15) + 
  scale_x_discrete(
    labels = c("48 native seeds", "200 native seeds", "400 native seeds")
  ) +
  ylim(0, 2) + 
  labs(
    x = "Native Seeding Density",
    y = NULL) + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
)

(plot_ab_bm_by_rich <- cld_abg_bm_rich %>%
    mutate(
      upper_se = emmean + SE, 
      lower_se = emmean - SE,
      cld_letters = case_when(
        .group == " 1 " ~ "A", 
        .group == " 12" ~ "AB",
        .group == "  2" ~ "B", 
        TRUE ~ .group)
    ) %>%
    ggplot(aes(x = Richness_ID, y = emmean)) + 
    geom_col() + 
    geom_errorbar(aes(ymax = upper_se, ymin = lower_se), width =.3) + 
    geom_text(aes(label = cld_letters), nudge_x = 0.15, nudge_y = 0.08) + 
    scale_x_discrete(
      labels = c("Monocultures", "2-species Mixture", "4-species Mixture")
    ) + 
    ylim(0, 2) + 
    labs(
      x = "Native Seeding Richness",
      y = "Aboveground biomass of Cirsium arvense (in grams)") + 
    theme_bw() + 
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
)
 
## belowground biomass ---------------------------------------------------------
(plot_bg_bm <- cld_bg_bm_dens %>%
    mutate(
      upper_se = emmean + SE, 
      lower_se = emmean - SE,
      cld_letters = case_when(
        .group == " 1 " ~ "B", 
        .group == "  2" ~ "A",
        TRUE ~ .group)
    ) %>%
    ggplot(aes(x = Density_ID, y = emmean)) + 
    geom_col() + 
    geom_errorbar(aes(ymax = upper_se, ymin = lower_se), width =.3) + 
    annotate("text", label = "B)", x = 0.6, y = 2.5) + 
    annotate(
     "text", 
     label = paste0("Seeding Density: ", ifelse(aov_bg_bm_rich_den$`Pr(>F)`[2] < 0.001, "p<0.001")),
     x = 2.7, y = 2.5,
     size = 3
    ) +
    annotate(
     "text", 
     label = paste0("Seeding Richness: p = ", round(aov_bg_bm_rich_den$`Pr(>F)`[1], digits = 2)),
     x = 2.7, y = 2.4,
     size = 3
    ) + 
    annotate(
     "text", 
     label = paste0("Density*Richness: p = ", round(aov_bg_bm_rich_den$`Pr(>F)`[3], digits = 2)),
     x = 2.7, y = 2.3,
     size = 3
    ) + 
    geom_text(aes(label = cld_letters), nudge_y = 0.1, nudge_x = 0.1) + 
    ylim(0, 2.5) + 
    scale_x_discrete(
      labels = c("48 Native Seeds", "200 Native Seeds", "400 Native Seeds")
    ) + 
    labs(
      x = "Native Seeding Density",
      y = "Belowground Invader Biomass (in grams)") + 
    theme_bw() + 
    theme(axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)))
)

## multipanel figures -----------------------------------------------------------

(plot_ab_bm <- plot_ab_bm_by_rich + plot_ab_bm_by_dens)

# save to disk -----------------------------------------------------------------
ggsave(
  plot = plot_ab_bm,
  filename = here("output", "results", "plot_ab_bm.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 8
)

ggsave(
  plot = mono_spp_id,
  filename = here("output", "results", "mono_ab_inv_bm.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 8
)
