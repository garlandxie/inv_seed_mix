# libraries --------------------------------------------------------------------
library(here)
library(dplyr)
library(car)
library(ggplot2)
library(emmeans)
library(mgcv)
library(piecewiseSEM)

# import -----------------------------------------------------------------------

inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))
rgr_height <- read.csv(here("data", "intermediate_data", "rgr_height.csv"))
cum_germ_perc <- read.csv(here("data", "intermediate_data", "c_germ_week.csv"))
res_bm <- read.csv(here("data", "intermediate_data", "resident_biomass.csv"))
cwm_res <- read.csv(here("data", "intermediate_data", "cwm_res.csv"))

# data cleaning ----------------------------------------------------------------

rgr_height_res <- rgr_height %>%
  dplyr::filter(Spp != "CIAR") %>%
  group_by(Richness_ID, Density_ID, Rep) %>%
  summarize(
    mean_rgr_height_res = mean(rgr_height, na.rm = TRUE),
    mean_rgr_height_res = round(mean_rgr_height_res, digits = 2)
    ) %>%
  ungroup() %>%
  janitor::clean_names()

rgr_height_ciar <- rgr_height %>%
  dplyr::filter(Spp == "CIAR") %>%
  rename(mean_rgr_height_ciar = rgr_height) %>%
  janitor::clean_names()

# perform inner joins ----------------------------------------------------------

multi_key_id <- c("richness_id", "density_id", "rep")

sem_df <- inv_bm %>%
  janitor::clean_names() %>%
  full_join(rgr_height_res, by = multi_key_id) %>%
  full_join(rgr_height_ciar, by = multi_key_id) %>%
  full_join(cum_germ_perc, by = multi_key_id) %>%
  full_join(janitor::clean_names(res_bm), by = multi_key_id) %>%
  full_join(janitor::clean_names(cwm_res), by = multi_key_id) %>%
  dplyr::select(
    richness_id, 
    density_id, 
    rep, 
    tot_inv_bm_g,
    num_invaders, 
    res_comm_biomass_g, 
    mean_rgr_height_res,
    mean_rgr_height_ciar, 
    sown_seeds_res = sown_seeds, 
    cum_germ_perc_ciar, 
    cum_germ_perc_res,
    mean_gsp, 
    wds_sla
  ) 
    
# statistical analyses ---------------------------------------------------------

## % germination (resident) <- richness * density + germination speed ----------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = richness_id, y = cum_germ_perc_res, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

sem_df %>%
  ggplot(aes(x = mean_gsp, y = cum_germ_perc_res)) + 
  geom_point() + 
  labs(
    x = "Germination Speed", 
    y = "% Germination") +
  theme_bw() 

# model fit
glm_res_germ <- glm(
  cum_germ_perc_res ~ richness_id + density_id + mean_gsp, 
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df)

# sanity checks
plot(glm_res_germ, 1)
plot(glm_res_germ, 2)
plot(glm_res_germ, 3)
plot(glm_res_germ, 4)

summary(glm_res_germ)

## germination speed (resident) <- richness * density --------------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = richness_id, y = mean_gsp, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Germination Speed") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

# model fit
glm_res_gsp <- glm(
  mean_gsp ~ density_id + richness_id, 
  family = Gamma(link = "log"), 
  data = sem_df)

# test global significance using one-way ANOVA
car::Anova(glm_res_gsp , type = "II")

## RGR height (resident) <- richness * density ---------------------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = richness_id, y = mean_rgr_height_res, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

# model fit 
lm_res_height <- lm(
  mean_rgr_height_res ~ richness_id + density_id, 
  data = sem_df)

# test global significance using one-way ANOVA
car::Anova(lm_res_height, type = "II")
summary(lm_res_height) 

# sanity checks
plot(lm_res_height, 1)
plot(lm_res_height, 2)
plot(lm_res_height, 3)
plot(lm_res_height, 4)

# check pairwise comparisons
pairs_height_rich_dens <- lm_res_height %>%
  emmeans::emmeans("density_id") %>% 
  pairs()

## biomass (resident) <- % germination (resident) + height (resident) ----------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = res_comm_biomass_g, y = cum_germ_perc_res)) + 
  geom_point() + 
  labs(
    x = "Percent Germination of Resident Community", 
    y = "Resident Community Biomass (in g)") + 
  theme_bw() 

sem_df %>%
  ggplot(aes(x = mean_rgr_height_res, y = res_comm_biomass_g)) + 
  geom_point() + 
  labs(
    x = "Relative Growth Rate (Height)", 
    y = "Resident Community Biomass (in mg)") + 
  theme_bw() 

# model fit
lm_res_bm <- lm(
  res_comm_biomass_g ~ cum_germ_perc_res + mean_rgr_height_res, 
  data = sem_df)

summary(lm_bm_res)

# sanity checks
plot(lm_res_bm, 1)
plot(lm_res_bm, 2)
plot(lm_res_bm, 3)
plot(lm_res_bm, 4)

## % germination (invader) <- community biomass (resident) ---------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = res_comm_biomass_mg, y = cum_germ_perc_ciar)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Resident Community Biomass")
  theme_bw() 

# model fit
glm_inv_germ <- glm(
  cum_germ_perc_ciar ~ res_comm_biomass_g, 
  family = binomial(link = "logit"),
  data = sem_df
  )

summary(glm_inv_germ)

# sanity checks
plot(glm_inv_germ, 1)
plot(glm_inv_germ, 2)
plot(glm_inv_germ, 3)
plot(glm_inv_germ, 4)

## RGR height (invader) <- community biomass (resident) ------------------------

# visualize data before running models 
lm_inv_height <- lm(
  mean_rgr_height_ciar ~ res_comm_biomass_g, 
  data = sem_df
)

summary(lm_inv_height)

# sanity checks
plot(lm_inv_height, 1)
plot(lm_inv_height, 2)
plot(lm_inv_height, 3)
plot(lm_inv_height, 4)

## biomass (invader) <- % germination (invader) + height (invader) + SLA -------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = cum_germ_perc_ciar, y =  tot_inv_bm_g)) + 
  geom_point() + 
  labs(x = "Percentage Germination (Invader)",
       y = "Total Invader Biomass") + 
  theme_bw() 

sem_df %>%
  ggplot(aes(x = mean_rgr_height_ciar, y =  tot_inv_bm_g)) + 
  geom_point() + 
  theme_bw() 

sem_df %>% 
  ggplot(aes(x = wds_sla, y = tot_inv_bm_g)) + 
  geom_point()+ 
  geom_smooth(method = "lm") + 
  labs(x = "Community weighted mean dissimilarity", 
       y = "Invader biomass (g)"
  ) + 
  theme_bw() 

# model fit
lm_inv_bm <- lm(
  log(tot_inv_bm_g) ~ cum_germ_perc_ciar + mean_rgr_height_ciar + wds_sla, 
  data = sem_df
)

summary(inv_bm)

# sanity checks
plot(lm_inv_bm, 1)
plot(lm_inv_bm, 2)
plot(lm_inv_bm, 3)
plot(lm_inv_bm, 4)

## CWM SLA (resident) <- community biomass (resident) --------------------------

# visualize data before running models
sem_df %>% 
  ggplot(aes(x = res_comm_biomass_g, y = wds_sla)) + 
  geom_point() + 
  labs(x = "Resident Community Biomass (mg)",
       y = "Mean Weighted Community Dissimilarity") + 
  theme_bw() 

# model fit
glm_res_sla <- glm(
  wds_sla ~ res_comm_biomass_g,
  family = Gamma(link = "log"),
  data = sem_df
)

summary(glm_res_sla)

# sanity checks
plot(glm_res_sla, 1)
plot(glm_res_sla, 2)
plot(glm_res_sla, 3)
plot(glm_res_sla, 4)

# run piecewise SEM ------------------------------------------------------------

sem_germ <- piecewiseSEM::psem(
  glm_res_germ,
  glm_res_gsp,
  lm_res_height,
  lm_res_bm,
  glm_res_sla,
  glm_inv_germ,
  lm_inv_height,
  lm_inv_bm
)

summary(sem_germ, conserve = TRUE)
