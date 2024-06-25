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

## % germination (resident) <- richness * density ------------------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = richness_id, y = cum_germ_perc_res, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

# model fit
glm_res_germ_rich_dens <- glm(
  cum_germ_perc_res ~ richness_id*density_id, 
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df)

# test global significance using one-way ANOVA
car::Anova(glm_res_germ_rich_dens, type = "II")
summary(glm_res_germ_rich_dens2)

# sanity checks
plot(glm_res_germ_rich_dens2, 1)
plot(glm_res_germ_rich_dens2, 2)
plot(glm_res_germ_rich_dens2, 3)
plot(glm_res_germ_rich_dens2, 4)

# check pairwise comparisons
pairs_res_germ_rich_dens2 <- glm_res_germ_rich_dens2 %>%
  emmeans::emmeans("richness_id") %>% 
  pairs()

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
glm_gsp_rich_dens <- glm(
  mean_gsp ~ density_id*richness_id, 
  family = Gamma(link = "log"), 
  data = sem_df)

# test global significance using one-way ANOVA
car::Anova(glm_gsp_rich_dens , type = "II")

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
lm_height_rich_dens <- lm(
  mean_rgr_height_res ~ richness_id*density_id, 
  data = sem_df)

# test global significance using one-way ANOVA
car::Anova(lm_height_rich_dens, type = "II")
summary(lm_height_rich_dens) 

# sanity checks
plot(lm_height_rich_dens, 1)
plot(lm_height_rich_dens, 2)
plot(lm_height_rich_dens, 3)
plot(lm_height_rich_dens, 4)

# check pairwise comparisons
pairs_height_rich_dens <- lm_height_rich_dens %>%
  emmeans::emmeans("density_id") %>% 
  pairs()


## % germination (resident) <- germination speed (resident) --------------------

# visualize data before running models 
sem_df %>%
  ggplot(aes(x = mean_gsp, y = cum_germ_perc_res)) + 
  geom_point() + 
  labs(
    x = "Germination Speed", 
    y = "% Germination") +
  theme_bw() 

# model fit 
glm_res_germ_gsp <- glm(cum_germ_perc_res ~ mean_gsp, 
    family = binomial(link = "logit"), 
    weights = sown_seeds_res, 
    data = sem_df
    )

summary(glm_res_germ_gsp)

# sanity check 
plot(glm_res_germ_gsp, 1)
plot(glm_res_germ_gsp, 2)
plot(glm_res_germ_gsp, 3)
plot(glm_res_germ_gsp, 4)

## RGR height (resident) <- % germination (resident) ---------------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = cum_germ_perc_res, y = mean_rgr_height_res)) + 
  geom_point() + 
  labs(
    x = "Percent Germination of Resident Community", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

# model fit
lm_height_germ_res <- lm(
  mean_rgr_height_res ~ cum_germ_perc_res, 
  data = sem_df)

summary(lm_height_germ_res)

# sanity checks
plot(lm_height_germ_res, 1)
plot(lm_height_germ_res, 2)
plot(lm_height_germ_res, 3)
plot(lm_height_germ_res, 4)

## community biomass (resident) <- % germination (resident) --------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = cum_germ_perc_res, y = res_comm_biomass_g)) + 
  geom_point() + 
  labs(
    x = "Percent Germination of Resident Community", 
    y = "Resident Community Biomass (in g)") + 
  theme_bw() 

# model fit
lm_bm_germ_res <- lm(
  res_comm_biomass_g ~ cum_germ_perc_res, 
  data = sem_df)

summary(lm_bm_germ_res)

# sanity checks
plot(lm_bm_germ_res, 1)
plot(lm_bm_germ_res, 2)
plot(lm_bm_germ_res, 3)
plot(lm_bm_germ_res, 4)

## community biomass (resident) <- RGR height (resident) -----------------------

# visualizing data before running models
sem_df %>%
  ggplot(aes(x = mean_rgr_height_res, y = res_comm_biomass_g)) + 
  geom_point() + 
  labs(
    x = "Relative Growth Rate (Height)", 
    y = "Resident Community Biomass (in mg)") + 
  theme_bw() 

# model fit
lm_bm_height_res <- lm(
  res_comm_biomass_g ~ mean_rgr_height_res,
  data = sem_df)

summary(lm_bm_height_res)

# sanity check
plot(lm_bm_height_res, 1)
plot(lm_bm_height_res, 2)
plot(lm_bm_height_res, 3)
plot(lm_bm_height_res, 4)

## % germination (invader) <- community biomass (resident) ---------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = res_comm_biomass_mg, y = cum_germ_perc_ciar)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(x = "Resident Community Biomass")
  theme_bw() 

# model fit
glm_res_bm_germ_inv <- glm(
  cum_germ_perc_ciar ~ res_comm_biomass_mg, 
  family = binomial(link = "logit"),
  data = sem_df
  )

summary(glm_res_bm_germ_inv)

# sanity checks
plot(glm_res_bm_germ_inv, 1)
plot(glm_res_bm_germ_inv, 2)
plot(glm_res_bm_germ_inv, 3)
plot(glm_res_bm_germ_inv, 4)

## community biomass (resident) <- RGR height (invader) ------------------------

# visualize data before running models 
sem_df %>%
  ggplot(aes(x = res_comm_biomass_mg, y =  mean_rgr_height_ciar)) + 
  geom_point() + 
  labs(x = "Resident Community Biomass",
       y = "Relative Growth Rate (Height of Invader)") + 
  theme_bw() 

lm_res_bm_height_inv <- lm(
  mean_rgr_height_ciar ~ res_comm_biomass_mg, 
  data = sem_df
)

summary(lm_res_bm_height_inv)

# sanity checks
plot(lm_res_bm_height_inv, 1)
plot(lm_res_bm_height_inv, 2)
plot(lm_res_bm_height_inv, 3)
plot(lm_res_bm_height_inv, 4)

## biomass (invader) <- % germination (invader) --------------------------------

# visualize data before running models
sem_df %>%
  ggplot(aes(x = cum_germ_perc_ciar, y =  tot_inv_bm_g)) + 
  geom_point() + 
  labs(x = "Percentage Germination (Invader)",
       y = "Total Invader Biomass") + 
  theme_bw() 

# model fit
lm_bm_inv_germ_inv <- lm(
  tot_inv_bm_g ~ cum_germ_perc_ciar, 
  data = sem_df
)

summary(lm_bm_inv_germ_inv)

# sanity checks
plot(lm_bm_inv_germ_inv, 1)
plot(lm_bm_inv_germ_inv, 2)
plot(lm_bm_inv_germ_inv, 3)
plot(lm_bm_inv_germ_inv, 4)

## biomass (invader) <- RGR height (invader) -----------------------------------

# visualizing data befoer running models 
sem_df %>%
  ggplot(aes(x = mean_rgr_height_ciar, y =  tot_inv_bm_g)) + 
  geom_point() + 
  theme_bw() 

# model fit
lm_bm_inv_height_inv <- lm(
  tot_inv_bm_g ~ mean_rgr_height_ciar, 
  data = sem_df
)

summary(lm_bm_inv_height_inv)

# sanity checks
plot(lm_bm_inv_height_inv, 1)
plot(lm_bm_inv_height_inv, 2)
plot(lm_bm_inv_height_inv, 3)
plot(lm_bm_inv_height_inv, 4)

## biomass (invader) <- % germination (invader) --------------------------

# visualize data before running the models
sem_df %>%
  ggplot(aes(x = cum_germ_perc_ciar, y =  tot_inv_bm_g)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() 

# model fit
lm_bm_inv_germ_inv <- lm(
  tot_inv_bm_g ~ cum_germ_perc_ciar, 
  data = sem_df
)

summary(lm_bm_inv_germ_inv)

# sanity checks
plot(lm_bm_inv_germ_inv, 1)
plot(lm_bm_inv_germ_inv, 2)
plot(lm_bm_inv_germ_inv, 3)
plot(lm_bm_inv_germ_inv, 4)

## biomass (invader) <- CWM SLA (resident) ---------------------------------------

# visualize data before running the models
sem_df %>% 
  ggplot(aes(x = wds_sla, y = tot_inv_bm_g)) + 
  geom_point()+ 
  geom_smooth(method = "lm") + 
  labs(x = "Community weighted mean dissimilarity", 
       y = "Invader biomass (g)"
  ) + 
  theme_bw() 

# model fit
lm_sla_inv_bm <- lm(
  tot_inv_bm_g ~ wds_sla, 
  data = sem_df
)

summary(lm_sla_inv_bm)

# sanity checks
plot(lm_sla_inv_bm, 1)
plot(lm_sla_inv_bm, 2)
plot(lm_sla_inv_bm, 3)
plot(lm_sla_inv_bm, 4)

## CWM SLA (resident) <- community biomass (resident) --------------------------

# a potential non-linear relationship seen in the raw data
sem_df %>% 
  ggplot(aes(x = res_comm_biomass_g, y = wds_sla)) + 
  geom_point() + 
  labs(x = "Resident Community Biomass (mg)",
       y = "Mean Weighted Community Dissimilarity") + 
  theme_bw() 

# fit using ordinary least squares 
lm_sla_res_bm <- lm(
  wds_sla ~ res_comm_biomass_g,
  data = sem_df
)

# fit using generalized additive model
gam_sla_res_bm <- mgcv::gam(
  wds_sla ~ s(res_comm_biomass_g), 
  family = Gamma(link = "log"),
  method = "ML",
  data = sem_df
)

# fit using a generalized linear model 
glm_sla_res_bm <- glm(
  wds_sla ~ res_comm_biomass_g,
  family = Gamma(link = "log"),
  data = sem_df
)

# compare these two models
# values are very similar between GLM and GAM, so just choose the simpler model
AIC(gam_sla_res_bm, glm_sla_res_bm, lm_sla_res_bm)

summary(glm_sla_res_bm)

# sanity checks
plot(glm_sla_res_bm, 1)
plot(glm_sla_res_bm, 2)
plot(glm_sla_res_bm, 3)
plot(glm_sla_res_bm, 4)

# run piecewise SEM ------------------------------------------------------------

piecewiseSEM::psem(
  
)