# libraries --------------------------------------------------------------------
library(here)
library(dplyr)
library(car)
library(ggplot2)
library(emmeans)

# import -----------------------------------------------------------------------

inv_bm <- read.csv(here("data", "intermediate_data", "invader_biomass.csv"))
rgr_height <- read.csv(here("data", "intermediate_data", "rgr_height.csv"))
cum_germ_perc <- read.csv(here("data", "intermediate_data", "c_germ_week12.csv"))
res_bm <- read.csv(here("data", "intermediate_data", "resident_biomass.csv"))

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
  dplyr::select(
    richness_id, 
    density_id, 
    rep, 
    tot_inv_bm_mg,
    num_invaders, 
    res_comm_biomass_mg, 
    mean_rgr_height_res,
    mean_rgr_height_ciar, 
    sown_seeds_res = sown_seeds, 
    cum_germ_perc_ciar, 
    cum_germ_perc_res
  ) 
    
# statistical analyses ---------------------------------------------------------

## % germination (resident) <- richness * density ------------------------------

sem_df %>%
  ggplot(aes(x = richness_id, y = cum_germ_perc_res, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Percentage Germination (Resident Community)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

lm_res_germ_rich_dens1 <- glm(
  cum_germ_perc_res ~ richness_id*density_id, 
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df)

lm_res_germ_rich_dens2 <- glm(
  cum_germ_perc_res ~ richness_id + density_id, 
  family = binomial(link = "logit"), 
  weights = sown_seeds_res, 
  data = sem_df)

AIC(lm_res_germ_rich_dens2, lm_res_germ_rich_dens1)

car::Anova(lm_res_germ_rich_dens2, type = "II")
summary(lm_res_germ_rich_dens2)
plot(lm_res_germ_rich_dens2)

pairs_res_germ_rich_dens2 <- lm_res_germ_rich_dens2 %>%
  emmeans::emmeans("richness_id") %>% 
  pairs()

## height <- richness * density ------------------------------------------------

sem_df %>%
  ggplot(aes(x = richness_id, y = mean_rgr_height_res, fill = density_id)) + 
  geom_boxplot() + 
  labs(
    x = "Seeding Richness", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

lm_height_rich_dens1 <- lm(
  mean_rgr_height_res ~ richness_id*density_id, 
  data = sem_df)

lm_height_rich_dens2 <- lm(
  mean_rgr_height_res ~ richness_id + density_id, 
  data = sem_df)

AIC(lm_height_rich_dens2, lm_height_rich_dens1)
car::Anova(lm_height_rich_dens2, type = "II")
summary(lm_height_rich_dens2) 
plot(lm_height_rich_dens2)

## height <- % germination (resident) ------------------------------------------

sem_df %>%
  ggplot(aes(x = cum_germ_perc_res, y = mean_rgr_height_res)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  labs(
    x = "Percent Germination of Resident Community", 
    y = "Relative Growth Rate (Height)") +
  scale_fill_discrete(name = "Seeding Density") + 
  theme_bw() 

lm_height_germ_res <- lm(
  mean_rgr_height_res ~ cum_germ_perc_res, 
  data = sem_df)

summary(lm_height_germ_res)
plot(lm_height_germ_res)

## resident biomass <- % germination (resident) --------------------------------

sem_df %>%
  ggplot(aes(x = cum_germ_perc_res, y = res_comm_biomass_mg)) + 
  geom_point() + 
  labs(
    x = "Percent Germination of Resident Community", 
    y = "Resident Community Biomass (in mg)") + 
  theme_bw() 

lm_bm_germ_res <- lm(
  res_comm_biomass_mg ~ cum_germ_perc_res, 
  data = sem_df)

summary(lm_bm_germ_res)
plot(lm_bm_germ_res)

## resident biomass <- height (resident) ---------------------------------------

sem_df %>%
  ggplot(aes(x = mean_rgr_height_res, y = res_comm_biomass_mg)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  labs(
    x = "Relative Growth Rate (Height)", 
    y = "Resident Community Biomass (in mg)") + 
  theme_bw() 

lm_bm_height_res <- lm(
  res_comm_biomass_mg ~ mean_rgr_height_res, 
  data = sem_df)

summary(lm_bm_height_res)

## % germination (invader) <- resident biomass ---------------------------------

sem_df %>%
  ggplot(aes(x = res_comm_biomass_mg, y = cum_germ_perc_ciar)) + 
  geom_point() + 
  geom_smooth(method = "lm") + 
  theme_bw() 

lm_res_bm_germ_inv <- lm(
  cum_germ_perc_ciar ~ res_comm_biomass_mg, 
  data = sem_df
  )

summary(lm_res_bm_germ_inv)
plot(lm_res_bm_germ_inv)


## % germination (invader) <- height (resident) --------------------------------

sem_df %>%
  ggplot(aes(x = res_comm_biomass_mg, y =  mean_rgr_height_ciar)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Resident Community Biomass",
       y = "Relative Growth Rate (Height of Invader)") + 
  theme_bw() 

lm_res_bm_height_inv <- lm(
  mean_rgr_height_ciar ~ res_comm_biomass_mg, 
  data = sem_df
)

summary(lm_res_bm_height_inv)
plot(lm_res_bm_germ_inv)

## biomass (invader) <- % germination (invader) --------------------------------

sem_df %>%
  ggplot(aes(x = cum_germ_perc_ciar, y =  tot_inv_bm_mg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  labs(x = "Percentage Germination (Invader)",
       y = "Total Invader Biomass") + 
  theme_bw() 

lm_bm_inv_germ_inv <- lm(
  tot_inv_bm_mg ~ cum_germ_perc_ciar, 
  data = sem_df
)

summary(lm_bm_inv_germ_inv)
plot(lm_bm_inv_germ_inv)

## biomass (invader)  <- height (invader) --------------------------------------

sem_df %>%
  ggplot(aes(x = mean_rgr_height_ciar, y =  tot_inv_bm_mg)) + 
  geom_point() + 
  geom_smooth(method = "lm", se = FALSE) + 
  theme_bw() 

lm_bm_inv_germ_inv <- lm(
  tot_inv_bm_mg ~ mean_rgr_height_ciar, 
  data = sem_df
)

summary(lm_bm_inv_germ_inv)
plot(lm_bm_inv_germ_inv)
