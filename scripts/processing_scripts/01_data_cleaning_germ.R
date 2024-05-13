# libraries -----
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data 
library(ggplot2)    # for visualizing data  

# import ----
c_germ <- read.csv(here("data", "input_data", "cumulative_germination.csv"))

# data cleaning ----

## cumulative germination ----
c_germ_tidy <- c_germ %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ, names_from = spp) %>%
  ungroup() %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep, 
    
    cum_germ_ciar = CIAR, 
    cum_germ_oebi = OEBI, 
    cum_germ_ruhi = RUHI, 
    cum_germ_ange = ANGE, 
    cum_germ_hehe = HEHE, 
    cum_germ_mofi = MOFI) %>%
  mutate(across(cum_germ_ciar:cum_germ_mofi, ~replace(., is.na(.), 0))) %>%
  mutate(
    cum_germ_res = 
      cum_germ_oebi + 
      cum_germ_ruhi + 
      cum_germ_ange + 
      cum_germ_hehe + 
      cum_germ_mofi
    ) 

## cumulative germination (as a percentage) ----

# invader-only
c_germ_perc_ciar <- c_germ %>%
  filter(spp == "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ_perc, names_from = spp) %>%
  ungroup() %>%
  rename(cum_germ_perc_ciar = CIAR)

c_germ_perc_res <- c_germ %>% 
  filter(spp != "CIAR") %>%
  group_by(week, richness_id, density_id, rep) %>%
  summarize(cum_germ_res = sum(cum_germ, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(
    sown_seeds = case_when(
      density_id == "D1" ~ 48,
      density_id == "D2" ~ 200, 
      density_id == "D3" ~ 400
      ), 
    cum_germ_perc_res = cum_germ_res/sown_seeds,
    cum_germ_perc_res = round(cum_germ_perc_res, digits = 2)
  ) 
  
# percentage for all native species within a tray 
c_germ_perc <- c_germ_perc_ciar %>% 
  inner_join(
    c_germ_perc_res, 
    by = c("week", "richness_id", "density_id", "rep")
    ) %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep,
    sown_seeds = sown_seeds.y, 
    cum_germ_perc_ciar, 
    cum_germ_perc_res
  )

# for data visualizing
c_germ_perc_tidy <- c_germ %>%
  group_by(week, richness_id, density_id, rep) %>%
  tidyr::pivot_wider(values_from = cum_germ, names_from = spp) %>%
  ungroup() %>%
  dplyr::select(
    week, 
    richness_id, 
    density_id, 
    rep, 
    
    cum_germ_perc_ciar = CIAR, 
    cum_germ_perc_oebi = OEBI, 
    cum_germ_perc_ruhi = RUHI, 
    cum_germ_perc_ange = ANGE, 
    cum_germ_perc_hehe = HEHE, 
    cum_germ_perc_mofi = MOFI) %>%
  
  mutate(across(
    cum_germ_perc_ciar:cum_germ_perc_mofi, 
    ~replace(., is.na(.), 0))
    ) %>%
  
  mutate(
    cum_germ_perc_res = 
      cum_germ_perc_oebi + 
      cum_germ_perc_ruhi + 
      cum_germ_perc_ange + 
      cum_germ_perc_hehe + 
      cum_germ_perc_mofi
  ) 

# data visualization ----

## cumulative germination ----
c_germ_tidy %>%
  group_by(week, richness_id, density_id) %>%
  summarize(
    mean_cum_germ_ciar = mean(cum_germ_ciar),
    mean_cum_germ_oebi = mean(cum_germ_oebi),
    mean_cum_germ_ruhi = mean(cum_germ_ruhi), 
    mean_cum_germ_ange = mean(cum_germ_ange), 
    mean_cum_germ_hehe = mean(cum_germ_hehe), 
    mean_cum_germ_mofi = mean(cum_germ_mofi)
  ) %>%
  ungroup() %>%
  dplyr::filter(
    richness_id %in% c("M1", "M2", "M4"), 
    density_id %in% c("D1", "D2", "D3")
    ) %>%
  ggplot() +
  geom_line(aes(x = week, y = mean_cum_germ_ciar), col = "purple") +
  geom_line(aes(x = week, y = mean_cum_germ_oebi), col = "orange") + 
  geom_line(aes(x = week, y = mean_cum_germ_ruhi), col = "red") + 
  geom_line(aes(x = week, y = mean_cum_germ_ange), col = "blue") + 
  geom_line(aes(x = week, y = mean_cum_germ_mofi), col = "green") + 
  geom_line(aes(x = week, y = mean_cum_germ_hehe), col = "grey") +
  labs(
    x = "Week of Observation", 
    y = "Cumulative germination") + 
  facet_wrap(richness_id~density_id) + 
  theme_bw()

## cumulative germination (as a percentage) -----

c_germ_perc_tidy %>%
  group_by(week, density_id, richness_id) %>%
  summarize(
    mean_cum_perc_germ_ciar = mean(cum_germ_perc_ciar),
    mean_cum_perc_germ_oebi = mean(cum_germ_perc_oebi),
    mean_cum_perc_germ_ruhi = mean(cum_germ_perc_ruhi), 
    mean_cum_perc_germ_ange = mean(cum_germ_perc_ange), 
    mean_cum_perc_germ_hehe = mean(cum_germ_perc_hehe), 
    mean_cum_perc_germ_mofi = mean(cum_germ_perc_mofi)
  ) %>%
  ungroup() %>%
  dplyr::filter(
    density_id %in% c("D1", "D2", "D3") & 
    richness_id %in% c("M1", "M2", "M4")
    ) %>%
  ggplot() +
  geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
  geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "grey") + 
  labs(x = "Week of Observation", y = "Cumulative germination") + 
  facet_wrap(density_id~richness_id) + 
  theme_bw() 

# some analyses ----
c_germ_w12 <- dplyr::filter(c_germ_perc, week == 12) 
glm_perc_ciar <- glm(
  cum_germ_perc_ciar ~ richness_id*density_id, 
  family = binomial(link = "logit"), 
  data = c_germ_w12
)

sum(residuals(glm_perc_ciar, "pearson")^2) / glm_perc_ciar$df.residual

glm_perc_ciar2 <- glm(
  cum_germ_perc_ciar ~ richness_id*density_id, 
  family = quasibinomial(link = "logit"), 
  data = c_germ_w12
)

car::Anova(glm_perc_ciar2, type = "II")

# save to disk -----

write.csv(
  c_germ_perc, 
  file = here("data", "intermediate_data", "c_germ_clean.csv")
)




