# libraries -----
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data 
library(ggplot2)    # for visualizing data  
library(patchwork)  # for making multi-panel figures 

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

## cumulative germination (as a percentage) -----

c_germ_perc_tidy2 <- c_germ_perc_tidy %>%
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
  ) 

## emergence time -----



# data visualization: cumulative percentage germination ------------------------

## M1 D1 -----------------------------------------------------------------------

(plot_m1_d1 <- c_germ_perc_tidy2 %>%
  dplyr::filter(richness_id == "M1" & density_id == "D1") %>% 
  ggplot() + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
  geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
  geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
  scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
  labs(x = NULL, y = NULL) + 
  ylim(0, 50) + 
  theme_bw() +
  theme(axis.text.x = element_blank())
)

## M1 D2 -----------------------------------------------------------------------

(plot_m1_d2 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M1" & density_id == "D2") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() + 
   theme(axis.text.x = element_blank())
)

## M1 D3 -----------------------------------------------------------------------

(plot_m1_d3 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M1" & density_id == "D3") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   labs(x = NULL, y = NULL) + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   ylim(0, 50) + 
   theme_bw() +
   theme(axis.text.x = element_blank())
)

## M2 D1 -----------------------------------------------------------------------

(plot_m2_d1 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M2" & density_id == "D1") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") +
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() + 
   theme(
     axis.text.x = element_blank(),
     axis.title.y = element_text(
       size = 20,
       margin = margin(t = 0, r = 10, b = 0, l = ))
     )
)

## M2 D2 -----------------------------------------------------------------------

(plot_m2_d2 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M2" & density_id == "D2") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() +
   theme(axis.text.x = element_blank())
)

## M2 D3 -----------------------------------------------------------------------

(plot_m2_d3 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M2" & density_id == "D3") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() + 
   theme(axis.text.x = element_blank())
)

## M4 D1 -----------------------------------------------------------------------

(plot_m4_d1 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M4" & density_id == "D1") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   labs(x = NULL, y = NULL) + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   ylim(0, 50) + 
   theme_bw() 
)

## M4 D2 -----------------------------------------------------------------------

(plot_m4_d2 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M4" & density_id == "D2") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() + 
   theme(axis.title.x = element_text(
     size = 20, 
     margin = margin(t = 10, r = 0, b = 0, l = 0)))
)

## M4 D3 -----------------------------------------------------------------------

(plot_m4_d3 <- c_germ_perc_tidy2 %>%
   dplyr::filter(richness_id == "M4" & density_id == "D3") %>% 
   ggplot() + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_oebi), col = "orange") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ruhi), col = "red") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_ange), col = "blue") + 
   geom_line(aes(x = week, y = mean_cum_perc_germ_mofi), col = "green") +
   geom_line(aes(x = week, y = mean_cum_perc_germ_hehe), col = "black") + 
   scale_x_continuous(breaks = c(2,4,6,8,10,12)) + 
   labs(x = NULL, y = NULL) + 
   ylim(0, 50) + 
   theme_bw() 
)

## multipanel figure -----------------------------------------------------------

(plot_cum_germ <- (plot_m1_d1 + plot_m1_d2 + plot_m1_d3) /
                  (plot_m2_d1 + plot_m2_d2 + plot_m2_d3) /
                  (plot_m4_d1 + plot_m4_d2 + plot_m4_d3)
)

# save to disk -----------------------------------------------------------------

ggsave(
  plot = plot_cum_germ, 
  filename = here("output", "results", "germ_curves.png"), 
  device = "png", 
  units = "in", 
  height = 6, 
  width = 8
)
