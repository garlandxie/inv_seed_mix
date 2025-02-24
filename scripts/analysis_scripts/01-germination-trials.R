# libraries ----
library(dplyr)
library(here)
library(ggplot2)
library(lubridate)

# import ----
germ_trials <- read.csv(here("data", "input_data", "germination_trials.csv"))

# visualize ----

(plot_germ <- germ_trials %>%
  mutate(
    Species = case_when(
      Species == "ASSY" ~ "A. syriaca", 
      Species == "HEHE" ~ "H. helianthoides",
      Species == "MOFI" ~ "M. fistulosa", 
      Species == "ANGE" ~ "A. gerardii",
      Species == "CIAR" ~ "C. arvense", 
      Species == "OEBI" ~ "O. biennis", 
      Species == "RUHI" ~ "R. hirta", 
      Species == "CONTROL" ~ "Control (No species)",
      TRUE ~ Species
    ) %>%
    factor(levels = c(
      "Control (No species)", 
      "A. syriaca",
      "M. fistulosa", 
      "H. helianthoides", 
      "A. gerardii", 
      "C. arvense",
      "R. hirta", 
      "O. biennis")
      ),
    Rep = factor(Rep),
    MonthDay = lubridate::date(Record_date_ymd) %>% format("%m-%d")
  ) %>%
  ggplot(aes(x = MonthDay, y = Germination, group = Rep)) + 
  geom_line((aes(col = Rep))) + 
  geom_point() + 
  ylim(0, 100) + 
  labs(
    x = "Record Date (YMD)",
    y = "Cumulative germination") + 
  facet_wrap(~Species, nrow = 2, ncol = 4) + 
  theme_bw() 
)

# save to disk ----
ggsave(
  plot = plot_germ, 
  filename = here("output", "data_appendix_output", "germ_trials.png"), 
  device = "png", 
  units = "in", 
  height = 5, 
  width = 10
)