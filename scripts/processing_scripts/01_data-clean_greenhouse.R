# libraries --------------------------------------------------------------------
library(here)       # for creating relative file-paths
library(dplyr)      # for manipulating data 
library(lubridate)  # for manipulating date-time formats

# import -----------------------------------------------------------------------

gr_df <- read.delim(here(
  "data",
  "input_data", 
  "greenhouse_conditions_March-1_to_June-30_2023.csv"
  ),
  sep = ";"
)

# data clean -------------------------------------------------------------------

gr_tidy <- gr_df %>%
  janitor::clean_names() %>%
  rename(record_date = x) %>%
  dplyr::filter(record_date != "") %>%
  mutate(
    record_date = lubridate::parse_date_time(record_date, "dmy HMS"),
    meas_gh_temp = as.numeric(meas_gh_temp),
    meas_gh_temp_1 = as.numeric(meas_gh_temp_1), 
    meas_rel_hum = as.numeric(meas_rel_hum), 
    meas_light = as.numeric(meas_light)
    ) %>%
  filter(record_date > lubridate::ymd("2023-03-28")) %>%
  filter(meas_gh_temp != 0.10) # remove outlier
  
# calculate summary statistics -------------------------------------------------

mean(gr_tidy$meas_gh_temp_1, na.rm = TRUE)
range(gr_tidy$meas_gh_temp_1, na.rm = TRUE)

mean(gr_tidy$meas_rel_hum, na.rm = TRUE)
range(gr_tidy$meas_rel_hum, na.rm = TRUE)

mean(gr_tidy$meas_light, na.rm = TRUE)
range(gr_tidy$meas_light, na.rm = TRUE)
