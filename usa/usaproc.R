# Process USA data

library(tidyverse)
library(here)

# Read in case data

usc = read_csv(here("data/cases_by_county.csv"))
usc %>% rename(cases = rolling_average_cases_per_100k_centered) -> usc

# Read in WW data

usww = read_csv(here("data/wastewater_by_county.csv"))
usww %>% rename(date = sampling_week, N1 = normalized_concentration_rolling_average) -> usww

# Join at each week

left_join(usww,usc,by=c("date","fipscode","name")) %>% 
  select(date,N1,cases,name) %>% 
  rename(county = name) -> usa

save(usa,file = here("data/usa.rda"))  
