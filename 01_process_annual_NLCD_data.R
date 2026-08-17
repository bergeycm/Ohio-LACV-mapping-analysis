#!/usr/bin/env Rscript

# Title: process_annual_NLCD_data.R (v1.0.0)
# Authors:
#   Morgan E. Chaney (Rutgers University)
#   Christina M. Bergey (Rutgers University) <christina.bergey@rutgers.edu>
# Created: 2026 (revision reanalysis data preparation)
# License: MIT (see LICENSE file)
# Description: Builds the county-year analysis panel for the Ohio La Crosse virus / Amish
#   reanalysis. Reads annual National Land Cover Database (NLCD) county tabulations and
#   computes per-county land-cover proportions, then joins LACV-ND case counts, Amish
#   population share, and total/child population to produce the combined panel (df_final)
#   saved as ForPart2.RData and consumed by Amish_LACV_full_analysis.R.
# Usage:
#   Rscript process_annual_NLCD_data.R      # edit setwd()/paths below for your machine
# Inputs:
#   - Annual NLCD land-cover county tabulations (Excel/CSV)
#   - CDC LACV-ND case counts by county and year
#   - U.S. Census total and child (ages 0-19) county population estimates
#   - County-level Amish adherent counts (US Religion Census 2000/2010/2020)
# Outputs:
#   - ForPart2.RData containing the assembled df_final county-year panel
# Dependencies: R >= 4.4; tidyverse, readxl, readr, dplyr, stringr, broom, sf, glmmTMB
# Reproducibility: see sessionInfo() in the downstream analysis script
# Citation: Please cite the associated article (Chaney & Bergey)

library(tidyverse)
library(readxl)
library(broom)
library(mgcv)
library(dplyr)
library(readr)
library(stringr)


setwd(dir = "./Amish/NLCD/CountyData")

ohio_counties <- data.frame(
  COUNTY = c("Adams", "Allen", "Ashland", "Ashtabula", "Athens", "Auglaize", "Belmont", "Brown", "Butler", "Carroll",
             "Champaign", "Clark", "Clermont", "Clinton", "Columbiana", "Coshocton", "Crawford", "Cuyahoga", "Darke", "Defiance",
             "Delaware", "Erie", "Fairfield", "Fayette", "Franklin", "Fulton", "Gallia", "Geauga", "Greene", "Guernsey",
             "Hamilton", "Hancock", "Hardin", "Harrison", "Henry", "Highland", "Hocking", "Holmes", "Huron", "Jackson",
             "Jefferson", "Knox", "Lake", "Lawrence", "Licking", "Logan", "Lorain", "Lucas", "Madison", "Mahoning",
             "Marion", "Medina", "Meigs", "Mercer", "Miami", "Monroe", "Montgomery", "Morgan", "Morrow", "Muskingum",
             "Noble", "Ottawa", "Paulding", "Perry", "Pickaway", "Pike", "Portage", "Preble", "Putnam", "Richland",
             "Ross", "Sandusky", "Scioto", "Seneca", "Shelby", "Stark", "Summit", "Trumbull", "Tuscarawas", "Union",
             "Van Wert", "Vinton", "Warren", "Washington", "Wayne", "Williams", "Wood", "Wyandot"),
  FIPS = c("39001", "39003", "39005", "39007", "39009", "39011", "39013", "39015", "39017", "39019",
           "39021", "39023", "39025", "39027", "39029", "39031", "39033", "39035", "39037", "39039",
           "39041", "39043", "39045", "39047", "39049", "39051", "39053", "39055", "39057", "39059",
           "39061", "39063", "39065", "39067", "39069", "39071", "39073", "39075", "39077", "39079",
           "39081", "39083", "39085", "39087", "39089", "39091", "39093", "39095", "39097", "39099",
           "39101", "39103", "39105", "39107", "39109", "39111", "39113", "39115", "39117", "39119",
           "39121", "39123", "39125", "39127", "39129", "39131", "39133", "39135", "39137", "39139",
           "39141", "39143", "39145", "39147", "39149", "39151", "39153", "39155", "39157", "39159",
           "39161", "39163", "39165", "39167", "39169", "39171", "39173", "39175")
)

data_dir <- "./Amish/NLCD/CountyData"
files <- list.files(pattern = "\\.xlsx$", full.names = TRUE)

process_file <- function(file_path) {
  
  fname <- basename(file_path)
  parts <- str_remove(fname, "\\.xlsx$") %>% str_split("_") %>% unlist()
  
  year <- as.integer(parts[1])
  fips <- parts[2]
  
  df <- read_excel(file_path)
  
  value_col <- names(df)[3]
  
  df_wide <- df %>%
    rename(
      class = CLASSNAME,
      area = !!sym(value_col)
    ) %>%
    select(class, area) %>%
    pivot_wider(
      names_from = class,
      values_from = area
    ) %>%
    mutate(
      year = year,
      FIPS = fips
    )
  
  return(df_wide)
}

# combine all files
all_data <- map_dfr(files, process_file)

# add county names
all_data <- all_data %>%
  left_join(ohio_counties, by = "FIPS") %>%
  relocate(year, COUNTY, FIPS)

all_data$county_name <- all_data$COUNTY

all_data <- all_data %>% select(-COUNTY)

# read in Amish population data from '00, '10, and '20
amish_pop <- read.csv("../../Ohio_Amish_census_years.csv")

df_long <- amish_pop %>%
  pivot_longer(
    cols = starts_with("amish_share_"),
    names_to = "year",
    names_prefix = "amish_share_",
    values_to = "amish_share"
  ) %>%
  mutate(year = as.numeric(year)) %>%
  select(c(county_name, fips, year, amish_share))

# fit linear models by county
models <- df_long %>%
  group_by(county_name) %>%
  nest() %>%
  mutate(
    model = map(data, ~ lm(amish_share ~ year, data = .x)),
    tidy = map(model, tidy)
  )

ggplot(df_long, aes(x = year, y = amish_share, 
                    group = county_name, color = county_name)) +
  geom_point(size = 2) +
  geom_line() +
  theme_minimal() 

# estimate population numbers for all years
df_long <- df_long %>%
  mutate(
    year = as.numeric(year),
    county_name = as.factor(county_name)
  ) %>%
  drop_na()

# piecewise linear interpolation
df_interp <- df_long %>%
  group_by(county_name) %>%
  complete(year = 2000:2024) %>%
  arrange(year) %>%
  mutate(
    amish_share = approx(
      x = year[!is.na(amish_share)],
      y = amish_share[!is.na(amish_share)],
      xout = year,
      method = "linear",
      rule = 1 # no extrapolation
    )$y
  ) %>%
  ungroup()

# Extrapolate Amish share for years after 2020
df_extrap <- df_interp %>%
  group_by(county_name) %>%
  
  # Add all desired years
  complete(year = 2000:2024) %>%
  arrange(year) %>%
  
  # Interpolate observed gaps
  mutate(
    amish_share = approx(
      x = year[!is.na(amish_share)],
      y = amish_share[!is.na(amish_share)],
      xout = year,
      method = "linear",
      rule = 1
    )$y
  ) %>%
  
  # Compute 2010→2020 slope
  mutate(
    slope_2010_2020 =
      (amish_share[year == 2020] -
         amish_share[year == 2010]) / 10
  ) %>%
  
  # Extrapolate beyond 2020
  mutate(
    amish_share = ifelse(
      year > 2020,
      amish_share[year == 2020] +
        slope_2010_2020 * (year - 2020),
      amish_share
    )
  ) %>%
  
  ungroup()

df_extrap <- df_extrap %>%
  mutate(
    amish_share = pmax(pmin(amish_share, 1), 0)
  )

# Graph Holmes to check things a little bit
ggplot(df_extrap %>% filter(county_name == "Holmes"),
       aes(year, amish_share)) +
  geom_line(linewidth = 1) +
  geom_point(
    data = df_extrap %>% filter(county_name == "Holmes"),
    size = 2
  ) +
  theme_minimal()

# Final dataframe with everything looking nice
ohio_counties <- rename(ohio_counties, county_name = COUNTY)
df_extrap <- rename(df_extrap, FIPS = fips) %>% 
  select(-c(slope_2010_2020, FIPS))

all_data$FIPS <- as.integer(all_data$FIPS)

df_final <- df_extrap %>% 
  right_join(
    all_data, by = c("year", "county_name")
  ) 

#remove(amish_pop, df_extrap, df_interp, df_long, models, ohio_counties)

# Trim away rows with no NLCD data
df_final <- df_final %>% 
  filter(year > 2002)

## Process yearly LACV-ND case data from CDC
state_lookup <- data.frame(
  STATE_ABBR = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA",
                 "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD",
                 "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ",
                 "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC",
                 "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY", "DC")
)

ohio_counties <- data.frame(
  county_name = c("Adams", "Allen", "Ashland", "Ashtabula", "Athens", "Auglaize", "Belmont", "Brown", "Butler", "Carroll",
             "Champaign", "Clark", "Clermont", "Clinton", "Columbiana", "Coshocton", "Crawford", "Cuyahoga", "Darke", "Defiance",
             "Delaware", "Erie", "Fairfield", "Fayette", "Franklin", "Fulton", "Gallia", "Geauga", "Greene", "Guernsey",
             "Hamilton", "Hancock", "Hardin", "Harrison", "Henry", "Highland", "Hocking", "Holmes", "Huron", "Jackson",
             "Jefferson", "Knox", "Lake", "Lawrence", "Licking", "Logan", "Lorain", "Lucas", "Madison", "Mahoning",
             "Marion", "Medina", "Meigs", "Mercer", "Miami", "Monroe", "Montgomery", "Morgan", "Morrow", "Muskingum",
             "Noble", "Ottawa", "Paulding", "Perry", "Pickaway", "Pike", "Portage", "Preble", "Putnam", "Richland",
             "Ross", "Sandusky", "Scioto", "Seneca", "Shelby", "Stark", "Summit", "Trumbull", "Tuscarawas", "Union",
             "Van Wert", "Vinton", "Warren", "Washington", "Wayne", "Williams", "Wood", "Wyandot"),
  FIPS = c("39001", "39003", "39005", "39007", "39009", "39011", "39013", "39015", "39017", "39019",
           "39021", "39023", "39025", "39027", "39029", "39031", "39033", "39035", "39037", "39039",
           "39041", "39043", "39045", "39047", "39049", "39051", "39053", "39055", "39057", "39059",
           "39061", "39063", "39065", "39067", "39069", "39071", "39073", "39075", "39077", "39079",
           "39081", "39083", "39085", "39087", "39089", "39091", "39093", "39095", "39097", "39099",
           "39101", "39103", "39105", "39107", "39109", "39111", "39113", "39115", "39117", "39119",
           "39121", "39123", "39125", "39127", "39129", "39131", "39133", "39135", "39137", "39139",
           "39141", "39143", "39145", "39147", "39149", "39151", "39153", "39155", "39157", "39159",
           "39161", "39163", "39165", "39167", "39169", "39171", "39173", "39175")
)

# Get all LACV files
files <- list.files(
  path = "../../CDC annual data",
  pattern = "^LACV_cases_\\d{4}\\.csv$",
  full.names = TRUE
)

# Process each file
all_lacv <- lapply(files, function(file) {
  
  df <- read_csv(file, show_col_types = FALSE)
  
  df %>%
    separate_wider_delim(
      Location,
      delim = ", ",
      names = c("STATE_ABBR", "county_name")
    ) %>%
    filter(STATE_ABBR == "OH") %>%
    left_join(ohio_counties, by = "county_name") %>%
    transmute(
      year = Year,
      county_name,
      FIPS,
      LACV_cases = `Total human disease cases`
    )
  
}) %>%
  bind_rows()

# Merge LACV-ND cases with earlier dataframe
all_lacv$FIPS <- as.integer(all_lacv$FIPS)

df_final <- df_final %>% 
  left_join(
    all_lacv, by = c("FIPS", "year", "county_name")
  ) %>% 
  mutate(LACV_cases = ifelse(is.na(LACV_cases), 0, LACV_cases)) %>% 
  filter(year > 2002)

# Transform the NLCD data into proportions with dplyr
df_final <- df_final %>%
  mutate(
    total_pixels = rowSums(across(5:19), na.rm = TRUE)
  ) %>%
  mutate(
    across(
      5:19,
      ~ .x / total_pixels,
      .names = "{.col}_prop"
    )
  ) %>% 
  select(-total_pixels) %>% 
  select(-(5:19))

# Read in and summarize yearly county population estimates for subadults
yearly_pops <- read.csv("../../cc-est2010-2024-summary.csv")

df_final <- df_final %>% 
  left_join(
    yearly_pops, by = c("county_name", "year")
  )
#remove(yearly_pops)

total_pops <- read.csv("../../cc-est2010-2024-tot-kid_pops.csv")

df_final <- total_pops %>% 
  select(-c(tot_pop, kid_proportion)) %>%
  right_join(
    df_final, by = c("county_name", "year")
  )

# ---- Save the assembled panel for Part 2 ----
save(df_final, file = "../../ForPart2.RData")
