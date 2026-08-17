#!/usr/bin/env Rscript
## =====================================================================
##  Fetch county-level socioeconomic covariates (ACS 5-year) via tidycensus
##  Output: socioeconomic.csv with
##    county_name, FIPS, year, median_income, pct_poverty, pct_college
##      median_income = median household income (USD)
##      pct_poverty   = % of population below the poverty line
##      pct_college   = % of adults 25+ with a bachelor's degree or higher
##
##  Requirements:
##    install.packages("tidycensus")
##    Free Census API key: https://api.census.gov/data/key_signup.html
##    One-time: tidycensus::census_api_key("YOUR_KEY", install = TRUE)   # then restart R
##  Run from the Amish_LACV directory.
##
##  Note: the detailed education table (B15003) begins with the 2012 ACS 5-year,
##  so we pull ACS endyears 2012-2023 and map study years <2012 to 2012 and
##  2024 to 2023 (SES changes slowly; ACS 5-year windows already overlap).
## =====================================================================

suppressMessages({library(tidycensus); library(dplyr)})
if (Sys.getenv("CENSUS_API_KEY") == "")
  stop("Set a Census API key first: tidycensus::census_api_key('YOUR_KEY', install = TRUE), then restart R.")

out_csv     <- "Additional_covariate_datasets/SES/socioeconomic.csv"
acs_years   <- 2012:2023          # available ACS 5-year endyears with B15003
study_years <- 2003:2024

vars <- c(
  median_income = "B19013_001",   # median household income
  pov_below     = "B17001_002",   # income below poverty in past 12 months
  pov_total     = "B17001_001",   # population for whom poverty status is determined
  edu_total     = "B15003_001",   # population 25+
  edu_ba        = "B15003_022",   # bachelor's
  edu_ma        = "B15003_023",   # master's
  edu_prof      = "B15003_024",   # professional degree
  edu_phd       = "B15003_025"    # doctorate
)

pull_year <- function(yr) {
  get_acs(geography = "county", state = "OH", variables = vars,
          year = yr, survey = "acs5", output = "wide") |>
    transmute(
      FIPS          = as.integer(GEOID),
      county_name   = sub(" County, Ohio$", "", NAME),
      acs_year      = yr,
      median_income = median_incomeE,
      pct_poverty   = 100 * pov_belowE / pov_totalE,
      pct_college   = 100 * (edu_baE + edu_maE + edu_profE + edu_phdE) / edu_totalE
    )
}
acs <- bind_rows(lapply(acs_years, function(y) { message("ACS 5-yr endyear ", y); pull_year(y) }))

## map each study year to the nearest available ACS endyear (clamp to [2012, 2023])
map_yr <- pmin(pmax(study_years, min(acs_years)), max(acs_years))
out <- bind_rows(lapply(seq_along(study_years), function(i) {
  acs |>
    filter(acs_year == map_yr[i]) |>
    transmute(county_name, FIPS, year = study_years[i],
              median_income = round(median_income),
              pct_poverty   = round(pct_poverty, 2),
              pct_college   = round(pct_college, 2))
})) |>
  arrange(county_name, year)

write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("Wrote %s : %d rows (%d counties x %d years)\n",
            out_csv, nrow(out), length(unique(out$county_name)), length(study_years)))
