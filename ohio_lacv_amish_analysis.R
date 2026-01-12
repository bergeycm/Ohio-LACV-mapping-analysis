#!/usr/bin/env Rscript

# Title: ohio_lacv_amish_analysis.R (v1.0.0)
# Authors: 
#   Morgan E. Chaney (Rutgers University)
#   Christina M. Bergey (Rutgers University) <christina.bergey@rutgers.edu>
# Created: 2026-01-12
# License: MIT (see LICENSE file)
# Description: Runs all analyses for the spatial analysis of La Crosse virus in Ohio.
# Usage:
#   Rscript ohio_lacv_amish_analysis.R
# Inputs:
#   - Amish population dataset from US Religious Council
#       2020_USRC_Amish_by_county.csv
#   - CDC LACV data
#       CDC_LACV_by_county.csv
#   - Census data for total human population for each county
#       cb_2018_us_county_500k/cb_2018_us_county_500k.shp
#   - Land-use data for Ohio counties
#       combined_LC_data.csv
#   - Incidence data pulled directly from the CDC
#       LACV_incidence_CDC_03-24.csv
# Dependencies: R >= 4.4; packages listed below
# Reproducibility: See sessionInfo() at end of script
# Citation: Please cite the associated preprint/article

library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(spdep)
library(spatialreg)
library(car)
library(ggrepel)
library(viridis)      # colorblind-safe
library(scales)
library(cowplot)
library(patchwork)
library(broom)
library(RColorBrewer)
library(ggnewscale)

# ----------------------------------------------------------------------------------------
# --- Helpful definitions
# ----------------------------------------------------------------------------------------

# Optional: set seed for reproducibility
set.seed(123)

pct <- function(x) scales::percent(x/100, accuracy = 1)
fmt_int <- label_number(big.mark = ",", accuracy = 1)

# ----------------------------------------------------------------------------------------
# --- Create lookup tables 
# ----------------------------------------------------------------------------------------

state_lookup <- data.frame(
  STATE_ABBR = c("AL", "AK", "AZ", "AR", "CA", "CO", "CT", "DE", "FL", "GA",
                 "HI", "ID", "IL", "IN", "IA", "KS", "KY", "LA", "ME", "MD",
                 "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ",
                 "NM", "NY", "NC", "ND", "OH", "OK", "OR", "PA", "RI", "SC",
                 "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY", "DC"),
  STATE_NAME = c("Alabama", "Alaska", "Arizona", "Arkansas", "California", "Colorado",
                 "Connecticut", "Delaware", "Florida", "Georgia", "Hawaii", "Idaho",
                 "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana",
                 "Maine", "Maryland", "Massachusetts", "Michigan", "Minnesota",
                 "Mississippi", "Missouri", "Montana", "Nebraska", "Nevada",
                 "New Hampshire", "New Jersey", "New Mexico", "New York",
                 "North Carolina", "North Dakota", "Ohio", "Oklahoma", "Oregon",
                 "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota",
                 "Tennessee", "Texas", "Utah", "Vermont", "Virginia", "Washington",
                 "West Virginia", "Wisconsin", "Wyoming", "District of Columbia")
)

ohio_counties <- data.frame(
  COUNTY = c("Adams", "Allen", "Ashland", "Ashtabula", "Athens", "Auglaize", "Belmont",
             "Brown", "Butler", "Carroll", "Champaign", "Clark", "Clermont", "Clinton",
             "Columbiana", "Coshocton", "Crawford", "Cuyahoga", "Darke", "Defiance",
             "Delaware", "Erie", "Fairfield", "Fayette", "Franklin", "Fulton", "Gallia",
             "Geauga", "Greene", "Guernsey", "Hamilton", "Hancock", "Hardin", "Harrison",
             "Henry", "Highland", "Hocking", "Holmes", "Huron", "Jackson", "Jefferson",
             "Knox", "Lake", "Lawrence", "Licking", "Logan", "Lorain", "Lucas", "Madison",
             "Mahoning", "Marion", "Medina", "Meigs", "Mercer", "Miami", "Monroe",
             "Montgomery", "Morgan", "Morrow", "Muskingum", "Noble", "Ottawa", "Paulding",
             "Perry", "Pickaway", "Pike", "Portage", "Preble", "Putnam", "Richland",
             "Ross", "Sandusky", "Scioto", "Seneca", "Shelby", "Stark", "Summit",
             "Trumbull", "Tuscarawas", "Union", "Van Wert", "Vinton", "Warren",
             "Washington", "Wayne", "Williams", "Wood", "Wyandot"),
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

# ----------------------------------------------------------------------------------------
# --- Read and parse data
# ----------------------------------------------------------------------------------------

# --- Read in Amish population dataset from US Religious Council
amish_data_USRC <- read.csv("2020_USRC_Amish_by_county.csv")

# --- Read in CDC LACV data
lacv_data <- read.csv("CDC_LACV_by_county.csv")
lacv_data <- lacv_data %>%
  rename(FIPS = County) %>%
  full_join(amish_data_USRC, by = "FIPS")
lacv_data <- subset(lacv_data, select = -c(Year, Activity))

# --- Set NA values for Total.human.disease.cases and Neuroinvasive.disease.cases to 0
lacv_data <- lacv_data %>%
  mutate(Total.human.disease.cases =
             ifelse(is.na(Total.human.disease.cases), 0, Total.human.disease.cases),
         Neuroinvasive.disease.cases =
             ifelse(is.na(Neuroinvasive.disease.cases), 0, Neuroinvasive.disease.cases),
         Adherents_PCT_totalAdh =
             ifelse(is.na(Adherents_PCT_totalAdh), 0, Adherents_PCT_totalAdh),
         Adherents_PCT_totalPop =
             ifelse(is.na(Adherents_PCT_totalPop), 0, Adherents_PCT_totalPop),
         Total_Adh =
             ifelse(is.na(Total_Adh), 0, Total_Adh)) %>%
  mutate(Amish_presence = ifelse(Total_Adh == 0, 0, 1)) %>%
  # Separate FullGeoName into STATE_ABBR and COUNTY
  separate(FullGeoName, into = c("STATE_ABBR", "COUNTY"), sep = ", ", extra = "merge") %>%
  # Remove leading and trailing spaces
  mutate(STATE_ABBR = trimws(STATE_ABBR), COUNTY = trimws(COUNTY)) %>% 
  ## Join with state lookup to get full state names
  left_join(state_lookup, by = c("STATE_ABBR" = "STATE_ABBR")) %>%
  select(-STATE_ABBR) %>%
  rename(STATE = STATE_NAME)

# --- Read in Census data for total human population for each county
lacv_data <- subset(lacv_data, select = -c(FIPS, State, County)) %>% 
  left_join(ohio_counties, by = "COUNTY")
lacv_data <- subset(lacv_data, STATE == "Ohio")

sf <- read_sf("cb_2018_us_county_500k/cb_2018_us_county_500k.shp")
sf <- sf %>% mutate(FIPS = paste0(STATEFP, COUNTYFP)) %>% 
  rename(COUNTY = NAME)
sf$FIPS <- as.integer(sf$FIPS)
lacv_sf <- right_join(lacv_data, subset(sf, STATEFP == 39), by = "COUNTY" ) %>% 
  select(-FIPS.y) %>% 
  st_set_geometry("geometry") %>% 
  rename(FIPS = FIPS.x)

lacv_sf <- lacv_sf %>%
  mutate(Total.human.disease.cases =
             ifelse(is.na(Total.human.disease.cases), 0, Total.human.disease.cases),
         Neuroinvasive.disease.cases =
             ifelse(is.na(Neuroinvasive.disease.cases), 0, Neuroinvasive.disease.cases),
         Adherents_PCT_totalAdh =
             ifelse(is.na(Adherents_PCT_totalAdh), 0, Adherents_PCT_totalAdh),
         Adherents_PCT_totalPop =
             ifelse(is.na(Adherents_PCT_totalPop), 0, Adherents_PCT_totalPop),
         Amish_presence =
             ifelse(is.na(Total_Adh), 0, 1)) %>% 
  mutate(Total_Adh = ifelse(is.na(Total_Adh), 0, Total_Adh))

# --- Add in land-use data for Ohio counties
nlcd_data <- read.csv("combined_LC_data.csv")
lacv_sf$FIPS <- as.integer(paste0(lacv_sf$STATEFP,lacv_sf$COUNTYFP))
lacv_sf <- lacv_sf %>% left_join(nlcd_data, by = "FIPS") 

# --- Add in incidence data pulled directly from the CDC
CDC_incidence <- read.csv("LACV_incidence_CDC_03-24.csv")
lacv_sf$FIPS <- as.integer(lacv_sf$FIPS)
lacv_sf <- lacv_sf %>% left_join(CDC_incidence, by = "FIPS") %>% 
  mutate(LACV_incidence = ifelse(is.na(LACV_incidence), 0, LACV_incidence)) %>% 
  st_set_geometry("geometry")

# --- Do correlation screen for land-use variables that we want to use downstream
env_vars <- c("OpenWater", "DevelopedOpenSpace", "DevelopedLowIntensity", 
              "DevelopedMediumIntensity", "DevelopedHighIntensity", 
              "BarrenLand", "DeciduousForest", "EvergreenForest", "MixedForest",
              "ShrubScrub", "GrasslandHerbaceous", "PastureHay",
              "CultivatedCrops", "WoodyWetlands", "EmergentHerbaceousWetlands")

spearman.res = sapply(env_vars, function(x) 
  cor.test(lacv_sf[[x]], lacv_sf$LACV_incidence, 
           method="spearman", use="complete.obs"))

# --- Do Spatial Durbin Model
nb <- poly2nb(lacv_sf)
lw <- nb2listw(nb, style = "W")
moran.test(lacv_sf$LACV_incidence, lw)
formula <- log(LACV_incidence+1) ~ Adherents_PCT_totalPop + DeciduousForest + PastureHay
lm <- lm(formula, data = lacv_sf)
vif(lm) # Some multicollinearity, but nothing severe
rs_tests <- lm.RStests(lm, listw = lw, test = "all") 
summary(rs_tests)
sdm <- lagsarlm(formula, data = lacv_sf, listw = lw, type = "mixed")
summary(sdm)
# Get summary of impact measures (with Holmes)
imp <- impacts(sdm, listw = lw, R = 2000)
summary(imp, short=TRUE)

# --- Do SDM for a dataset that doesn't include Holmes County
no_holmes <- subset(lacv_sf, COUNTY != "Holmes")
lm_noholmes <- lm(formula, data = no_holmes)
vif(lm_noholmes) # Some multicollinearity, but nothing severe
no_holmes <- st_set_geometry(no_holmes, "geometry")
nb_noholmes <- poly2nb(no_holmes)
lw_noholmes <- nb2listw(nb_noholmes, style = "W")
rs_tests_noHolmes <- lm.RStests(lm_noholmes, listw = lw_noholmes, test = "all")
summary(rs_tests_noHolmes)
sdm_noholmes <- lagsarlm(formula, data = no_holmes, listw = lw_noholmes,
                         type = "mixed")
summary(sdm_noholmes)
# Get summary of impact measures (without Holmes)
imp_noholmes <- impacts(sdm_noholmes, listw = lw_noholmes, R = 2000)  # Monte Carlo SEs
summary(imp_noholmes, short = TRUE)

# Define df_all (+ Holmes flag) used later by the scatter overlay
df_all <- lacv_sf %>%
    st_drop_geometry() %>%
    mutate(is_holmes = COUNTY == "Holmes")

# ========================================================================================
# ---  Graphing 
# ========================================================================================

# ----------------------------------------------------------------------------------------
# --- Define color palette, themes
# ----------------------------------------------------------------------------------------

my.cols    = c("#8ECAE6", "#219EBC", "#023047", "#FFB703", "#FB8500")
my.cols.bw = c("#8ECAE6", "#219EBC", "#023047", "#FFB703", "#FB8500", "#000000", "#FFFFFF")

#                red        blue       yellow     darkred    darkblue
my.cols     = c("#C1133D", "#001C5A", "#FAEB08", "#78001C", "#00113A")
my.cols.bw  = c("#C1133D", "#001C5A", "#FAEB08", "#78001C", "#00113A", "#000000", "#FFFFFF")

# Red and blue are from Ohio flag
# Yellow is from Skyline Chili (#FAEB08)

# --- Discrete palettes (auto-extend if >5 categories)
my_pal_discrete <- function(n) {
  if (n <= length(my.cols)) my.cols[1:n]
  else grDevices::colorRampPalette(my.cols)(n)
}

scale_color_amishlacv_d <- function(..., na.value = "grey80") {
  ggplot2::discrete_scale("colour", "amishlacv_d",
                          palette = my_pal_discrete, na.value = na.value, ...)
}
scale_fill_amishlacv_d <- function(..., na.value = "grey80") {
  ggplot2::discrete_scale("fill", "amishlacv_d",
                          palette = my_pal_discrete, na.value = na.value, ...)
}

# --- Continuous palettes
# multi-stop gradient using full palette
scale_color_amishlacv_c <- function(..., reverse = FALSE) {
  ggplot2::scale_color_gradientn(colors = if (reverse) rev(my.cols) else my.cols, ...)
}
scale_fill_amishlacv_c <- function(..., reverse = FALSE) {
  ggplot2::scale_fill_gradientn(colors = if (reverse) rev(my.cols) else my.cols, ...)
}

# --- Diverging (for residuals centered at 0, etc.)
# pick low/mid/high from set; mid can be white or a mid color
scale_color_amishlacv_div <- function(midpoint = 0, low = my.cols[1],
                                   mid = "white", high = my.cols[2], ...) {
  ggplot2::scale_color_gradient2(low = low, mid = mid, high = high,
                                 midpoint = midpoint, ...)
}
scale_fill_amishlacv_div <- function(midpoint = 0, low = my.cols[1],
                                  mid = "white", high = my.cols[2], ...) {
  ggplot2::scale_fill_gradient2(low = low, mid = mid, high = high,
                                midpoint = midpoint, ...)
}

# --- 2-color continuous gradient between any two indices of my.cols.bw
scale_fill_amishlacv_pair <- function(a, b, ..., reverse = FALSE) {
  cols <- my.cols.bw[c(a, b)]
  if (reverse) cols <- rev(cols)
  ggplot2::scale_fill_gradientn(colors = cols, ...)
}
scale_color_amishlacv_pair <- function(a, b, ..., reverse = FALSE) {
  cols <- my.cols.bw[c(a, b)]
  if (reverse) cols <- rev(cols)
  ggplot2::scale_color_gradientn(colors = cols, ...)
}

# --- Some theme defaults

theme_set(theme_minimal(base_size = 12))
update_geom_defaults("text", list(family = "Helvetica"))
update_geom_defaults("label", list(family = "Helvetica"))

# ----------------------------------------------------------------------------------------
# --- 2x2 Amish x LACV contingency chart
# ----------------------------------------------------------------------------------------

tab_df <- lacv_sf %>%
  st_drop_geometry() %>%
  transmute(
    Amish = ifelse(Total_Adh > 0, "Present", "Absent"),
    LACV  = ifelse(LACV_incidence > 0, "Present", "Absent")
  ) %>%
  count(Amish, LACV, name = "Count")

tab_df2 <- tab_df %>%
  mutate(
    Amish = factor(Amish, levels = c("Absent","Present")),
    LACV  = factor(LACV,  levels = c("Absent","Present"))
  ) %>%
  complete(Amish, LACV, fill = list(Count = 0))   # adds Amish=Present & LACV=Absent if missing

p_con <- ggplot(tab_df2, aes(x = Amish, y = Count, fill = LACV)) +
  geom_col(position = position_dodge2(width = 0.8, preserve = "single"), width = 0.7) +
  scale_fill_manual(
    values = c("Absent" = my.cols[1], "Present" = my.cols[5]),
    drop = FALSE, name = "LACV"
  ) +
  labs(x = "Amish presence", y = "Number of counties",
       title = "Amish presence vs. recorded LACV-ND incidence") +
  theme(legend.position = "top")

#print(p_con)
ggsave("fig_S2_contingency_amish_lacv.pdf", p_con, width = 6, height = 4)

# ----------------------------------------------------------------------------------------
# --- Maps of LACV Incidence and Amish share
# ----------------------------------------------------------------------------------------

# Incidence
p_inc <- ggplot(lacv_sf) +
  geom_sf(aes(fill = LACV_incidence), color = "grey20", linewidth = 0.25) +
  scale_fill_amishlacv_pair(7, 2, na.value = "grey90",
                            name = "LACV-ND incidence\n(per 100k)") +
  labs(title = "LACV-ND incidence by county in Ohio") +
  geom_sf(data = subset(lacv_sf, COUNTY == "Holmes"),
            fill = NA, color = "white", linewidth = 0.5) +
  coord_sf(datum = NA)

#print(p_inc)
ggsave("fig_S1_map_incidence.pdf", p_inc, width = 6.5, height = 5)

# Amish share map (Holmes masked in black)

# Holmes value + pretty label
holmes_sf  <- subset(lacv_sf, COUNTY == "Holmes")
holmes_pct <- holmes_sf$Adherents_PCT_totalPop[1]
holmes_lab <- paste0("Holmes County (", percent(holmes_pct, accuracy = 1), ")")

# Base map: continuous gradient for all counties
p_amish <- lacv_sf %>%
  mutate(fill_amish = ifelse(COUNTY == "Holmes", NA, Adherents_PCT_totalPop)) %>%
  ggplot() +
  geom_sf(aes(fill = fill_amish), color = "grey20", linewidth = 0.25) +
  # continuous gradient; show legend ticks as percentages
  scale_fill_amishlacv_pair(
    7, 1,
    labels = percent_format(accuracy = 1),
    name   = NULL,
    guide  = guide_colorbar(order = 2)
  ) +
  labs(title = "Amish population share by county in Ohio") +
  coord_sf(datum = NA) +
  theme(legend.title = element_text())

# Start a NEW fill scale for the Holmes overlay + legend key
p_amish <- p_amish +
  ggnewscale::new_scale_fill() +
  geom_sf(data = holmes_sf, aes(fill = "holmes_key"), color = "grey20", linewidth = 0.2) +
  coord_sf(datum = NA) +
  scale_fill_manual(
    values = c(holmes_key = "black"),
    labels = c(holmes_key = holmes_lab),
    name   = "Amish population\n(% of total)",
    guide  = guide_legend(order = 1, override.aes = list(color = "grey20"))
  ) +
  # Red (%-rounded) label on top of Holmes
  geom_sf_text(
    data = holmes_sf,
    aes(label = percent(Adherents_PCT_totalPop, accuracy = 1)),
    color = "#f6a7b9", fontface = "bold", size = 4
  ) +
  # Outline around Holmes
  geom_sf(data = subset(lacv_sf, COUNTY == "Holmes"),
            fill = NA, color = "white", linewidth = 0.5) +
  coord_sf(datum = NA) +
  theme(
    axis.title  = element_blank(),
    axis.text   = element_blank(),
    axis.ticks  = element_blank()
  )

#print(p_amish)
ggsave("fig_S3_map_amish_share.pdf", p_amish, width = 6.5, height = 5)

# ----------------------------------------------------------------------------------------
# --- Scatterplot with and without Holmes County
# ----------------------------------------------------------------------------------------

# Build “with” and “without” data sources for the overlaid lines
df_with  <- df_all  %>% mutate(set = "All counties")
df_wo    <- df_all  %>% filter(!is_holmes) %>% mutate(set = "Excluding Holmes")

# Plot log1p, or the natural logarithm of 1 + x)
p_overlay <- ggplot(df_all,
                    aes(x = Adherents_PCT_totalPop,
                        y = log1p(LACV_incidence))) +
  # Fit line WITH Holmes
  geom_smooth(data = df_with, method = "lm", se = TRUE,
              aes(linetype = set), color = "black", fill=my.cols[3], alpha=0.2, lwd=0.5) +
  # Fit line WITHOUT Holmes
  geom_smooth(data = df_wo, method = "lm", se = TRUE,
              aes(linetype = set), color = "black", fill=my.cols[3], alpha=0.2, lwd=0.5) +
  # Points with white outline for clarity
  geom_point(data=df_with[df_with$is_holmes,], aes(color = is_holmes),
    size = 3, 
    shape = 21,
    fill = "white",
    stroke = 1.2,
    na.rm = TRUE
  ) +
  geom_point(data=df_with[df_with$is_holmes,], aes(color = is_holmes),
    size = 2.2, 
    shape = 16,
    na.rm = TRUE
  ) +
  geom_point(data=df_with[!df_with$is_holmes,], aes(color = is_holmes),
    size = 2.4, alpha = 0.5) +
  scale_color_manual(
    values = c(`TRUE` = my.cols[1], `FALSE` = my.cols[2]),
    guide  = "none"
  ) +
  # Line-type legend only, with clean labels
  scale_linetype_manual(
    values = c("All counties" = "solid", "Excluding Holmes" = "dotted"),
    breaks = c("All counties","Excluding Holmes"),
    labels = c("All counties","Excluding Holmes"),
    name   = NULL
  ) +
  guides(
    linetype = guide_legend(
      order = 1,
      override.aes = list(color = "black", linewidth = 0.4)  # show lines, not points
    )
  ) +
  scale_x_continuous(name = "Amish population share",
                     labels = scales::percent_format(accuracy = 1),
                     breaks = scales::pretty_breaks()) +
  scale_y_continuous(name = "LACV-ND incidence",
                     labels = function(b) format(expm1(b), trim = TRUE, big.mark = ",", digits = 3),
                     breaks = log1p(c(1:6))) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor     = element_blank(),
        legend.position      = c(0.03, 0.97),  # inset: upper-left
        legend.justification = c(0, 1),
        legend.background    = element_rect(fill = scales::alpha("white", 0.85), color = "black"),
        legend.key           = element_blank(),
        legend.margin        = margin(4, 6, 4, 6),
        legend.direction     = "vertical") #+
  #ggtitle("Amish share vs. LACV-ND incidence (log1p) with overlaid fits")

# Add label for Holmes
p_overlay <- p_overlay +
  geom_label_repel(
    data = subset(df_all, is_holmes),
    aes(label = paste0(COUNTY, " County")),
    label.padding = unit(0.25, "lines"),
    label.r = 0,
    fill = "white", color = my.cols[4],
    segment.color = my.cols[1], segment.size = 0.6,
    nudge_y = -0.4, nudge_x = -0.02, size = 3.5, seed = 123, max.overlaps = 100
  )

#print(p_overlay)
ggsave("fig_1_scatter_amish_vs_incidence_overlay.pdf", p_overlay, width = 6.5, height = 4)

# ----------------------------------------------------------------------------------------
# --- Bivariate choropleth showing co-distribution of Amish share and LACV-ND
# ----------------------------------------------------------------------------------------

# Prepare bivariate classes

df_biv <- lacv_sf %>%
  st_drop_geometry() %>%
  transmute(
    COUNTY,
    amish      = Adherents_PCT_totalPop,
    inc_log1p  = log1p(LACV_incidence),

    # Robust 3-bin categories that handle ties/zeros cleanly
    amish_i = if (n_distinct(amish, na.rm = TRUE) >= 3)
                ntile(amish, 3) else
                as.integer(cut(amish, breaks = pretty(amish, n = 3), include.lowest = TRUE)),

    inc_i   = if (n_distinct(inc_log1p, na.rm = TRUE) >= 3)
                ntile(inc_log1p, 3) else
                as.integer(cut(inc_log1p, breaks = pretty(inc_log1p, n = 3), include.lowest = TRUE)),

    amish_cat = factor(amish_i, levels = 1:3, labels = c("low","med","high")),
    inc_cat   = factor(inc_i,   levels = 1:3, labels = c("low","med","high")),
    biv_key   = paste0(amish_i, "-", inc_i)  # e.g., "1-3" = low Amish, high incidence
)

# Join classes back to sf
ohio_biv <- lacv_sf %>%
  left_join(df_biv, by = "COUNTY")

# Define a 3×3 bivariate color matrix
#    x = Amish share (left -> right low -> high)
#    y = Incidence  (bottom -> top low -> high)

# Ohio colors, muted
biv_pal <- c(
  "1-1" = "#E5E5EA",  # low var_pop,  low disease_incidence (soft gray)
  "2-1" = "#D9B5BC",  # mid var_pop,  low disease_incidence (dusty rose)
  "3-1" = "#B5697A",  # high var_pop, low disease_incidence (muted Ohio crimson)
  "1-2" = "#A8B5CD",  # low var_pop,  mid disease_incidence (muted slate blue)
  "2-2" = "#9D8A99",  # mid var_pop,  mid disease_incidence (muted mauve)
  "3-2" = "#8E5869",  # high var_pop, mid disease_incidence (dusty burgundy)
  "1-3" = "#4A6B8A",  # low var_pop,  high disease_incidence (muted Ohio navy)
  "2-3" = "#5D4A5E",  # mid var_pop,  high disease_incidence (dusty purple)
  "3-3" = "#4A3344"   # high var_pop, high disease_incidence (deep muted purple)
)

# Get Holmes County centroid for the point
holmes_centroid <- holmes_sf %>%
  st_centroid() %>%
  st_coordinates() %>%
  as.data.frame() %>%
  rename(lon = X, lat = Y)

# Create map
p_map <- ggplot(ohio_biv) +
  geom_sf(aes(fill = biv_key), color = "white", size = 0.2) +
  scale_fill_manual(values = biv_pal, na.value = "grey95",
                    guide = "none") +
  coord_sf() +
  
  # Outline around Holmes
  geom_sf(data = holmes_sf,
          fill = NA, color = my.cols[3], linewidth = 0.8) +
  
  # Point with outline for clarity (on top of the polygon outline)
  geom_point(data = holmes_centroid,
             aes(x = lon, y = lat),
             size = 3, 
             shape = 21,
             fill = "#4A3344",
             color = my.cols[3],
             stroke = 1.2,
             na.rm = TRUE) +
  geom_point(data = holmes_centroid,
             aes(x = lon, y = lat),
             size = 2.2, 
             shape = 16,
             color = my.cols[3],
             na.rm = TRUE) +
  
  theme_void(base_size = 12, base_family = "sans") +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 8)),
    plot.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(10, 10, 10, 10)
  )

# Create legend (small 3x3 tile grid)
legend_df <- expand.grid(
  amish_i = 1:3,
  inc_i   = 1:3
) %>%
  mutate(
    biv_key = paste0(amish_i, "-", inc_i),
    fill = biv_pal[biv_key]
  )

# Enhanced legend
p_legend <- ggplot(legend_df, aes(x = amish_i, y = inc_i)) +
  geom_tile(aes(fill = biv_key), color = "white", linewidth = 0.5) +
  scale_fill_manual(values = biv_pal, guide = "none") +
  scale_x_continuous(expand = c(0, 0.1), breaks = 1:3, 
                     labels = c("Low", "Med", "High"),
                     position = "top") +
  scale_y_continuous(expand = c(0, 0.1), breaks = 1:3, 
                     labels = c("Low", "Med", "High")) +
  labs(x = "Amish share", y = "Incidence") +
  theme_minimal(base_size = 9, base_family = "sans") +
  theme(
    axis.title.x = element_text(size = 9, margin = margin(b = 2)),
    axis.title.y = element_text(size = 9, margin = margin(r = 2)),
    axis.text = element_text(size = 8, color = "grey30"),
    panel.grid = element_blank(),
    plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
    plot.margin = margin(5, 5, 5, 5)
  )

# Combine with adjusted positioning
final_plot <- cowplot::ggdraw() +
  cowplot::draw_plot(p_map, 0, 0, 1, 1) +
  cowplot::draw_plot(p_legend, 0.75, 0.06, 0.22, 0.28)

#print(final_plot)
ggsave("fig_4_bivariate_choropleth_amish_incidence.pdf", final_plot, 
       width = 7.5, height = 5.5, bg = "white")

# ----------------------------------------------------------------------------------------
# --- LISA/Local Moran cluster map of LACV-ND
# ----------------------------------------------------------------------------------------

library(sf)
library(spdep)

# --- Neighbors and weights

# Queen contiguity; zero.policy=TRUE avoids errors if any islands exist
nb <- poly2nb(lacv_sf, queen = TRUE)
lw <- nb2listw(nb, style = "W", zero.policy = TRUE)

# --- Variable & spatial lag
y  <- log1p(lacv_sf$LACV_incidence)
z  <- as.numeric(scale(y))                                  # standardize variable
lag_z <- lag.listw(lw, z, zero.policy = TRUE)              # lag of standardized var

# --- Local Moran's I + pvals

lm_out <- localmoran(y, lw, zero.policy = TRUE)            # normal-approx p-values
# Robustly pick the p-value column (name varies by spdep version)
p_col  <- grep("^Pr", colnames(lm_out), value = TRUE)[1]
pvals  <- as.numeric(lm_out[, p_col])

# ---LISA cluster categories

cluster <- case_when(
  z >  0 & lag_z >  0 & pvals <= 0.05 ~ "High-High",
  z <  0 & lag_z <  0 & pvals <= 0.05 ~ "Low-Low",
  z >  0 & lag_z <  0 & pvals <= 0.05 ~ "High-Low",
  z <  0 & lag_z >  0 & pvals <= 0.05 ~ "Low-High",
  TRUE                                ~ "Not significant"
)

lisa_sf <- lacv_sf %>%
  mutate(
    LISA_cluster = factor(cluster,
      levels = c("High-High","Low-Low","High-Low","Low-High","Not significant")
    ),
    LISA_p = pvals,
    z_value = z,
    lag_z_value = lag_z
  )

# --- Plot LISA cluster map

# Queen contiguity; p <= 0.05
# High-High: high county surrounded by high neighbors;
# Not present: Low-Low: low with low neighbors; High-Low / Low-High: spatial outliers"

pal <- c(
  "High-High"       = "#e31a1c",  # red
  "Low-Low"         = "#1f78b4",  # blue
  "High-Low"        = "#fb9a99",  # light red
  "Low-High"        = "#a6cee3",  # light blue
  "Not significant" = "grey90"
)

p_lisa <- ggplot(lisa_sf) +
  geom_sf(aes(fill = LISA_cluster), color = "white", size = 0.25) +
  geom_sf(
    data = subset(lisa_sf, COUNTY == "Holmes"),
    fill = NA, color = "black", linewidth = 0.7
  ) +
  scale_fill_manual(values = pal, drop = FALSE) +
  guides(fill = guide_legend(title = "Local Moran cluster", ncol = 1)) +
  labs(
    title = "Local Moran (LISA) cluster map of LACV-ND incidence (log1p)",
  ) +
  theme_void(base_size = 12) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    plot.subtitle = element_text(margin = margin(b = 6))
  )

#print(p_lisa)
ggsave("fig_S4_LISA_LACV_incidence.pdf", p_lisa, width = 7.5, height = 5.5)

# Add Moran scatter
moran_df <- data.frame(z = z, lag_z = lag_z)
p_moran <- ggplot(moran_df, aes(x = z, y = lag_z)) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.7) +
  labs(x = "Standardized LACV-ND (log1p)",
       y = "Spatial lag of standardized LACV-ND",
       title = "Moran scatterplot for LACV-ND (Ohio counties)") +
  theme(panel.grid.minor = element_blank())

ggsave("fig_S5_LISA_LACV_incidence_w_moran_scatter.pdf", p_moran, width = 6, height = 4.8)

# ----------------------------------------------------------------------------------------
# --- Impacts (direct/indirect/total) decomposition plot from the SDM
# ----------------------------------------------------------------------------------------

# --- Impacts decomposition with bootstrap SE bars

library(spatialreg)

# Recompute with bootstrap if you haven't already (R = 2000 is good for CIs)
imp <- spatialreg::impacts(sdm, listw = lw, R = 2000)

library(tibble)

# Works when names(imp) == c("res", "sres"); falls back to imp$res if no sims.

# Get fallback term names from the model (non-lag, non-intercept)
coefs <- coef(sdm)
native_terms <- names(coefs)[!grepl("^lag\\.", names(coefs)) & names(coefs) != "(Intercept)"]

# Helper: coerce any object to a draws-by-terms matrix with column names
coerce_draws_matrix <- function(S, fallback_names) {
  M <- as.matrix(S)
  # If S collapsed to a vector (single term), make it a column matrix
  if (is.null(dim(M))) M <- matrix(M, ncol = 1)
  # Orient so that *rows = draws*, *cols = terms* (R is typically >> p)
  if (nrow(M) < ncol(M)) M <- t(M)
  # Attach term names if missing or empty
  if (is.null(colnames(M)) || any(colnames(M) == "") || any(is.na(colnames(M)))) {
    if (ncol(M) == length(fallback_names)) {
      colnames(M) <- fallback_names
    } else {
      colnames(M) <- paste0("var", seq_len(ncol(M)))
    }
  }
  M
}

# Extract simulated draws if present; else compute CIs from point estimates/SE
make_impacts_df <- function(imp, term_fallback) {
  # Prefer simulations from sres (newer spatialreg)
  if (!is.null(imp$sres) && is.list(imp$sres)) {
    effects_present <- intersect(names(imp$sres), c("direct","indirect","total"))
    stopifnot(length(effects_present) > 0)

    bind_rows(lapply(effects_present, function(effect_name) {
      M <- coerce_draws_matrix(imp$sres[[effect_name]], term_fallback)
      tibble(
        effect   = effect_name,
        term     = colnames(M),
        estimate = colMeans(M, na.rm = TRUE),
        se       = apply(M, 2, sd, na.rm = TRUE),
        lower    = apply(M, 2, quantile, probs = 0.025, na.rm = TRUE),
        upper    = apply(M, 2, quantile, probs = 0.975, na.rm = TRUE)
      )
    }))
  } else {
    # Fallback: try point estimates and (if available) SEs from imp$res
    # Structure varies; we try common patterns
    r <- imp$res
    # Expect r$direct, r$indirect, r$total as named numeric vectors
    effects_present <- intersect(names(r), c("direct","indirect","total"))
    stopifnot(length(effects_present) > 0)

    out <- bind_rows(lapply(effects_present, function(effect_name) {
      v <- r[[effect_name]]
      # Coerce to named vector with sensible names
      if (is.null(names(v)) || any(names(v) == "")) {
        if (length(v) == length(term_fallback)) {
          names(v) <- term_fallback
        } else {
          names(v) <- paste0("var", seq_along(v))
        }
      }
      tibble(effect = effect_name, term = names(v), estimate = as.numeric(v))
    }))

    # Try to find SEs (various packages store as r$se.direct, r$se_indirect, etc.)
    # If not found, leave CI as NA
    add_se <- function(df, r, eff) {
      se_candidates <- r[grep(paste0("^(se\\.?|sd\\.?|SE_)?", eff), names(r), ignore.case = TRUE)]
      if (length(se_candidates) == 0) return(df %>% mutate(se = NA_real_, lower = NA_real_, upper = NA_real_))
      se_vec <- unlist(se_candidates[[1]])
      if (is.null(names(se_vec)) || any(names(se_vec) == "")) {
        if (length(se_vec) == nrow(df)) names(se_vec) <- df$term
        else names(se_vec) <- df$term  # best effort
      }
      df %>%
        left_join(tibble(term = names(se_vec), se = as.numeric(se_vec)), by = "term") %>%
        mutate(lower = estimate - 1.96*se, upper = estimate + 1.96*se)
    }

    out <- bind_rows(
      add_se(filter(out, effect == "direct"),   r, "direct"),
      add_se(filter(out, effect == "indirect"), r, "indirect"),
      add_se(filter(out, effect == "total"),    r, "total")
    )
    out
  }
}

imp_df <- make_impacts_df(imp, native_terms)

# Drop any intercept that might sneak in
if ("(Intercept)" %in% imp_df$term) {
  imp_df <- dplyr::filter(imp_df, term != "(Intercept)")
}

# Order terms by total effect magnitude for readability
order_terms <- imp_df %>%
  filter(effect == "total") %>%
  arrange(desc(estimate)) %>%
  pull(term) %>%
  unique()

imp_df <- imp_df %>%
  mutate(
    term   = factor(term, levels = rev(order_terms)),
    effect = factor(effect, levels = c("direct","indirect","total"))
  )

# Niceify labels
lab_map <- c(
  Adherents_PCT_totalPop = "Amish population share (%)",
  DeciduousForest        = "Deciduous forest cover",
  PastureHay             = "Pasture/Hay cover"
)

# Apply nicified labels FIRST, before any ordering
imp_df <- imp_df %>%
  mutate(term_label = case_when(
    as.character(term) == "Adherents_PCT_totalPop dy/dx" ~ "Amish population share (%) dy/dx",
    as.character(term) == "DeciduousForest dy/dx" ~ "Deciduous forest cover dy/dx",
    as.character(term) == "PastureHay dy/dx" ~ "Pasture/Hay cover dy/dx",
    TRUE ~ as.character(term)
  ))

# THEN order by total effect magnitude using the new labels
order_terms_labeled <- imp_df %>%
  filter(effect == "total") %>%
  arrange(desc(estimate)) %>%
  pull(term_label) %>%
  unique()

imp_df <- imp_df %>%
  mutate(
    term_label = factor(term_label, levels = rev(order_terms_labeled)),
    effect = factor(effect, levels = c("direct","indirect","total"))
  )

# Plot with term_label instead of term
p_impacts <- ggplot(imp_df, aes(x = estimate, y = term_label, color = effect)) +
  # Reference line
  geom_vline(xintercept = 0, linetype = "solid", color = "grey60", linewidth = 0.4) +
  
  # Error bars with improved styling
  geom_errorbarh(
    aes(xmin = lower, xmax = upper),
    height = 0, 
    position = position_dodge(width = 0.65), 
    linewidth = 0.6,
    alpha = 0.8,
    na.rm = TRUE
  ) +
  
  # Points with white outline for clarity
  geom_point(
    position = position_dodge(width = 0.65), 
    size = 3, 
    shape = 21,
    fill = "white",
    stroke = 1.2,
    na.rm = TRUE
  ) +
  geom_point(
    position = position_dodge(width = 0.65), 
    size = 2.2, 
    shape = 16,
    na.rm = TRUE
  ) +
  
  # Color scale
  scale_color_manual(
    values = my.cols[c(1, 2, 3)], 
    name = "Effect type",
    labels = c("Direct", "Indirect (spillover)", "Total")
  ) +
  
  # Clean axis labels
  scale_x_continuous(
    breaks = scales::pretty_breaks(n = 6),
    expand = expansion(mult = c(0.05, 0.05))
  ) +
  
  labs(
    x = "Impact estimate (95% CI)", 
    y = NULL,
    title = "Spatial Durbin Model impacts",
    subtitle = "Effects on log(1 + LACV-ND incidence)"
  ) +
  
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(
    # Grid
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
    
    # Axes
    axis.text.y = element_text(size = 10, color = "grey20", hjust = 1),
    axis.text.x = element_text(size = 9, color = "grey30"),
    axis.title.x = element_text(face = "bold", size = 10, 
                                margin = margin(t = 10), color = "grey20"),
    axis.line.y = element_line(color = "grey40", linewidth = 0.4),
    axis.ticks.y = element_line(color = "grey40", linewidth = 0.3),
    axis.ticks.length.y = unit(0.1, "cm"),
    
    # Legend
    legend.position = "bottom",
    legend.title = element_text(face = "bold", size = 10, color = "grey20"),
    legend.text = element_text(size = 9, color = "grey20"),
    legend.key.height = unit(0.4, "cm"),
    legend.key.width = unit(1.2, "cm"),
    legend.margin = margin(t = 10),
    legend.box.spacing = unit(0.3, "cm"),
    
    # Title
    plot.title = element_text(face = "bold", size = 13, color = "grey20",
                             margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10, color = "grey40",
                                 margin = margin(b = 15)),
    
    # Overall plot
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA),
    plot.margin = margin(15, 15, 10, 10)
  )

print(p_impacts)
ggsave("fig_3_sdm_impacts_forest.pdf", p_impacts, 
       width = 7.5, height = 5.5, bg = "white")

# ----------------------------------------------------------------------------------------
# --- Partial-residual plots for key predictors (component+residual plots)
# ----------------------------------------------------------------------------------------

# Extract non-lag, non-intercept coefficient names present in the SDM
coefs <- coef(sdm)
native_terms <- names(coefs)[!grepl("^lag\\.", names(coefs)) & names(coefs) != "(Intercept)"]

# Make sure these columns exist in the data
X <- model.matrix(formula, data = st_drop_geometry(lacv_sf))
available_terms <- intersect(native_terms, colnames(X))

# Choose which to plot (defaults to all available RHS terms)
key_terms <- available_terms

fhat  <- fitted.values(sdm)      # linear predictor on log(1 + incidence) scale
resid <- residuals(sdm)

# Build a long data frame of partial residuals
pr_df <- bind_rows(lapply(key_terms, function(term) {
  beta_j <- coefs[term]
  xj     <- X[, term]
  # Component-plus-residual (partial residual)
  pr_j   <- resid + beta_j * xj
  tibble(
    COUNTY = lacv_sf$COUNTY,
    is_holmes = lacv_sf$COUNTY == "Holmes",
    term   = term,
    x      = xj,
    pr     = pr_j
  )
}))

# Nicify labels
lab_map <- c(
  Adherents_PCT_totalPop = "Amish population share (%)",
  DeciduousForest        = "Deciduous forest cover",
  PastureHay             = "Pasture/Hay cover"
)
pr_df$term_lab <- ifelse(pr_df$term %in% names(lab_map), lab_map[pr_df$term], pr_df$term)

# Plot one panel per term
p_pr <- ggplot(pr_df, aes(x = x, y = pr, color = is_holmes)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 0.6) +
  scale_color_manual(
    values = c(`TRUE` = my.cols[1], `FALSE` = my.cols[2]),
    guide  = "none"
  ) +
  facet_wrap(~ term_lab, scales = "free_x", ncol = 2) +
  labs(
    x = NULL,
    y = "Partial residual (log(1 + incidence))",
    title = "Partial-residual plots for SDM covariates",
    subtitle = "Points = counties; black line = linear fit; Holmes highlighted in red"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave("fig_S6_sdm_partial_residuals.pdf", p_pr, width = 8, height = 5.5)

# ----------------------------------------------------------------------------------------
# --- Residual map to visualize remaining spatial structure
# --- (SDM residual map + Moran's I on residuals)
# ----------------------------------------------------------------------------------------

lacv_sf$resid_sdm <- residuals(sdm)  # residuals on the model scale (log1p)

# Moran's I on residuals (row-standardized weights; safe for isolates)
mi_resid <- moran.test(lacv_sf$resid_sdm, lw, zero.policy = TRUE)

# --- Robust extraction of I, Z, p across versions
# Z-score lives in $statistic
mi_z <- tryCatch(unname(as.numeric(mi_resid$statistic[1])), error = function(e) NA_real_)
# p-value
mi_p <- tryCatch(as.numeric(mi_resid$p.value), error = function(e) NA_real_)
# Moran's I lives in $estimate (first element, or the one named like "Moran I statistic")
est_names <- names(mi_resid$estimate)
i_idx <- if (length(est_names)) {
  idx <- grep("Moran\\s*I", est_names, ignore.case = TRUE)
  if (length(idx)) idx[1] else 1
} else 1
mi_I <- tryCatch(unname(as.numeric(mi_resid$estimate[i_idx])), error = function(e) NA_real_)

mi_label <- sprintf("Moran's I = %.3f (Z = %.3f, p = %.3g)", mi_I, mi_z, mi_p)

# Centered diverging palette
rng <- max(abs(range(lacv_sf$resid_sdm, na.rm = TRUE)))
p_resid_map <- ggplot(lacv_sf) +
  geom_sf(aes(fill = resid_sdm), color = "white", linewidth = 0.25) +
  geom_sf(data = subset(lacv_sf, COUNTY == "Holmes"),
          fill = NA, color = "black", linewidth = 0.7) +
  geom_sf(data = subset(lacv_sf, COUNTY == "Holmes"),
          fill = NA, color = my.cols[3], linewidth = 0.5) +
  scale_fill_gradient2(
    low = my.cols[2], mid = "white", high = my.cols[1],
    midpoint = 0, limits = c(-rng, rng), name = "Residual"
  ) +
  labs(
    title = "Residual map from Spatial Durbin Model",
    subtitle = mi_label
  ) +
  coord_sf(datum = NA) +
  theme_void(base_size = 12) +
  theme(legend.position = "right")

ggsave("fig_S7_sdm_residual_map.pdf", p_resid_map, width = 7, height = 5.2)

# ----------------------------------------------------------------------------------------
# --- Observed vs. model-predicted scatter (with 1:1 line and R^2)
# ----------------------------------------------------------------------------------------

# Row-standardized W, identity
W   <- spdep::listw2mat(lw)
I_n <- diag(nrow(W))

# Spatial parameter (rho for lag/Durbin; fallback if needed)
rho <- if (!is.null(sdm$rho)) sdm$rho else as.numeric(coef(sdm)["rho"])

# Model matrix for your formula (includes (Intercept) column)
X_all <- model.matrix(formula, data = sf::st_drop_geometry(lacv_sf))

coef_all <- coef(sdm)

# --- Native (non-lag) coefficients: align strictly to columns in X_all
beta_all <- coef_all[!grepl("^lag\\.", names(coef_all))]

# Exclude any stray spatial params that might appear in coef() (e.g., "rho", "lambda")
beta <- beta_all[names(beta_all) %in% colnames(X_all)]

# Subset X to *exactly* those columns (keeps order consistent)
Xnat <- X_all[, names(beta), drop = FALSE]

# --- Lagged part: coefficients named like "lag.VarName" -> match to columns without intercept
theta_all <- coef_all[grepl("^lag\\.", names(coef_all))]
theta_vals <- as.numeric(theta_all)
names(theta_vals) <- sub("^lag\\.", "", names(theta_all))      # strip "lag."

# Keep only lag terms that have matching columns in X_all (and never the intercept)
theta_vals <- theta_vals[names(theta_vals) %in% colnames(X_all) &
                         names(theta_vals) != "(Intercept)"]

# Build WX for those lag columns (same order as theta_vals)
WXnat <- W %*% X_all[, names(theta_vals), drop = FALSE]

# --- Linear predictor on log1p scale via SAR-Durbin solution
eta_hat <- solve(I_n - rho * W, Xnat %*% beta + WXnat %*% theta_vals)

# Back-transform to incidence
yhat <- expm1(eta_hat)
obs <- lacv_sf$LACV_incidence
df_pred <- tibble(
  COUNTY = lacv_sf$COUNTY,
  observed = obs,
  predicted = as.numeric(yhat),
  is_holmes = lacv_sf$COUNTY == "Holmes"
)

# R^2 on log1p scale (in-sample), plus RMSE on original scale
R2_log  <- cor(log1p(obs), as.numeric(eta_hat))^2
RMSE    <- sqrt(mean((obs - yhat)^2, na.rm = TRUE))

p_obs_pred <- ggplot(df_pred, aes(x = observed, y = predicted, color = is_holmes)) +
  geom_point(, alpha = 0.8, size = 2) +
  scale_color_manual(values = c(`TRUE` = my.cols[1], `FALSE` = "grey30"), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = label_number(big.mark = ",")) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(
    x = "Observed LACV-ND incidence",
    y = "Predicted LACV-ND incidence",
    title = "Observed vs. Predicted (Spatial Durbin Model)",
    subtitle = sprintf("In-sample R² (log1p) = %.3f; RMSE = %.3f", R2_log, RMSE)
  ) +
  theme_minimal(base_size = 12)

ggsave("fig_S8_durbin_obs_vs_pred.pdf", p_obs_pred, width = 7, height = 5.2)

# ----------------------------------------------------------------------------------------
# --- Leave-one-out for Holmes: fit on all other counties, then predict Holmes - CUT
# ----------------------------------------------------------------------------------------

# Use the sdm_noholmes from above
# Apply its parameters to the full Ohio graph to get a Holmes prediction.

# Full W and model matrix for the *full* dataset
W_full <- spdep::listw2mat(lw)
I_n    <- diag(nrow(W_full))
X_full <- model.matrix(formula, data = sf::st_drop_geometry(lacv_sf))

# Parameters from the no-Holmes fit
coef_all <- coef(sdm_noholmes)
rho_hat  <- if (!is.null(sdm_noholmes$rho)) sdm_noholmes$rho else as.numeric(coef_all["rho"])

# Native (non-lag) part aligned to columns in X_full
beta_all <- coef_all[!grepl("^lag\\.", names(coef_all))]
beta     <- beta_all[names(beta_all) %in% colnames(X_full)]
Xnat     <- X_full[, names(beta), drop = FALSE]

# Lagged part: "lag.Var" -> Var; align to X_full (never use the intercept)
theta_all  <- coef_all[grepl("^lag\\.", names(coef_all))]
theta_vals <- as.numeric(theta_all)
names(theta_vals) <- sub("^lag\\.", "", names(theta_all))
theta_vals <- theta_vals[names(theta_vals) %in% colnames(X_full) & names(theta_vals) != "(Intercept)"]

WXnat <- W_full %*% X_full[, names(theta_vals), drop = FALSE]

# Predict (log1p), then back-transform to incidence for *all* counties using the no-Holmes params
eta_noHolmes_params <- solve(I_n - rho_hat * W_full, Xnat %*% beta + WXnat %*% theta_vals)
yhat_noHolmes_params <- expm1(eta_noHolmes_params)

# Read off Holmes
holmes_idx <- which(lacv_sf$COUNTY == "Holmes")
cat(sprintf("Holmes observed incidence = %.3f\n", lacv_sf$LACV_incidence[holmes_idx]))
cat(sprintf("Predicted (trained without Holmes) = %.3f\n",
            as.numeric(yhat_noHolmes_params[holmes_idx])))
cat(sprintf("Prediction error (pred - obs) = %.3f\n",
            as.numeric(yhat_noHolmes_params[holmes_idx]) - lacv_sf$LACV_incidence[holmes_idx]))

# ----------------------------------------------------------------------------------------
# --- Refit without Amish covariate, then show Obs-vs-Pred with Holmes highlighted
# ----------------------------------------------------------------------------------------

# Remove Amish term from the formula
fml_noAmish <- update(formula, . ~ . - Adherents_PCT_totalPop)

# Fit Spatial Durbin (mixed) without Amish
sdm_noAmish <- lagsarlm(fml_noAmish, data = lacv_sf, listw = lw, type = "mixed", zero.policy = TRUE)

# Robust prediction (same mechanics as above, but for the no-Amish model)
W   <- spdep::listw2mat(lw)
I_n <- diag(nrow(W))
X_all <- model.matrix(fml_noAmish, data = sf::st_drop_geometry(lacv_sf))

coef_all <- coef(sdm_noAmish)
rho2     <- if (!is.null(sdm_noAmish$rho)) sdm_noAmish$rho else as.numeric(coef_all["rho"])

beta_all <- coef_all[!grepl("^lag\\.", names(coef_all))]
beta     <- beta_all[names(beta_all) %in% colnames(X_all)]
Xnat     <- X_all[, names(beta), drop = FALSE]

theta_all  <- coef_all[grepl("^lag\\.", names(coef_all))]
theta_vals <- as.numeric(theta_all)
names(theta_vals) <- sub("^lag\\.", "", names(theta_all))
theta_vals <- theta_vals[names(theta_vals) %in% colnames(X_all) & names(theta_vals) != "(Intercept)"]

WXnat <- W %*% X_all[, names(theta_vals), drop = FALSE]

eta_hat_noAmish <- solve(I_n - rho2 * W, Xnat %*% beta + WXnat %*% theta_vals)
yhat_noAmish    <- expm1(eta_hat_noAmish)

df_noAmish <- tibble::tibble(
  COUNTY = lacv_sf$COUNTY,
  observed = lacv_sf$LACV_incidence,
  predicted = as.numeric(yhat_noAmish),
  is_holmes = lacv_sf$COUNTY == "Holmes"
)

R2_log_noAmish <- cor(log1p(df_noAmish$observed), as.numeric(eta_hat_noAmish))^2

p_obs_pred_noAmish <- ggplot(df_noAmish, aes(x = observed, y = predicted)) +
  geom_point(aes(color = is_holmes), size = 2, alpha = 0.9) +
  scale_color_manual(values = c(`TRUE` = my.cols[1], `FALSE` = "grey30"), guide = "none") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  scale_x_continuous(labels = label_number(big.mark = ",")) +
  scale_y_continuous(labels = label_number(big.mark = ",")) +
  labs(
    x = "Observed LACV-ND incidence",
    y = "Predicted incidence (model without Amish term)",
    title = "Observed vs. Predicted without Amish covariate",
    subtitle = sprintf("In-sample R² (log1p) = %.3f; Holmes highlighted", R2_log_noAmish)
  ) +
  theme_minimal(base_size = 12)

print(p_obs_pred_noAmish)

ggsave("fig_S9_durbin_obs_vs_pred_noAmish.pdf", p_obs_pred_noAmish, width = 6.5, height = 5)

# Print out Holmes under-prediction numerically
h <- which(df_noAmish$is_holmes)
cat(sprintf("Holmes observed = %.3f; predicted (no Amish) = %.3f; error = %.3f\n",
            df_noAmish$observed[h], df_noAmish$predicted[h],
            df_noAmish$predicted[h] - df_noAmish$observed[h]))

# ----------------------------------------------------------------------------------------
# --- Incidence by Amish-share quartile (violin + box, with counts)
# ----------------------------------------------------------------------------------------

pal <- my.cols

df0 <- sf::st_drop_geometry(lacv_sf) %>%
  mutate(is_zero = Adherents_PCT_totalPop == 0)

# Make quartiles only on the non-zero subset (ntile is tie-safe)
df_nz <- df0 %>%
  filter(!is_zero) %>%
  mutate(q = dplyr::ntile(Adherents_PCT_totalPop, 4)) %>%
  select(COUNTY, q)

# Join quartiles back; build 5-category factor
df_qnz <- df0 %>%
  left_join(df_nz, by = "COUNTY") %>%
  mutate(
    amish_cat = case_when(
      is_zero           ~ "0%",
      q == 1L           ~ "Q1",
      q == 2L           ~ "Q2",
      q == 3L           ~ "Q3",
      q == 4L           ~ "Q4",
      TRUE              ~ NA_character_
    ),
    amish_cat = factor(amish_cat, levels = c("0%","Q1","Q2","Q3","Q4"))
  )

# Optional: labels with counts
lvl <- levels(df_qnz$amish_cat)
counts <- df_qnz %>% count(amish_cat, name = "n")
base_labels <- c("0% Amish", "Q1", "Q2", "Q3", "Q4")
lab_with_n <- setNames(
  paste0(base_labels, "\n(n=", counts$n[match(lvl, counts$amish_cat)], ")"),
  lvl
)
df_qnz <- df_qnz %>% mutate(amish_lab = factor(amish_cat, levels = lvl, labels = lab_with_n[lvl]))

# (Optional) print bounds for the non-zero quartiles
nz_bounds <- df_qnz %>%
  filter(amish_cat %in% c("Q1","Q2","Q3","Q4")) %>%
  group_by(amish_cat) %>%
  summarise(
    min = min(Adherents_PCT_totalPop, na.rm = TRUE),
    max = max(Adherents_PCT_totalPop, na.rm = TRUE),
    .groups = "drop"
  )
cat(sprintf("Zero-Amish counties: n = %d\n", sum(df_qnz$is_zero, na.rm = TRUE)))
nz_bounds %>%
  arrange(amish_cat) %>%
  mutate(
    min = percent(min, accuracy = 0.1),
    max = percent(max, accuracy = 0.1)
  ) %>%
  { for (i in seq_len(nrow(.))) cat(sprintf("%s: %s to %s\n", .$amish_cat[i], .$min[i], .$max[i])) }

# Plot
p_quart_nz <- ggplot(df_qnz, aes(x = amish_lab, y = LACV_incidence)) +
  # Fatter violins
  geom_violin(trim = TRUE, width = 3, fill = pal[2], color = "black",
              alpha = 0.6, linewidth = 0.5) +
  # Narrower boxplots on top
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = pal[2], color = pal[1],
               alpha = 0.6, linewidth = 0.6) +
  # Keep jitter modest so it doesn't hide the violin edges
  geom_jitter(data=df_qnz[df_qnz$LACV_incidence < 2,],
              width = 0.08, height = 0, fill = "white", color = "white",
              alpha = 0.75, size = 0.6) +
  scale_y_continuous(labels = scales::label_number(big.mark = ",")) +
  labs(x = "Amish population share", y = "LACV-ND incidence") +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3)
  )

# Label Holmes
holmes_pt <- df_qnz %>%
  dplyr::filter(COUNTY == "Holmes") %>%
  dplyr::mutate(label = "Holmes County")

# add a highlighted point and a short labeled pointer
p_quart_nz <- p_quart_nz +
  # Point with white outline for clarity
  geom_point(data=holmes_pt, aes(x = amish_lab, y = LACV_incidence),
    size = 3 * 0.75, 
    shape = 21,
    fill = "white", col=pal[1],
    stroke = 1.2 * 0.75,
    na.rm = TRUE
  ) +
  geom_point(data=holmes_pt, aes(x = amish_lab, y = LACV_incidence),
    size = 2.2 * 0.6, 
    shape = 16,
    fill = pal[1], color = pal[1],
    na.rm = TRUE
  ) +
  ggrepel::geom_label_repel(
    data = holmes_pt,
    aes(x = amish_lab, y = LACV_incidence, label = label),
    nudge_y = -0.06,
    nudge_x = -1.15,
    size = 2.5,
    label.size = 0.2,                        # thin label border
    label.padding = unit(0.25, "lines"),
    label.r = 0,
    fill = "white", color = pal[4],
    segment.color = pal[4], segment.size = 0.6,
    min.segment.length = 0,                  # always draw a short line
    seed = 123, max.overlaps = Inf
  )

ggsave("fig_2_incidence_by_amish_zero_vs_quartiles_nonzero.pdf", p_quart_nz, width = 5, height = 4)

save.image("end_of_analysis.Rdata")

if (requireNamespace("sessioninfo", quietly = TRUE)) {
  capture.output(sessioninfo::session_info(pkgs = "attached"),
                 file = "output/sessionInfo_detailed.txt")
}

sessionInfo()
#  R version 4.4.2 (2024-10-31)
#  Platform: aarch64-apple-darwin20
#  Running under: macOS Sequoia 15.6.1
#  
#  locale:
#  [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
#  
#  attached base packages:
#  [1] stats     graphics  grDevices utils     datasets  methods   base     
#  
#  other attached packages:
#   [1] tibble_3.2.1       readr_2.1.5        reshape2_1.4.4     ggnewscale_0.5.2   RColorBrewer_1.1-3 broom_1.0.10      
#   [7] patchwork_1.3.2    cowplot_1.2.0      scales_1.4.0       viridis_0.6.5      viridisLite_0.4.2  ggrepel_0.9.6     
#  [13] car_3.1-3          carData_3.0-5      spatialreg_1.4-2   Matrix_1.7-1       spdep_1.4-1        spData_2.3.4      
#  [19] sf_1.0-21          ggplot2_4.0.0      tidyr_1.3.1        dplyr_1.1.4       
 