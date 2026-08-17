#!/usr/bin/env Rscript
## =====================================================================
##  Fetch annual county-level climate (PRISM) for the Ohio LACV analysis
##  Outputs climate.csv with: county_name, FIPS, year, temp_mean, precip, vpd
##    temp_mean = annual mean temperature (deg C)
##    precip    = annual total precipitation (mm)
##    vpd       = annual mean vapor pressure deficit (hPa; mean of vpdmin & vpdmax)
##
##  Requirements (run once):
##    install.packages(c("prism", "terra", "sf", "tigris"))
##  Needs internet: downloads PRISM annual rasters (~4 vars x 22 yrs).
##  Run from the Amish_LACV directory (paths below are relative to it).
## =====================================================================

suppressMessages({library(prism); library(terra); library(sf)})

## ---- settings ----
years   <- 2003:2024
dl_dir  <- "Additional_covariate_datasets/climate/prism_downloads"   # raster cache
out_csv <- "Additional_covariate_datasets/climate/climate.csv"
dir.create(dl_dir, recursive = TRUE, showWarnings = FALSE)
prism_set_dl_dir(dl_dir)

## ---- 1. download PRISM annual grids (skips anything already downloaded) ----
## Annual products aggregate correctly per variable: tmean = mean, ppt = sum, vpd = mean.
## (function is get_prism_annual() in current prism; get_prism_annuals() in older versions)
dl_annual <- if ("get_prism_annual" %in% ls("package:prism")) get_prism_annual else get_prism_annuals
for (v in c("tmean", "ppt", "vpdmin", "vpdmax")) {
  message("Downloading PRISM annual: ", v)
  dl_annual(type = v, years = years, resolution = "4km", keepZip = FALSE)
}

## ---- 2. Ohio county polygons (name + FIPS), matched to PRISM's NAD83 lon/lat ----
## tigris NAME is bare county name ("Adams", "Holmes") to match df_final$county_name.
oh <- tigris::counties(state = "OH", cb = TRUE, year = 2020, progress_bar = FALSE)
oh <- sf::st_transform(oh, 4269)                 # NAD83 geographic, PRISM's CRS
oh_v <- terra::vect(oh)

## ---- 3. zonal county means for one variable across all years -> counties x years ----
county_annual <- function(type) {
  sapply(years, function(yr) {
    f <- pd_to_file(prism_archive_subset(type, "annual", years = yr, resolution = "4km"))
    r <- terra::rast(f)
    terra::extract(r, oh_v, fun = mean, na.rm = TRUE, ID = FALSE)[, 1]  # order = oh rows
  })
}
tmean <- county_annual("tmean"); ppt <- county_annual("ppt")
vmin  <- county_annual("vpdmin"); vmax <- county_annual("vpdmax")

## ---- 4. assemble long table and write ----
out <- do.call(rbind, lapply(seq_along(years), function(j) {
  data.frame(county_name = oh$NAME,
             FIPS = as.integer(oh$GEOID),
             year = years[j],
             temp_mean = round(tmean[, j], 3),
             precip    = round(ppt[, j], 1),
             vpd       = round((vmin[, j] + vmax[, j]) / 2, 3))
}))
out <- out[order(out$county_name, out$year), ]
write.csv(out, out_csv, row.names = FALSE)
cat(sprintf("Wrote %s : %d rows (%d counties x %d years)\n",
            out_csv, nrow(out), length(unique(out$county_name)), length(years)))
