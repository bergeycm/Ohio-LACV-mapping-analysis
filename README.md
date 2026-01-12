# Spatial correlates of La Crosse virus (LACV) in Ohio: Amish share & incidence

**Script:** `ohio_lacv_amish_analysis.R` (v1.0.0)  
**Authors:** Morgan E. Chaney (Rutgers), Christina M. Bergey (Rutgers)  
**License:** MIT (see `LICENSE`)  
**Citation:** See `CITATION.cff` (preprint/article link TBD)  

## Overview

This repository contains a single R workflow to reproduce the spatial analyses and figures relating county-level La Crosse virus (LACV) incidence in Ohio to Amish population share and land-use covariates. The script uses public datasets, fits spatial models, and generates figures (scatter/smoothers, bivariate choropleths, LISA cluster maps, etc.)

## Data inputs

Place these files in the working directory (or update paths in the script):

- Amish population (US Religious Census)
   - `2020_USRC_Amish_by_county.csv`
   - This file contains data available from the [U.S. Religion Census](https://www.usreligioncensus.org/). From that, we made this CSV file containing: State, FIPS, County, Total_Adh (total Amish adherents), Adherents_PCT_total_Adh (Amish adherents as proportion of total adherents included in survey), and Adherents_PCT_totalPop (adherents as a proportion of population).
- CDC LACV counts by county
   - `CDC_LACV_by_county.csv`
   - Included in repository.
- County boundaries (Census cartographic shapefile)
   - `cb_2018_us_county_500k/cb_2018_us_county_500k.shp` (with sidecar files)
   - Downloadable from the [U.S. Census website](https://www.census.gov/geographies/mapping-files/time-series/geo/carto-boundary-file.html)
- NLCD land-use summary by county
   - `combined_LC_data.csv`
   - Summaries by county available through the [IPUMS National Historical Geographic Information System](https://www.nhgis.org/environmental-summaries)
- CDC incidence table (per county)
   - `LACV_incidence_CDC_03-24.csv`
   - Included in repository.

## Software requirements

- R ≥ 4.4 (tested on 4.4.2)
- Packages used:
   - `dplyr`, `tidyr`, `ggplot2`, `sf`, `spdep`, `spatialreg`, `car`,
`ggrepel`, `viridis`, `scales`, `cowplot`, `patchwork`, `broom`,
`RColorBrewer`, `ggnewscale`

## Quick install
```
pkgs <- c(
  "dplyr","tidyr","ggplot2","sf","spdep","spatialreg","car",
  "ggrepel","viridis","scales","cowplot","patchwork","broom",
  "RColorBrewer","ggnewscale"
)
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

## How to run

From the repo root (where the script and data live):
```
Rscript ohio_lacv_amish_analysis.R
```

The script assumes a project root working directory. If you prefer another layout, edit the path lines near the top (e.g., where CSVs and the shapefile are read).

## Contact

Christina M. Bergey — christina.bergey@rutgers.edu  
Morgan E. Chaney

Lab website: https://www.bergey-lab.org/

## Changelog

v1.0.0 (2026-01-12): Initial public release.
