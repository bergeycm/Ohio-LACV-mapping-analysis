# Spatial correlates of La Crosse virus (LACV) in Ohio: Amish share & incidence

[![Preprint DOI](https://zenodo.org/badge/1133031822.svg)](https://doi.org/10.5281/zenodo.18228033)

**Main Scripts:** `01_process_annual_NLCD_data.R` and `02_ohio_lacv_amish_analysis.R` (v2.0.0)  
**Authors:** Morgan E. Chaney (Rutgers), Christina M. Bergey (Rutgers)  
**License:** MIT (see `LICENSE`)  
**Citation:** See `CITATION.cff` (preprint/article link TBD)


## Overview

This repository contains a workflow to reproduce the analyses and figures relating county-level La Crosse virus (LACV) incidence in Ohio to Amish population share and other covariates.

The workflow is two parts:

1. **`01_process_annual_NLCD_data.R`** builds the county-year analysis panel
   (`df_final`) from raw public sources and saves it as `ForPart2.RData`.
2. **`02_ohio_lacv_amish_analysis.R`** loads that panel plus the covariate tables and
   a county shapefile, then runs all models, sensitivity analyses, diagnostics, figures,
   and tables.


## Reproducing the results

### Fast path (recommended): run Part 2 from the provided panel
Everything Part 2 needs is in the repository (the county boundary shapefile downloads
itself on first run). From the repo root:
```
Rscript 02_ohio_lacv_amish_analysis.R
```

Session used is recorded in `sessionInfo_detailed.txt`.

### Full path: rebuild the panel from raw data
To regenerate `ForPart2.RData` from scratch, obtain the raw NLCD/CDC/Census inputs
listed below, set the paths at the top of `01_process_annual_NLCD_data.R`, and run it;
then run Part 2. The covariate tables can likewise be regenerated with the fetch scripts
under `Additional_covariate_datasets/` (see below).

## Data inputs

### Included in this repository
- `ForPart2.RData` - the assembled county-year panel `df_final` (output of Part 1; input to Part 2).
- `Additional_covariate_datasets/RUCC/RUCC_Ohio.csv` - USDA ERS Rural-Urban Continuum Codes, 2003/2013/2023 releases (FIPS, RUCC_2003, RUCC_2013, RUCC_2023).
- `Additional_covariate_datasets/climate/climate.csv` - annual county climate from PRISM (county_name, FIPS, year, temp_mean, precip, vpd). Regenerate with `fetch_prism_climate.R`.
- `Additional_covariate_datasets/SES/socioeconomic.csv` - annual county SES from the ACS 5-year (county_name, FIPS, year, median_income, pct_poverty, pct_college). Regenerate with `fetch_acs_ses.R`.

### Downloaded automatically
- **County boundary shapefile** `cb_2018_us_county_500k` - U.S. Census cartographic boundary file (1:500k). Part 2 downloads and unzips it on first run from the [Census Bureau](https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip); it is not stored in the repository.

### Required to rebuild the panel (Part 1 raw sources, not redistributed here)
- **Annual NLCD land-cover county tabulations** - [MRLC](https://www.mrlc.gov/) or [IPUMS NHGIS environmental summaries](https://www.nhgis.org/environmental-summaries).
- **CDC LACV-ND case counts by county-year** - CDC ArboNET.
- **U.S. Census total and child (ages 0-19) county population** - Census Bureau population estimates (CC-EST series).
- **County-level Amish adherent counts, 2000 / 2010 / 2020** - [U.S. Religion Census](https://www.usreligioncensus.org/). See the manuscript Methods.

## Software requirements

- R >= 4.4 (tested on 4.4.2; see `sessionInfo_detailed.txt`)
- Part 1: `tidyverse`, `readxl`, `readr`, `dplyr`, `stringr`, `broom`, `mgcv`, `sf`
- Part 2: `glmmTMB`, `sf`, `spdep`, `conleyreg`, `compositions`, `DHARMa`, `ggplot2`, `scales`, `ggrepel`, `cowplot`, `dplyr`, `tidyr`
- Covariate fetch scripts (optional): `tidycensus` (+ Census API key), `prism`, `terra`, `tigris`

```r
pkgs <- c("tidyverse","readxl","readr","broom","mgcv","sf","spdep",
          "glmmTMB","conleyreg","compositions","DHARMa","ggplot2",
          "scales","ggrepel","cowplot","tidyr")
install.packages(setdiff(pkgs, rownames(installed.packages())))
```

## AI-assistance disclosure

Analysis code was developed with assistance from Claude (Anthropic; Claude Opus 4 family,
including model identifier `claude-opus-4-8`). All code was reviewed and verified by the
authors, and every reported statistic was confirmed.

## Contact

Christina M. Bergey — christina.bergey@rutgers.edu  
Morgan E. Chaney

Lab website: https://www.bergey-lab.org/

## Changelog

- **v2.0.0** (2026-08): Reanalysis. Replaced the single-script spatial-model workflow with a
  two-part NB-GLMM pipeline (`01_process_annual_NLCD_data.R`, `02_ohio_lacv_amish_analysis.R`)
  on the 2003-2024 county-year panel; added compositional (ILR) land-use, SES, spatially-robust
  (Conley) SE, leave-one-out, and DHARMa sensitivity/diagnostics; the boundary shapefile now
  downloads on demand; included the assembled panel and covariate tables for direct reproducibility.
- **v1.0.0** (2026-01): Initial public release.
