#!/usr/bin/env Rscript

# Title: Amish_LACV_full_analysis.R (v2.0.0)
# Authors:
#   Christina M. Bergey (Rutgers University) <christina.bergey@rutgers.edu>
#   Morgan E. Chaney (Rutgers University)
# Created: 2026-07 (revision reanalysis; supersedes ohio_lacv_amish_analysis.R v1.0.0)
# License: MIT (see LICENSE file)
# Usage:
#   Rscript Amish_LACV_full_analysis.R
# Inputs:
#   - Prepared county-year panel (df_final) from the NLCD/CDC/Census processing pipeline
#       as Rdata file with output of process_annual_NLCD_data.R
#   - County boundary shapefile (maps and spatial weights)
#       cb_2018_us_county_500k/*.shp
#   - Confounder covariates
#       Additional_covariate_datasets/{RUCC,climate,SES}/*.csv
# Dependencies: R >= 4.4; glmmTMB, sf, spdep, conleyreg, compositions, DHARMa,
#   ggplot2, scales, ggrepel, cowplot, dplyr, tidyr
# Reproducibility: set.seed(123) below; sessionInfo() printed at end of script
# Prepared with assistance of Anthropic Claude Opus 4.8, model identifier claude-opus-4-8

## =====================================================================
##  Ohio Amish population share & LACV-ND: full reanalysis
##  Negative-binomial GLMMs (county-year panel), sensitivity + diagnostics
##  Christina Bergey / Morgan Chaney — revision analysis
## =====================================================================

suppressMessages({library(glmmTMB); library(sf); library(spdep)})
set.seed(123)   # reproducibility (figure jitter / label repel)

## ---- paths ----------------------------------------------------------
rdata_path <- "ForPart2.RData" # Included in GitHub repo so you can start here

## ---- county boundary shapefile: download from U.S. Census if absent ----
shapefile <- "cb_2018_us_county_500k/cb_2018_us_county_500k.shp"
if (!file.exists(shapefile)) {
  url <- "https://www2.census.gov/geo/tiger/GENZ2018/shp/cb_2018_us_county_500k.zip"
  zip <- tempfile(fileext = ".zip")
  download.file(url, zip, mode = "wb")
  unzip(zip, exdir = "cb_2018_us_county_500k")
  unlink(zip)
}

## =====================================================================
## 0. Setup & data prep
## =====================================================================å
load(rdata_path)                          # provides df_final
df_final$county_name <- factor(df_final$county_name)
df_final$LACV_cases  <- as.integer(df_final$LACV_cases)
df_final$year        <- as.integer(df_final$year)
df_final$year_c      <- df_final$year - mean(df_final$year)

bt <- function(v) paste0("`", v, "`")     # backtick land-use names with spaces/slashes

## Collapse the four developed-intensity classes into one term due to near-collinearity (r ~ 0.9-0.99; VIF 35-100)
df_final$Developed_total <- with(df_final,
  `Developed Open Space_prop` + `Developed Low Intensity_prop` +
  `Developed Medium Intensity_prop` + `Developed High Intensity_prop`)

## Kitchen-sink "full >5%" set = every NLCD category exceeding 5% of max county land area,
## with the developed-intensity classes lumped into Developed_total and Cultivated Crops dropped
## as the compositional reference.
lu_full <- c("Open Water_prop", "Emergent Herbaceous Wetlands_prop", "Woody Wetlands_prop",
             "Pasture/Hay_prop", "Deciduous Forest_prop", "Mixed Forest_prop", "Developed_total")

## Vector larval ecology (VLE) informed set (a priori)
lu_vle <- c("Deciduous Forest_prop", "Mixed Forest_prop")

## Split developed-intensity classes, retained ONLY for the Section 2 VIF demonstration
lu_dev_split <- c("Developed Open Space_prop", "Developed Low Intensity_prop",
                  "Developed Medium Intensity_prop", "Developed High Intensity_prop")

## helper: fit an NB GLMM given a land-use string, offset variable, extras, and data
fit_nb <- function(lu_terms, offset_var = "tot_pop", extra = NULL, data = df_final) {
  rhs <- paste(c("year_c", "amish_share", lu_terms, extra,
                 "(1 | county_name)", sprintf("offset(log(%s))", offset_var)), collapse = " + ")
  glmmTMB(as.formula(paste("LACV_cases ~", rhs)), family = nbinom2, data = data)
}
## helper: pull the amish_share row
amish_row <- function(fit) {
  co <- summary(fit)$coefficients$cond; i <- grep("amish", rownames(co), ignore.case = TRUE)
  data.frame(beta = co[i, "Estimate"], se = co[i, "Std. Error"], p = co[i, "Pr(>|z|)"])
}

## ---- Confounders for the PRIMARY model (rurality + climate), merged HERE
conf <- character(0)
## Rurality (RUCC): NEAREST decennial vintage by year. Releases are 2003/2013/2023
rucc_path <- "Additional_covariate_datasets/RUCC/RUCC_Ohio.csv"
if (file.exists(rucc_path)) {
  rt <- read.csv(rucc_path); names(rt)[1] <- "FIPS"
  df_final <- merge(df_final, rt[, c("FIPS", "RUCC_2003", "RUCC_2013", "RUCC_2023")], by = "FIPS", all.x = TRUE)
  df_final$rucc <- with(df_final, ifelse(year <= 2008, RUCC_2003, ifelse(year <= 2018, RUCC_2013, RUCC_2023)))
  conf <- c(conf, "rucc")
}
## Climate (PRISM annual): mean temperature, total precipitation, mean VPD.
clim_path <- "Additional_covariate_datasets/climate/climate.csv"
if (file.exists(clim_path)) {
  cl <- read.csv(clim_path)
  df_final <- merge(df_final, cl[, c("FIPS", "year", "temp_mean", "precip", "vpd")], by = c("FIPS", "year"), all.x = TRUE)
  conf <- c(conf, "temp_mean", "precip", "vpd")
}
if (!length(conf)) conf <- NULL

## =====================================================================
## 1. The five models
##    1: per-capita, full LU
##    2: per-capita, VLE (vector larval ecology) set
##    3: per-capita, full LU + kid_prop
##    4: per-child,  full LU
##    5: per-child,  VLE set
## =====================================================================
m1 <- fit_nb(bt(lu_full),   "tot_pop")
m2 <- fit_nb(bt(lu_vle),    "tot_pop")
m3 <- fit_nb(bt(lu_full),   "tot_pop", extra = "kid_prop")
m4 <- fit_nb(bt(lu_full),   "kid_pop")
m5 <- fit_nb(bt(lu_vle),    "kid_pop")
models <- list(`1 percap/full` = m1, `2 percap/VLE` = m2, `3 percap/full+kid` = m3,
               `4 perchild/full` = m4, `5 perchild/VLE` = m5)

cat("\n=== Section 1: Amish coefficient by model ===\n")
amish_tab <- do.call(rbind, lapply(models, amish_row))
amish_tab <- cbind(model = names(models), round(amish_tab, 3)); print(amish_tab, row.names = FALSE)

cat("\n=== Section 1b: full model summaries ===\n")
for (nm in names(models)) { cat("\n---------------- Model", nm, "----------------\n"); print(summary(models[[nm]])) }

## =====================================================================
## 2. VIFs
## =====================================================================
vif_manual <- function(fit) {
  X <- model.matrix(fit)[, -1, drop = FALSE]
  v <- sapply(seq_len(ncol(X)), function(j) 1 / (1 - summary(lm(X[, j] ~ X[, -j]))$r.squared))
  round(setNames(v, colnames(X)), 2)
}
cat("\n=== Section 2: VIFs ===\n")
for (nm in names(models)) { cat("\n--", nm, "--\n"); print(sort(vif_manual(models[[nm]]), decreasing = TRUE)) }

## VIF demonstration: why the developed-intensity classes are lumped
cat("\n-- demonstration: kitchen sink with developed split into 4 intensity classes --\n")
m_devsplit <- fit_nb(bt(c(setdiff(lu_full, "Developed_total"), lu_dev_split)), "tot_pop")
print(sort(vif_manual(m_devsplit), decreasing = TRUE))

## =====================================================================
## 3. Zero-inflation check (is plain NB + county RE enough?)
## =====================================================================
zero_check <- function(fit, data = df_final) {
  exp0 <- sum(dnbinom(0, mu = fitted(fit), size = sigma(fit)))
  obs0 <- sum(data$LACV_cases == 0)
  cat(sprintf("observed zeros %d | NB-expected %.0f | ratio %.2f\n",
              obs0, exp0, obs0 / exp0))
}
cat("\n=== Section 3: Zero-inflation (per-capita primary model) ===\n")
zero_check(fit_nb(bt(lu_vle), "tot_pop", extra = conf))   # primary spec (VLE + climate + rurality)

## =====================================================================
## 4. Sensitivity / multiverse table for the Amish coefficient
##    offset x land-use set x Holmes, all on the PRIMARY spec (+ climate + rurality),
## =====================================================================
cat("\n=== Section 4: Amish coefficient across specifications (+ climate + rurality) ===\n")
lu_full_noWater <- setdiff(lu_full, "Open Water_prop")
specs <- expand.grid(offset = c("tot_pop", "kid_pop"),
                     lu = c("full", "vle", "full_noWater"),
                     holmes = c("in", "out"), stringsAsFactors = FALSE)
multi <- do.call(rbind, Map(function(off, lu, hol) {
  d  <- if (hol == "out") droplevels(subset(df_final, county_name != "Holmes")) else df_final
  tt <- switch(lu, full = lu_full, vle = lu_vle, full_noWater = lu_full_noWater)
  cbind(offset = off, lu = lu, holmes = hol, round(amish_row(fit_nb(bt(tt), off, extra = conf, data = d)), 3))
}, specs$offset, specs$lu, specs$holmes))
print(multi, row.names = FALSE)

## =====================================================================
## 5. Leverage checks — all on the per-capita PRIMARY model (VLE land use [Deciduous + Mixed
##    Forest] + climate + rurality), so the influence diagnostics match the primary specification.
## =====================================================================
cat("\n=== Section 5a: drop top-k Open Water (Lake Erie shore) counties (per-capita primary) ===\n")
ow <- sort(tapply(df_final[["Open Water_prop"]], df_final$county_name, mean), decreasing = TRUE)
cat("Top open-water counties:\n"); print(round(head(ow, 6), 4))
for (k in c(0, 1, 2, 3, 5)) {
  d <- if (k > 0) droplevels(subset(df_final, !(county_name %in% names(ow)[1:k]))) else df_final
  a <- amish_row(fit_nb(bt(lu_vle), "tot_pop", extra = conf, data = d))
  cat(sprintf("drop %d: amish beta %.2f (p=%.3f)\n", k, a$beta, a$p))
}

cat("\n=== Section 5b: leave-one-county-out, all 88 counties (per-capita primary) ===\n")
base <- amish_row(fit_nb(bt(lu_vle), "tot_pop", extra = conf))
cat(sprintf("base (all 88): amish beta %.2f (p=%.4f)\n", base$beta, base$p))
loco <- t(sapply(levels(df_final$county_name), function(cty) {
  a <- amish_row(fit_nb(bt(lu_vle), "tot_pop", extra = conf, data = droplevels(subset(df_final, county_name != cty))))
  c(beta = a$beta, p = a$p)
}))
loco <- data.frame(county = rownames(loco), beta = round(loco[, 1], 2), p = round(loco[, 2], 4),
                   delta = round(loco[, 1] - base$beta, 2))
cat("Most influential counties:\n"); print(head(loco[order(-abs(loco$delta)), ], 6), row.names = FALSE)
cat(sprintf("Holmes delta = %+.2f; max |delta| among the other 87 counties = %.2f\n",
    loco$delta[loco$county == "Holmes"], max(abs(loco$delta[loco$county != "Holmes"]))))

cat("\n=== Section 5c: leave-one-county-out among the 87 counties EXCLUDING Holmes (per-capita primary) ===\n")
d87    <- droplevels(subset(df_final, county_name != "Holmes"))
base87 <- amish_row(fit_nb(bt(lu_vle), "tot_pop", extra = conf, data = d87))
cat(sprintf("base (87, Holmes excluded): amish beta %.2f (p=%.4f)\n", base87$beta, base87$p))
loco87 <- t(sapply(levels(d87$county_name), function(cty) {
  a <- amish_row(fit_nb(bt(lu_vle), "tot_pop", extra = conf, data = droplevels(subset(d87, county_name != cty))))
  c(beta = a$beta, p = a$p)
}))
loco87 <- data.frame(county = rownames(loco87), beta = round(loco87[, 1], 2), p = round(loco87[, 2], 4))
cat(sprintf("coefficient positive in %d of %d fits; significant (p<0.05) in %d of %d; largest p = %.4f\n",
    sum(loco87$beta > 0), nrow(loco87), sum(loco87$p < 0.05), nrow(loco87), max(loco87$p)))
cat("least-significant removals:\n"); print(head(loco87[order(-loco87$p), ], 4), row.names = FALSE)

## =====================================================================
## 6. Figures
## =====================================================================
suppressMessages({library(ggplot2); library(scales); library(ggrepel); library(dplyr); library(tidyr); library(cowplot)})
dir.create("figures", showWarnings = FALSE)

## ---- palette + theme (inspired by Ohio flag and Skyline chili yellow) ----
##             red        blue       yellow     darkred    darkblue
my.cols    <- c("#C1133D", "#001C5A", "#FAEB08", "#78001C", "#00113A")
my.cols.bw <- c(my.cols, "#000000", "#FFFFFF")
scale_fill_amishlacv_pair <- function(a, b, ..., reverse = FALSE) {
  cols <- my.cols.bw[c(a, b)]; if (reverse) cols <- rev(cols)
  scale_fill_gradientn(colors = cols, ...)
}
theme_set(theme_minimal(base_size = 12))

## ---- county-level per-child aggregation (single row per county) ----
figdat <- aggregate(cbind(LACV_cases, kid_pop, tot_pop) ~ county_name + FIPS, df_final, sum)
figdat <- merge(figdat, aggregate(amish_share ~ county_name, df_final, mean))   # mean share over years
figdat$inc_child <- 1e5 * figdat$LACV_cases / figdat$kid_pop                    # per 100k child-years
figdat$inc_cap   <- 1e5 * figdat$LACV_cases / figdat$tot_pop                    # per 100k person-years
figdat$is_holmes <- figdat$county_name == "Holmes"

## ---- descriptive: county-level Spearman (incidence vs mean Amish share) ----
cat("\n=== Section 6 descriptive: Spearman county incidence vs mean Amish share ===\n")
for (v in c("inc_cap", "inc_child")) {
  ct <- cor.test(figdat[[v]], figdat$amish_share, method = "spearman", exact = FALSE)
  cat(sprintf("  %-9s: rho = %.3f, p = %.3g\n", v, ct$estimate, ct$p.value))
}

## ---- Fig 2: Amish share vs per-capita incidence, fits with & without Holmes ----
p1 <- ggplot(figdat, aes(amish_share, log1p(inc_cap))) +
  geom_smooth(aes(linetype = "All counties"), method = "lm", se = TRUE,
              color = "black", fill = my.cols[3], alpha = 0.2, linewidth = 0.5) +
  geom_smooth(data = subset(figdat, !is_holmes), aes(linetype = "Excluding Holmes"),
              method = "lm", se = TRUE, color = "black", fill = my.cols[3], alpha = 0.2, linewidth = 0.5) +
  geom_point(aes(color = is_holmes), size = 2.4, alpha = 0.6) +
  scale_color_manual(values = c(`TRUE` = my.cols[1], `FALSE` = my.cols[2]), guide = "none") +
  scale_linetype_manual(values = c("All counties" = "solid", "Excluding Holmes" = "dotted"), name = NULL) +
  scale_x_continuous("Amish population share", labels = percent_format(accuracy = 1)) +
  scale_y_continuous("LACV-ND incidence (per 100k population)",
                     breaks = log1p(c(0, 1, 2, 3)), labels = function(b) round(expm1(b))) +
  geom_label_repel(data = subset(figdat, is_holmes), aes(label = paste0(county_name, " County")),
                   fill = "white", color = my.cols[4], segment.color = my.cols[1],
                   nudge_y = -0.4, size = 3.2, seed = 123) +
  theme(legend.position = c(0.02, 0.98), legend.justification = c(0, 1),
        legend.background = element_rect(fill = alpha("white", 0.85), color = "black"))
ggsave("figures/fig_2_scatter_amish_vs_incidence.pdf", p1, width = 6.5, height = 4)

## ---- Fig 3: per-capita incidence by Amish-share quartile (0% + Q1-Q4) ----
figdat$is_zero <- figdat$amish_share == 0
nz <- figdat %>% filter(!is_zero) %>% mutate(q = ntile(amish_share, 4))
figdat <- figdat %>% left_join(select(nz, county_name, q), by = "county_name") %>%
  mutate(amish_cat = factor(case_when(is_zero ~ "0%", q == 1 ~ "Q1", q == 2 ~ "Q2",
                                      q == 3 ~ "Q3", q == 4 ~ "Q4"),
                            levels = c("0%", "Q1", "Q2", "Q3", "Q4")))
## x-axis labels with per-group county counts (computed from current data, not hard-coded)
grp <- c("0%", "Q1", "Q2", "Q3", "Q4"); cnt <- as.integer(table(figdat$amish_cat)[grp])
lab_n <- paste0(c("0% Amish", "Q1", "Q2", "Q3", "Q4"), "\n(n=", cnt, ")")
figdat$amish_lab <- factor(figdat$amish_cat, levels = grp, labels = lab_n)
p2 <- ggplot(figdat, aes(amish_lab, inc_cap)) +
  geom_violin(width = 1.1, fill = my.cols[2], color = "black", alpha = 0.5, trim = TRUE) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = my.cols[2], color = my.cols[1], alpha = 0.6) +
  ## semi-transparent white points for individual counties (Holmes shown separately below)
  geom_jitter(data = subset(figdat, !is_holmes), width = 0.08, height = 0,
              fill = "white", color = "white", alpha = 0.75, size = 0.9) +
  geom_point(data = subset(figdat, is_holmes), color = my.cols[1], size = 2.5) +
  geom_label_repel(data = subset(figdat, is_holmes), aes(label = "Holmes County"),
                   nudge_x = -1, size = 3, color = my.cols[4], seed = 123) +
  labs(x = "Amish population share", y = "LACV-ND incidence (per 100k population)")
ggsave("figures/fig_3_incidence_by_amish_quartile.pdf", p2, width = 5, height = 4)

## ---- Fig S3: Amish-presence x LACV-presence contingency ----
tab <- figdat %>% transmute(Amish = ifelse(amish_share > 0, "Present", "Absent"),
                            LACV  = ifelse(inc_child   > 0, "Present", "Absent")) %>%
  count(Amish, LACV, name = "Count") %>%
  complete(Amish = c("Absent", "Present"), LACV = c("Absent", "Present"), fill = list(Count = 0))
pS2 <- ggplot(tab, aes(Amish, Count, fill = LACV)) +
  geom_col(position = position_dodge2(preserve = "single"), width = 0.7) +
  scale_fill_manual(values = c("Absent" = my.cols[1], "Present" = my.cols[5])) +
  labs(x = "Amish presence", y = "Number of counties") + theme(legend.position = "top")
ggsave("figures/fig_S3_contingency.pdf", pS2, width = 6, height = 4)

## ---- Fig S4: Holmes leave-one-out, observed vs predicted from the non-Holmes model ----
## (per-capita PRIMARY model: VLE land use + climate + rurality, matching Section 5 and Section 8)
fitNoH <- fit_nb(bt(lu_vle), "tot_pop", extra = conf, data = droplevels(subset(df_final, county_name != "Holmes")))
df_final$pred_noH <- predict(fitNoH, newdata = df_final, type = "response",
                             re.form = NA, allow.new.levels = TRUE)
pp <- aggregate(cbind(LACV_cases, pred_noH, tot_pop) ~ county_name, df_final, sum)
pp$obs_inc  <- 1e5 * pp$LACV_cases / pp$tot_pop
pp$pred_inc <- 1e5 * pp$pred_noH  / pp$tot_pop
pp$is_holmes <- pp$county_name == "Holmes"
pLOO <- ggplot(pp, aes(log1p(obs_inc), log1p(pred_inc), color = is_holmes)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
  geom_point(size = 2.4, alpha = 0.8) +
  scale_color_manual(values = c(`TRUE` = my.cols[1], `FALSE` = "grey30"), guide = "none") +
  geom_label_repel(data = subset(pp, is_holmes), aes(label = "Holmes County"),
                   color = my.cols[4], seed = 123, size = 3) +
  scale_x_continuous("Observed incidence (per 100k population)",
                     breaks = log1p(c(0, 1, 3, 10)), labels = function(b) round(expm1(b))) +
  scale_y_continuous("Predicted from non-Holmes model",
                     breaks = log1p(c(0, 1, 3, 10)), labels = function(b) round(expm1(b))) +
  labs(title = "Leave-one-out: Holmes predicted from all other counties")
ggsave("figures/fig_S4_holmes_leave_one_out.pdf", pLOO, width = 6.5, height = 5)

## ---- Maps + LISA ----
if (file.exists(shapefile)) {
  sfc <- subset(read_sf(shapefile), STATEFP == "39")
  sfc$FIPS <- as.integer(paste0(sfc$STATEFP, sfc$COUNTYFP))
  cty_sf <- left_join(sfc, figdat, by = "FIPS")
  holmes_sf <- subset(cty_sf, county_name == "Holmes")

  ggsave("figures/fig_S1_map_incidence.pdf",
    ggplot(cty_sf) + geom_sf(aes(fill = inc_cap), color = "grey20", linewidth = 0.25) +
      scale_fill_amishlacv_pair(7, 2, name = "Incidence\n(per 100k pop.)") +
      geom_sf(data = holmes_sf, fill = NA, color = "white", linewidth = 0.5) +
      coord_sf(datum = NA) + labs(title = "LACV-ND incidence (per-capita) by county"),
    width = 6.5, height = 5)

  ggsave("figures/fig_S2_map_amish_share.pdf",
    ggplot(cty_sf) + geom_sf(aes(fill = amish_share), color = "grey20", linewidth = 0.25) +
      scale_fill_amishlacv_pair(7, 1, labels = percent_format(accuracy = 1), name = "Amish share") +
      geom_sf(data = holmes_sf, fill = NA, color = "white", linewidth = 0.5) +
      coord_sf(datum = NA) + labs(title = "Amish population share by county"),
    width = 6.5, height = 5)

  biv <- st_drop_geometry(cty_sf) %>%
    mutate(ai = ntile(amish_share, 3), ii = ntile(log1p(inc_cap), 3), biv_key = paste0(ai, "-", ii))
  cty_sf$biv_key <- biv$biv_key[match(cty_sf$county_name, biv$county_name)]
  ## 3x3 bivariate palette (x = Amish share low->high; y = incidence low->high)
  biv_pal <- c("1-1"="#E5E5EA","2-1"="#D9B5BC","3-1"="#B5697A","1-2"="#A8B5CD","2-2"="#9D8A99",
               "3-2"="#8E5869","1-3"="#4A6B8A","2-3"="#5D4A5E","3-3"="#4A3344")
  holmes_centroid <- as.data.frame(st_coordinates(st_centroid(st_geometry(holmes_sf))))
  names(holmes_centroid) <- c("lon", "lat")
  p_map <- ggplot(cty_sf) +
    geom_sf(aes(fill = biv_key), color = "white", linewidth = 0.2) +
    scale_fill_manual(values = biv_pal, na.value = "grey95", guide = "none") +
    geom_sf(data = holmes_sf, fill = NA, color = my.cols[3], linewidth = 0.8) +
    geom_point(data = holmes_centroid, aes(lon, lat), size = 3, shape = 21,
               fill = "#4A3344", color = my.cols[3], stroke = 1.2) +
    geom_point(data = holmes_centroid, aes(lon, lat), size = 2.2, shape = 16, color = my.cols[3]) +
    coord_sf(datum = NA) + theme_void(base_size = 12) +
    labs(title = "Amish share x LACV-ND incidence (per-capita)") +
    theme(plot.title = element_text(face = "bold", size = 13, margin = margin(b = 6)),
          plot.background = element_rect(fill = "white", color = NA), plot.margin = margin(10, 10, 10, 10))
  ## small 3x3 tile legend, overlaid via cowplot
  legend_df <- expand.grid(amish_i = 1:3, inc_i = 1:3)
  legend_df$biv_key <- paste0(legend_df$amish_i, "-", legend_df$inc_i)
  p_legend <- ggplot(legend_df, aes(amish_i, inc_i)) +
    geom_tile(aes(fill = biv_key), color = "white", linewidth = 0.5) +
    scale_fill_manual(values = biv_pal, guide = "none") +
    scale_x_continuous(expand = c(0, 0.1), breaks = 1:3, labels = c("Low", "Med", "High"), position = "top") +
    scale_y_continuous(expand = c(0, 0.1), breaks = 1:3, labels = c("Low", "Med", "High")) +
    labs(x = "Amish share", y = "Incidence") +
    theme_minimal(base_size = 9) +
    theme(axis.title.x = element_text(size = 9, margin = margin(b = 2)),
          axis.title.y = element_text(size = 9, margin = margin(r = 2)),
          axis.text = element_text(size = 8, color = "grey30"), panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", color = "grey70", linewidth = 0.5),
          plot.margin = margin(5, 5, 5, 5))
  p4 <- cowplot::ggdraw() + cowplot::draw_plot(p_map, 0, 0, 1, 1) +
        cowplot::draw_plot(p_legend, 0.74, 0.05, 0.23, 0.29)
  ggsave("figures/fig_1_bivariate_choropleth.pdf", p4, width = 7, height = 5.5, bg = "white")

  nb <- poly2nb(cty_sf, queen = TRUE); lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  y <- log1p(cty_sf$inc_cap); z <- as.numeric(scale(y)); lagz <- lag.listw(lw, z, zero.policy = TRUE)  # per-capita (Figs S5/S6)
  ## global Moran's I of incidence (motivates the LISA), reported for both offsets
  for (v in c("inc_cap", "inc_child")) {
    gmi <- moran.test(log1p(cty_sf[[v]]), lw, zero.policy = TRUE)
    cat(sprintf("Global Moran's I, LACV-ND incidence (log1p %s, queen): I = %.3f, p = %.3g\n",
                v, gmi$estimate["Moran I statistic"], gmi$p.value))
  }
} else {
  cat("\nSection 6: non-spatial figures written to figures/ (maps skipped: shapefile not found).\n")
}

## =====================================================================
## 7. Misc. sensitivity/influence analyses
##    7a <5% land-area threshold sensitivity
##    7b ILR compositional sensitivity
##    7c Residual/influence diagnostics (DHARMa)
## =====================================================================
prop_cols <- setdiff(grep("_prop$", names(df_final), value = TRUE), "kid_prop")  # land-cover props only (exclude child proportion which annoyingly has prop in name)
amish_of <- function(fit) {
  co <- summary(fit)$coefficients$cond; i <- grep("amish", rownames(co), ignore.case = TRUE)
  sprintf("%.2f (SE %.2f, p=%.3f)", co[i, "Estimate"], co[i, "Std. Error"], co[i, "Pr(>|z|)"])
}

## ---- 7a. Inclusion-threshold sensitivity (developed lumped; per-capita & per-child; Table 3) ----
## Candidate land classes exclude Cultivated Crops (compositional reference) and the four
## developed-intensity classes, which enter as the single lumped Developed_total term. Climate +
## rurality (conf) are included.
cat("\n=== 7a. Inclusion-threshold sensitivity (developed lumped, Cultivated Crops ref; + climate + rurality) ===\n")
lc_thr <- setdiff(prop_cols, c("Cultivated Crops_prop", "Developed Open Space_prop",
            "Developed Low Intensity_prop", "Developed Medium Intensity_prop", "Developed High Intensity_prop"))
for (off in c("tot_pop", "kid_pop")) {
  for (thr in c(0.05, 0.02, 0)) {
    keep <- c(bt(lc_thr[sapply(df_final[lc_thr], max) >= thr]), "Developed_total")
    fit  <- try(fit_nb(keep, off, extra = conf), silent = TRUE)
    cat(sprintf("  per-%-6s  >%-4s land-use set (%2d LU terms):  amish = %s\n",
                ifelse(off == "tot_pop", "capita", "child"), paste0(100 * thr, "%"), length(keep),
                if (inherits(fit, "try-error")) "did not converge" else amish_of(fit)))
  }
}

## ---- 7b. ILR (isometric log-ratio) compositional sensitivity, ACROSS inclusion thresholds ----
cat("\n=== 7b. ILR compositional treatment of land use, across inclusion thresholds ===\n")
suppressMessages(library(compositions))
ilr_amish <- function(thr) {
  lu   <- prop_cols[sapply(df_final[prop_cols], max) >= thr]
  comp <- as.matrix(df_final[lu]); comp[comp == 0] <- min(comp[comp > 0]) / 2   # simple zero replacement
  comp <- comp / rowSums(comp)
  ic   <- as.data.frame(unclass(ilr(acomp(comp)))); names(ic) <- paste0("ilr", seq_len(ncol(ic)))
  d    <- cbind(df_final, ic)
  m    <- glmmTMB(as.formula(paste("LACV_cases ~ year_c + amish_share +", paste(names(ic), collapse = " + "),
                  "+", paste(conf, collapse = " + "), "+ (1|county_name) + offset(log(tot_pop))")),
                  family = nbinom2, data = d)
  co <- summary(m)$coefficients$cond["amish_share", ]
  cat(sprintf("  ILR >%-4s (%2d classes, %2d coords): amish = %.2f (SE %.2f, p=%.4f)\n",
              paste0(100 * thr, "%"), length(lu), ncol(ic), co[1], co[2], co[4]))
  invisible(m)
}
m_ilr <- ilr_amish(0.05)          # >5% ILR retained as m_ilr for any downstream reference
for (thr in c(0.02, 0)) ilr_amish(thr)

## ---- 7c. Residual & influence diagnostics (DHARMa), primary model (VLE + climate + rurality) ----
## Report the full DHARMa panel and also refit EXCLUDING Holmes
cat("\n=== 7c. DHARMa simulated-residual diagnostics (per-capita primary model) ===\n")
if (requireNamespace("DHARMa", quietly = TRUE)) {
  dharma_panel <- function(fit, tag) {
    set.seed(123); s <- DHARMa::simulateResiduals(fit, n = 500)   # seeded: reproducible p-values (cited in paper)
    cat(sprintf("  %-18s KS(uniformity)=%.3f  quantile=%.3f  dispersion=%.3f  outlier=%.3f  zeroinfl=%.3f\n",
        tag, DHARMa::testUniformity(s, plot = FALSE)$p.value,
        tryCatch(DHARMa::testQuantiles(s, plot = FALSE)$p.value, error = function(e) NA_real_),
        DHARMa::testDispersion(s, plot = FALSE)$p.value,
        DHARMa::testOutliers(s, plot = FALSE)$p.value,
        DHARMa::testZeroInflation(s, plot = FALSE)$p.value))
    s
  }
  m_prim_diag <- fit_nb(bt(lu_vle), "tot_pop", extra = conf)   # == primary spec (Section 8 m_primary)
  sim_all <- dharma_panel(m_prim_diag, "all 88 counties")
  ## Holmes-excluded refit (tighter optimizer to avoid an intermittent step-failure warning), same spec
  frm_diag   <- as.formula(paste("LACV_cases ~", paste(c("year_c", "amish_share", bt(lu_vle), conf,
                  "(1 | county_name)", "offset(log(tot_pop))"), collapse = " + ")))
  m_diag_noH <- glmmTMB(frm_diag, family = nbinom2,
                        data = droplevels(subset(df_final, county_name != "Holmes")),
                        control = glmmTMBControl(optCtrl = list(iter.max = 1e4, eval.max = 1e4)))
  sim_noH <- dharma_panel(m_diag_noH, "excluding Holmes")

  ## ---- Fig S5: styled 2x2 DHARMa figure (manuscript palette; Holmes included top, excluded bottom) ----
  pfmt <- function(p) if (p < 0.001) "< 0.001" else sprintf("= %.3f", p)
  qq_panel <- function(s, title) {                                    # uniformity QQ (residuals vs uniform)
    r <- sort(s$scaledResiduals); ks <- DHARMa::testUniformity(s, plot = FALSE)$p.value
    ggplot(data.frame(exp = ppoints(length(r)), obs = r), aes(exp, obs)) +
      geom_abline(slope = 1, intercept = 0, color = my.cols[1], linewidth = 0.6) +
      geom_point(color = my.cols[2], alpha = 0.25, size = 0.9) +
      labs(x = "Expected (uniform)", y = "Observed residual", title = title,
           subtitle = bquote(KS ~ italic(p) ~ .(pfmt(ks)) ~ .(if (ks < 0.05) " (significant)" else " (n.s.)"))) +
      coord_equal() + theme_minimal(base_size = 11) +
      theme(plot.subtitle = element_text(color = if (ks < 0.05) my.cols[1] else "grey30"), panel.grid.minor = element_blank())
  }
  res_panel <- function(s, title) {                                   # residual vs rank-transformed prediction
    x <- rank(s$fittedPredictedResponse, ties.method = "average") / length(s$scaledResiduals)
    qp <- DHARMa::testQuantiles(s, plot = FALSE)$p.value
    ggplot(data.frame(x = x, y = s$scaledResiduals), aes(x, y)) +
      geom_hline(yintercept = c(.25, .5, .75), linetype = "dashed", color = "grey70", linewidth = 0.3) +
      geom_point(color = my.cols[2], alpha = 0.22, size = 0.9) +
      geom_quantile(quantiles = c(.25, .5, .75), formula = y ~ splines::ns(x, 4), method = "rq",
                    color = my.cols[1], linewidth = 0.7) +
      labs(x = "Model predictions (rank-transformed)", y = "DHARMa residual", title = title,
           subtitle = bquote(Quantile ~ test ~ italic(p) ~ .(pfmt(qp)) ~ .(if (qp < 0.05) " (significant)" else " (n.s.)"))) +
      theme_minimal(base_size = 11) +
      theme(plot.subtitle = element_text(color = if (qp < 0.05) my.cols[1] else "grey30"), panel.grid.minor = element_blank())
  }
  s5 <- cowplot::plot_grid(qq_panel(sim_all, "QQ plot (Holmes included)"),  res_panel(sim_all, "Residual vs. predicted (Holmes included)"),
                           qq_panel(sim_noH, "QQ plot (Holmes excluded)"),  res_panel(sim_noH, "Residual vs. predicted (Holmes excluded)"),
                           ncol = 2, labels = c("A", "B", "C", "D"), label_size = 12)
  s5 <- cowplot::plot_grid(cowplot::ggdraw() + cowplot::draw_label("DHARMa residual diagnostics, primary model", fontface = "bold", size = 13),
                           s5, ncol = 1, rel_heights = c(0.05, 1))
  ggsave("figures/fig_S5_dharma.pdf", s5, width = 9, height = 8, bg = "white")
} else cat("Install 'DHARMa' for simulated-residual diagnostics.\n")

## =====================================================================
## 8. Primary model + child-pathway mediation + land-use sensitivity
##     MEDIATION = child pathway removed two ways:
##       (A) per-child offset,
##       (B) kid_prop covariate.
##     SENSITIVITY = kitchen-sink full >5% set (adds Open Water, wetlands, Pasture/Hay, and
##       lumped Developed_total).
##     Confounders: rurality (RUCC) + climate (PRISM). SES as a robustness check.
## =====================================================================

## Primary confounders (rurality [RUCC, nearest-year vintage] + climate [PRISM]) are merged in
## Section 0, so `conf` already exists here. Only the robustness-only SES covariates are added now.
ses_path <- "Additional_covariate_datasets/SES/socioeconomic.csv"
if (file.exists(ses_path)) {
  se <- read.csv(ses_path)
  df_final <- merge(df_final, se[, c("FIPS", "year", "median_income", "pct_poverty", "pct_college")],
                    by = c("FIPS", "year"), all.x = TRUE)
  df_final$income10k <- df_final$median_income / 10000
}

fit_lu <- function(lu, offset_var, extra = NULL) {
  rhs <- paste(c("year_c", "amish_share", lu, extra, "(1 | county_name)",
                 paste0("offset(log(", offset_var, "))")), collapse = " + ")
  glmmTMB(as.formula(paste("LACV_cases ~", rhs)), family = nbinom2, data = df_final)
}
m_primary <- fit_lu(bt(lu_vle),  "tot_pop", conf)               # PRIMARY: total association
m_medA    <- fit_lu(bt(lu_vle),  "kid_pop", conf)               # mediation: child pool-size pathway (denominator)
m_medB    <- fit_lu(bt(lu_vle),  "tot_pop", c(conf, "kid_prop"))# mediation: child pathway (covariate)
m_fullLU  <- fit_lu(bt(lu_full), "tot_pop", conf)               # sensitivity: kitchen-sink full >5% set

cat("\n=== 8. Primary (vector larval ecology set) + mediation + land-use sensitivity ===\n")
show <- function(lbl, f) { co <- summary(f)$coefficients$cond; i <- grep("amish", rownames(co), ignore.case = TRUE)
  cat(sprintf("%-52s amish = %.2f (SE %.2f, p=%.3f)\n", lbl, co[i,1], co[i,2], co[i,4])) }
show("PRIMARY  per-capita, no kid_prop (total assoc.):", m_primary)
show("MED A    per-child offset:", m_medA)
show("MED B    per-capita + kid_prop covariate:", m_medB)
show("SENS     full >5% land-use set, per-capita:", m_fullLU)

cat("\n-- Primary model: full fixed-effects summary --\n")
print(round(summary(m_primary)$coefficients$cond, 5))

Xp <- model.matrix(m_primary)[, -1, drop = FALSE]
vifp <- sapply(seq_len(ncol(Xp)), function(j) 1 / (1 - summary(lm(Xp[, j] ~ Xp[, -j]))$r.squared))
cat("\nPrimary-model VIFs:\n"); print(round(setNames(vifp, colnames(Xp)), 2))
cat(sprintf("Amish incidence-rate ratio across observed 0-46%% Amish range (primary): %.1fx\n",
    exp(summary(m_primary)$coefficients$cond["amish_share", "Estimate"] * 0.46)))

## SES robustness: genuine SES confounders.
if (file.exists(ses_path)) {
  cat("\n-- SES robustness --\n")
  cat("SES-Amish correlations (Pearson, county-year; to justify treating SES as robustness-only):\n")
  for (v in c("median_income", "pct_poverty", "pct_college"))
    cat(sprintf("  %-14s r = %+.3f\n", v, cor(df_final[[v]], df_final$amish_share, use = "complete.obs")))
  show("+ SES confounders (income/poverty/education):", fit_lu(bt(lu_vle), "tot_pop", c(conf, "income10k", "pct_poverty", "pct_college")))
}

## ---- Fig 4: modeled incidence-rate ratio across each predictor's observed range ----
## Predictors are on incomparable scales (land-use proportions vs deg C vs mm), so instead of raw
## coefficients we plot the multiplicative change in incidence across the full observed range of
## each predictor: IRR = exp(beta * (max - min)).
co  <- summary(m_primary)$coefficients$cond
Xp  <- model.matrix(m_primary)[, -1, drop = FALSE]
rng <- apply(Xp, 2, function(z) diff(range(z))); names(rng) <- gsub("`", "", names(rng))
fe  <- data.frame(term = rownames(co), est = co[, 1], se = co[, 2])
fe  <- subset(fe, term != "(Intercept)")
fe$rng <- rng[gsub("`", "", fe$term)]
fe$irr <- exp(fe$est * fe$rng)
fe$lo  <- exp((fe$est - 1.96 * fe$se) * fe$rng); fe$hi <- exp((fe$est + 1.96 * fe$se) * fe$rng)
labmap <- c(amish_share = "Amish population share", year_c = "Year (2003-2024)",
            rucc = "Rurality (RUCC)", temp_mean = "Mean temperature",
            precip = "Precipitation", vpd = "Vapor pressure deficit")
fe$label <- ifelse(fe$term %in% names(labmap), labmap[fe$term], gsub("`|_prop", "", fe$term))
fe$hl <- fe$term == "amish_share"
fe <- fe[order(fe$irr), ]; fe$label <- factor(fe$label, levels = fe$label)
p3 <- ggplot(fe, aes(irr, label, color = hl)) +
  geom_vline(xintercept = 1, linetype = "solid", color = "grey60", linewidth = 0.4) +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.6, alpha = 0.8, na.rm = TRUE) +
  geom_point(size = 3, shape = 21, fill = "white", stroke = 1.2, na.rm = TRUE) +
  geom_point(size = 2.2, shape = 16, na.rm = TRUE) +
  scale_color_manual(values = c(`TRUE` = my.cols[1], `FALSE` = my.cols[2]), guide = "none") +
  scale_x_log10(breaks = c(0.1, 0.3, 1, 3, 10, 30, 100),
                labels = c("0.1", "0.3", "1", "3", "10", "30", "100")) +
  labs(x = "Incidence-rate ratio across observed range of predictor (log scale)", y = NULL,
       title = "Modeled effect of each predictor on LACV-ND incidence",
       subtitle = "Incidence-rate ratio across the observed range of each predictor (95% CI)") +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.major.x = element_line(color = "grey92", linewidth = 0.3),
        axis.text.y = element_text(size = 10, color = "grey20", hjust = 1),
        axis.text.x = element_text(size = 9, color = "grey30"),
        axis.title.x = element_text(face = "bold", size = 10, margin = margin(t = 10), color = "grey20"),
        axis.line.y = element_line(color = "grey40", linewidth = 0.4),
        axis.ticks.y = element_line(color = "grey40", linewidth = 0.3),
        axis.ticks.length.y = unit(0.1, "cm"),
        plot.title = element_text(face = "bold", size = 13, color = "grey20", margin = margin(b = 4)),
        plot.subtitle = element_text(size = 10, color = "grey40", margin = margin(b = 15)),
        plot.background = element_rect(fill = "white", color = NA),
        panel.background = element_rect(fill = "white", color = NA),
        plot.margin = margin(15, 15, 10, 10))
ggsave("figures/fig_4_glmm_coefficients.pdf", p3, width = 7, height = 4.5)

## =====================================================================
## 9. Spatial robustness for the primary.
##     (a) residual spatial autocorrelation under alternative weights (Moran's I)
##     (b) Conley spatial-HAC SEs -> spatially-robust inference on the estimate
##         (conleyreg's poisson is cross-sectional, so we aggregate to county).
## =====================================================================
suppressMessages(library(conleyreg))
sfc <- subset(read_sf(shapefile), STATEFP == "39")
  sfc$FIPS <- as.integer(paste0(sfc$STATEFP, sfc$COUNTYFP))
  ct  <- st_coordinates(st_centroid(st_transform(st_geometry(sfc), 4326)))
  ## (a) residual Moran's I on the PRIMARY, under queen / rook / KNN weights
  df_final$res_p <- residuals(m_primary, type = "pearson")
  rc <- aggregate(res_p ~ FIPS, df_final, mean)
  sfc$res_p <- rc$res_p[match(sfc$FIPS, rc$FIPS)]
  cat("\n=== 9. Spatial robustness ===\n(a) Residual Moran's I on the primary, by weights:\n")
  wl <- list(queen = nb2listw(poly2nb(sfc, queen = TRUE),  style = "W", zero.policy = TRUE),
             rook  = nb2listw(poly2nb(sfc, queen = FALSE), style = "W", zero.policy = TRUE),
             knn5  = nb2listw(knn2nb(knearneigh(ct, k = 5, longlat = TRUE)), style = "W", zero.policy = TRUE))
  for (nm in names(wl)) { mt <- moran.test(sfc$res_p, wl[[nm]], zero.policy = TRUE)
    cat(sprintf("    %-6s: I = %.3f, p = %.4f\n", nm, mt$estimate["Moran I statistic"], mt$p.value)) }
  ## (b) Conley spatial-HAC SEs on a county-level Poisson cross-section
  df_final$decid   <- df_final[["Deciduous Forest_prop"]]
  df_final$mixedf  <- df_final[["Mixed Forest_prop"]]
  cov_syn <- c("amish_share", "decid", "mixedf", conf)   # mirror the VLE primary (Decid + Mixed + conf)
  cty <- merge(merge(aggregate(cbind(LACV_cases, tot_pop) ~ county_name + FIPS, df_final, sum),
                     aggregate(df_final[cov_syn], list(county_name = df_final$county_name), mean)),
               data.frame(FIPS = sfc$FIPS, lon = ct[, 1], lat = ct[, 2]), by = "FIPS")
  cty$off <- log(cty$tot_pop)
  f_c <- as.formula(paste("LACV_cases ~", paste(cov_syn, collapse = " + "), "+ offset(off)"))
  a0 <- summary(glm(f_c, poisson, cty))$coefficients["amish_share", ]
  cat(sprintf("(b) Conley SEs (county-level Poisson cross-section, keeps the point estimate):\n    naive : amish = %.2f (SE %.2f, p=%.3f)\n", a0[1], a0[2], a0[4]))
  for (dc in c(100, 150, 200)) {
    ac <- conleyreg::conleyreg(f_c, data = cty, dist_cutoff = dc, model = "poisson", lat = "lat", lon = "lon")["amish_share", ]
    cat(sprintf("    %3d km: amish = %.2f (SE %.2f, p=%.3f)\n", dc, ac[1], ac[2], ac[4]))
  }
## =====================================================================

## ---- Reproducibility: session information ----
cat("\n=== sessionInfo() ===\n")
print(sessionInfo())
