### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 17, 2026
### This code runs local bivariate Moran's I (bivariate LISA)
### to identify statistically significant counties where:
### high lung cancer rates are associated with nearby high exposure values.
###
### Pairs analyzed:
### 1. Lung Cancer vs Tmin
### 2. Lung Cancer vs Tmax
### 3. Lung Cancer vs HI
### 4. Lung Cancer vs WBGT
### 5. Lung Cancer vs Low Income Risk

rm(list = ls())

## --- Packages ---
library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(spdep)

## --- Paths and Settings ---
setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/bivariate_lisa"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

set.seed(1234)

## --- Read and Clean Data ---
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income),
    Tmin = as.numeric(Tmin),
    Tmax = as.numeric(Tmax),
    HI = as.numeric(HI),
    WBGT = as.numeric(WBGT)
  )

## --- Build County Geometry ---
counties_sf <- df %>%
  select(GEOID, County, geometry) %>%
  distinct() %>%
  st_as_sf(wkt = "geometry", crs = 3857) %>%
  st_make_valid()

## --- Build County-Year Data ---
county_year <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    Tmin = mean(Tmin, na.rm = TRUE),
    Tmax = mean(Tmax, na.rm = TRUE),
    HI = mean(HI, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lung_Cancer_Rate = ifelse(is.nan(Lung_Cancer_Rate), NA, Lung_Cancer_Rate),
    Median_Household_Income = ifelse(is.nan(Median_Household_Income), NA, Median_Household_Income),
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    Tmax = ifelse(is.nan(Tmax), NA, Tmax),
    HI = ifelse(is.nan(HI), NA, HI),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

## --- County Means for 2018-2023 ---
county_mean <- county_year %>%
  group_by(GEOID, County) %>%
  summarise(
    Lung_Cancer_Rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
    Median_Household_Income = mean(Median_Household_Income, na.rm = TRUE),
    Tmin = mean(Tmin, na.rm = TRUE),
    Tmax = mean(Tmax, na.rm = TRUE),
    HI = mean(HI, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lung_Cancer_Rate = ifelse(is.nan(Lung_Cancer_Rate), NA, Lung_Cancer_Rate),
    Median_Household_Income = ifelse(is.nan(Median_Household_Income), NA, Median_Household_Income),
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    Tmax = ifelse(is.nan(Tmax), NA, Tmax),
    HI = ifelse(is.nan(HI), NA, HI),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

county_sf <- counties_sf %>%
  left_join(county_mean, by = c("GEOID", "County")) %>%
  mutate(
    Low_Income_Risk = max(Median_Household_Income, na.rm = TRUE) - Median_Household_Income
  )

## --- Helper: robustly find p-value column from localmoran_bv output ---
pick_p_col <- function(df_cols) {
  preferred <- c(
    "Pr(folded) Sim",
    "Pr(z != E(Ii)) Sim",
    "Pr(z != E(Ibvi)) Sim",
    "Pr(z > 0)",
    "Pr(z != E(Ii))"
  )
  
  hit <- preferred[preferred %in% df_cols]
  
  if (length(hit) > 0) {
    return(hit[1])
  }
  
  generic_hit <- grep("^Pr", df_cols, value = TRUE)
  
  if (length(generic_hit) == 0) {
    stop("No p-value column found in localmoran_bv output.")
  }
  
  generic_hit[1]
}

## --- Bivariate LISA Function ---
run_bivariate_lisa <- function(sf_obj, outcome_var, exposure_var, title_text, out_stub,
                               nsim = 999, sig_level = 0.05) {
  
  sf_pair <- sf_obj %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(.data[[exposure_var]])) %>%
    st_make_valid()
  
  ## Neighbor structure
  nb <- poly2nb(sf_pair, queen = TRUE)
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  
  ## Outcome at county i; exposure in neighboring counties
  x <- sf_pair[[outcome_var]]
  y <- sf_pair[[exposure_var]]
  
  ## Standardized values for quadrant interpretation
  x_z <- as.numeric(scale(x))
  y_z <- as.numeric(scale(y))
  lag_y_z <- lag.listw(lw, y_z, zero.policy = TRUE)
  
  ## Local bivariate Moran's I
  lisa_raw <- localmoran_bv(
    x = x,
    y = y,
    listw = lw,
    nsim = nsim,
    zero.policy = TRUE
  )
  
  lisa_df <- as.data.frame(lisa_raw)
  
  stat_col <- names(lisa_df)[1]
  p_col <- pick_p_col(names(lisa_df))
  
  z_candidates <- grep("^Z", names(lisa_df), value = TRUE)
  z_vals <- if (length(z_candidates) > 0) lisa_df[[z_candidates[1]]] else NA_real_
  
  sf_pair <- sf_pair %>%
    mutate(
      outcome_z = x_z,
      exposure_z = y_z,
      lag_exposure_z = lag_y_z,
      lisa_stat = lisa_df[[stat_col]],
      lisa_z = z_vals,
      p_raw = lisa_df[[p_col]],
      p_fdr = p.adjust(p_raw, method = "fdr"),
      significant = p_fdr <= sig_level,
      lisa_cluster = case_when(
        significant & outcome_z > 0 & lag_exposure_z > 0 ~ "High-High",
        significant & outcome_z < 0 & lag_exposure_z < 0 ~ "Low-Low",
        significant & outcome_z > 0 & lag_exposure_z < 0 ~ "High-Low",
        significant & outcome_z < 0 & lag_exposure_z > 0 ~ "Low-High",
        TRUE ~ "Not significant"
      )
    )
  
  ## Create table outputs
  result_table <- sf_pair %>%
    st_drop_geometry() %>%
    transmute(
      pair = out_stub,
      GEOID,
      County,
      outcome_variable = outcome_var,
      exposure_variable = exposure_var,
      outcome_value = .data[[outcome_var]],
      exposure_value = .data[[exposure_var]],
      outcome_z,
      lag_exposure_z,
      lisa_stat,
      lisa_z,
      p_raw,
      p_fdr,
      significant,
      lisa_cluster
    )
  
  significant_table <- result_table %>%
    filter(significant) %>%
    arrange(p_fdr, desc(lisa_stat))
  
  ## Cluster map
  p <- ggplot(sf_pair) +
    geom_sf(aes(fill = lisa_cluster), color = "white", linewidth = 0.2) +
    geom_sf(
      data = sf_pair %>% filter(significant),
      fill = NA, color = "black", linewidth = 0.35
    ) +
    scale_fill_manual(
      values = c(
        "High-High" = "#b2182b",
        "Low-Low" = "#2166ac",
        "High-Low" = "#ef8a62",
        "Low-High" = "#67a9cf",
        "Not significant" = "grey90"
      ),
      drop = FALSE,
      name = "Bivariate LISA"
    ) +
    labs(
      title = title_text,
      subtitle = paste0(
        "Significant if FDR-adjusted p ≤ ", sig_level
      ),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0(out_stub, "_bivariate_LISA_map.png")),
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
  
  return(list(
    sf = sf_pair,
    results = result_table,
    significant = significant_table
  ))
}

## --- Define Pairs ---
pair_list <- list(
  list(
    outcome_var = "Lung_Cancer_Rate",
    exposure_var = "Tmin",
    title = "Local Bivariate Moran's I: Lung Cancer vs Tmin",
    out_stub = "LungCancer_Tmin"
  ),
  list(
    outcome_var = "Lung_Cancer_Rate",
    exposure_var = "Tmax",
    title = "Local Bivariate Moran's I: Lung Cancer vs Tmax",
    out_stub = "LungCancer_Tmax"
  ),
  list(
    outcome_var = "Lung_Cancer_Rate",
    exposure_var = "HI",
    title = "Local Bivariate Moran's I: Lung Cancer vs HI",
    out_stub = "LungCancer_HI"
  ),
  list(
    outcome_var = "Lung_Cancer_Rate",
    exposure_var = "WBGT",
    title = "Local Bivariate Moran's I: Lung Cancer vs WBGT",
    out_stub = "LungCancer_WBGT"
  ),
  list(
    outcome_var = "Lung_Cancer_Rate",
    exposure_var = "Low_Income_Risk",
    title = "Local Bivariate Moran's I: Lung Cancer vs Low Income Risk",
    out_stub = "LungCancer_LowIncome"
  )
)

## --- Run All Pairs ---
all_runs <- lapply(pair_list, function(p) {
  run_bivariate_lisa(
    sf_obj = county_sf,
    outcome_var = p$outcome_var,
    exposure_var = p$exposure_var,
    title_text = p$title,
    out_stub = p$out_stub,
    nsim = 999,
    sig_level = 0.05
  )
})

## --- Export CSVs ---
all_results <- bind_rows(lapply(all_runs, function(x) x$results))
all_significant <- bind_rows(lapply(all_runs, function(x) x$significant))

write_csv(
  county_mean,
  file.path(output_dir, "county_mean_values_2018_2023.csv")
)

write_csv(
  all_results,
  file.path(output_dir, "county_bivariate_LISA_results_2018_2023.csv")
)

write_csv(
  all_significant,
  file.path(output_dir, "county_bivariate_LISA_significant_2018_2023.csv")
)

## --- Print Significant Counties ---
for (p in pair_list) {
  cat("\n====================================================\n")
  cat(p$out_stub, "\n")
  cat("\nStatistically significant counties (FDR-adjusted):\n")
  
  print(
    all_significant %>%
      filter(pair == p$out_stub) %>%
      select(
        County, GEOID, outcome_value, exposure_value,
        lisa_stat, p_raw, p_fdr, lisa_cluster
      )
  )
}