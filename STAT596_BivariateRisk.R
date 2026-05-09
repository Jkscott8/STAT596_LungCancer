### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 17, 2026
### This code creates bivariate risk maps for compound risk:
### 1) Tmin + Lung Cancer
### 2) Tmax + Lung Cancer
### 3) HI + Lung Cancer
### 4) WBGT + Lung Cancer
### 5) Low Income + Lung Cancer

rm(list = ls())

## ---Packages---

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(cowplot)

## ---Paths and Settings---

setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/bivariate_risk_maps"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## ---Read and Clean Data---

df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income),
    Tmin = as.numeric(Tmin),
    Tmax = as.numeric(Tmax),   # <-- add this
    HI = as.numeric(HI),
    WBGT = as.numeric(WBGT)
  )

## ---Build County Geometry---

counties_sf <- df %>%
  select(GEOID, County, geometry) %>%
  distinct() %>%
  st_as_sf(wkt = "geometry", crs = 3857) %>%
  st_make_valid()

## ---Clean Dataframe---

county_year <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    Tmin = mean(Tmin, na.rm = TRUE),
    Tmax = mean(Tmax, na.rm = TRUE),   # <-- add this
    HI = mean(HI, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    Tmax = ifelse(is.nan(Tmax), NA, Tmax),
    HI = ifelse(is.nan(HI), NA, HI),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

## ---County Means for 2018-2023---

county_mean <- county_year %>%
  group_by(GEOID, County) %>%
  summarise(
    Lung_Cancer_Rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
    Median_Household_Income = mean(Median_Household_Income, na.rm = TRUE),
    Tmin = mean(Tmin, na.rm = TRUE),
    Tmax = mean(Tmax, na.rm = TRUE),   # <-- add this
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
    # Reverse income so that higher values = greater income-related risk
    # (i.e., lower income counties become higher-risk counties)
    Low_Income_Risk = max(Median_Household_Income, na.rm = TRUE) - Median_Household_Income
  )

## ---Bivariate Palette (3 x 3)---

bi_pal <- c(
  "1-1" = "#e8e8e8",
  "2-1" = "#dfb0d6",
  "3-1" = "#be64ac",
  "1-2" = "#ace4e4",
  "2-2" = "#a5add3",
  "3-2" = "#8c62aa",
  "1-3" = "#5ac8c8",
  "2-3" = "#5698b9",
  "3-3" = "#3b4994"
)

## ---Legend Function---

make_bi_legend <- function(xlab, ylab) {
  
  legend_df <- expand.grid(x = 1:3, y = 1:3) %>%
    mutate(bi_class = paste0(x, "-", y))
  
  ggplot(legend_df, aes(x = x, y = y, fill = bi_class)) +
    geom_tile() +
    scale_fill_manual(values = bi_pal, drop = FALSE) +
    coord_equal() +
    labs(
      x = xlab,
      y = ylab
    ) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_text(face = "bold"),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

## ---Bivariate Map Function---

make_bivariate_map <- function(sf_obj, xvar, yvar, title_text, xlab, ylab, out_stub) {
  
  sf_pair <- sf_obj %>%
    filter(!is.na(.data[[xvar]]), !is.na(.data[[yvar]])) %>%
    mutate(
      x_q = ntile(.data[[xvar]], 3),
      y_q = ntile(.data[[yvar]], 3),
      bi_class = paste0(x_q, "-", y_q),
      risk_level = case_when(
        bi_class == "3-3" ~ "Highest compound risk",
        bi_class %in% c("3-2", "2-3") ~ "High compound risk",
        bi_class %in% c("2-2") ~ "Moderate compound risk",
        TRUE ~ "Lower / mixed risk"
      )
    )
  
  p_map <- ggplot(sf_pair) +
    geom_sf(aes(fill = bi_class), color = "white", linewidth = 0.2) +
    scale_fill_manual(values = bi_pal, drop = FALSE) +
    labs(
      title = title_text,
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold")
    )
  
  p_legend <- make_bi_legend(xlab = xlab, ylab = ylab)
  
  final_plot <- ggdraw() +
    draw_plot(p_map, 0, 0, 0.82, 1) +
    draw_plot(p_legend, 0.68, 0.08, 0.22, 0.22)
  
  print(final_plot)
  
  ggsave(
    filename = file.path(output_dir, paste0(out_stub, "_bivariate_map.png")),
    plot = final_plot,
    width = 12,
    height = 9,
    dpi = 300,
    bg = "white"
  )
  
  summary_tbl <- sf_pair %>%
    st_drop_geometry() %>%
    transmute(
      pair = out_stub,
      GEOID,
      County,
      x_variable = xvar,
      x_value = .data[[xvar]],
      y_variable = yvar,
      y_value = .data[[yvar]],
      x_q,
      y_q,
      bi_class,
      risk_level
    )
  
  highest_risk_tbl <- summary_tbl %>%
    filter(bi_class == "3-3") %>%
    arrange(desc(y_value), desc(x_value))
  
  return(list(
    map_data = sf_pair,
    summary_table = summary_tbl,
    highest_risk_table = highest_risk_tbl
  ))
}

## ---Define Pairs---

pair_list <- list(
  list(
    xvar = "Tmin",
    yvar = "Lung_Cancer_Rate",
    title = "Bivariate Risk Map: Tmin and Lung Cancer Rate",
    xlab = "Higher Tmin",
    ylab = "Higher Lung Cancer",
    out_stub = "Tmin_LungCancer"
  ),
  list(
    xvar = "Tmax",
    yvar = "Lung_Cancer_Rate",
    title = "Bivariate Risk Map: Tmax and Lung Cancer Rate",
    xlab = "Higher Tmax",
    ylab = "Higher Lung Cancer",
    out_stub = "Tmax_LungCancer"
  ),
  list(
    xvar = "HI",
    yvar = "Lung_Cancer_Rate",
    title = "Bivariate Risk Map: HI and Lung Cancer Rate",
    xlab = "Higher HI",
    ylab = "Higher Lung Cancer",
    out_stub = "HI_LungCancer"
  ),
  list(
    xvar = "WBGT",
    yvar = "Lung_Cancer_Rate",
    title = "Bivariate Risk Map: WBGT and Lung Cancer Rate",
    xlab = "Higher WBGT",
    ylab = "Higher Lung Cancer",
    out_stub = "WBGT_LungCancer"
  ),
  list(
    xvar = "Low_Income_Risk",
    yvar = "Lung_Cancer_Rate",
    title = "Bivariate Risk Map: Low Income and Lung Cancer Rate",
    xlab = "Lower Income",
    ylab = "Higher Lung Cancer",
    out_stub = "LowIncome_LungCancer"
  )
)

## ---Run All Bivariate Maps---

bivariate_results <- lapply(pair_list, function(p) {
  make_bivariate_map(
    sf_obj = county_sf,
    xvar = p$xvar,
    yvar = p$yvar,
    title_text = p$title,
    xlab = p$xlab,
    ylab = p$ylab,
    out_stub = p$out_stub
  )
})

## ---Export Summary Tables---

bivariate_summary <- bind_rows(lapply(bivariate_results, function(x) x$summary_table))
bivariate_highest_risk <- bind_rows(lapply(bivariate_results, function(x) x$highest_risk_table))

write_csv(
  bivariate_summary,
  file.path(output_dir, "county_bivariate_summary_2018_2023.csv")
)

write_csv(
  bivariate_highest_risk,
  file.path(output_dir, "county_highest_compound_risk_2018_2023.csv")
)

## ---Print Highest-Risk Counties---

for (p in pair_list) {
  cat("\n====================================================\n")
  cat(p$out_stub, "\n")
  cat("\nCounties in highest compound-risk class (3-3):\n")
  
  print(
    bivariate_highest_risk %>%
      filter(pair == p$out_stub) %>%
      select(County, GEOID, x_variable, x_value, y_variable, y_value, bi_class, risk_level)
  )
}