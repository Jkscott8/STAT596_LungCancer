### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 18, 2026
### This code performs a targeted temporal analysis for overlap counties:
### counties that are both
### 1. High-High in the bivariate LISA maps, and
### 2. High Heat + High Income in the county-group maps

rm(list = ls())

## --- Packages ---
library(dplyr)
library(readr)
library(ggplot2)
library(ggrepel)

## --- Paths and Settings ---
setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/overlap_temporal_analysis"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## --- Hard-code overlap counties from your maps ---
## Use exact names as they appear in df$County
## Example format: "Alameda County, CA"

tmin_overlap_counties <- c(
   "San Diego County, CA",
   "Orange County, CA",
   "Riverside County, CA",
   "San Bernardino County, CA"
)

wbgt_overlap_counties <- c(
   "San Diego County, CA",
   "Riverside County, CA"
)

## --- Read and Clean Data ---
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Month_num = as.integer(Month_num),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income),
    Tmin = as.numeric(Tmin),
    WBGT = as.numeric(WBGT)
  ) %>%
  filter(Year >= 2018, Year <= 2023)

## --- Add season labels ---
## December is assigned to the following winter year
## so DJF 2018 will drop because Dec 2017 is not in the dataset
monthly <- df %>%
  mutate(
    season = case_when(
      Month_num %in% c(1, 2, 3)   ~ "JFM",
      Month_num %in% c(4, 5, 6)   ~ "AMJ",
      Month_num %in% c(7, 8, 9)   ~ "JAS",
      Month_num %in% c(10, 11, 12) ~ "OND"
    ),
    season_year = Year,
    season = factor(season, levels = c("JFM", "AMJ", "JAS", "OND"))
  ) %>%
  filter(season_year >= 2018, season_year <= 2023)

## --- Build county-season-year table ---
seasonal_df <- monthly %>%
  group_by(GEOID, County, season_year, season) %>%
  summarise(
    n_months = n_distinct(Month_num),
    Tmin = mean(Tmin, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    .groups = "drop"
  ) %>%
  filter(n_months == 3) %>%
  rename(Year = season_year) %>%
  mutate(
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

## --- Annual lung cancer table ---
annual_lung_df <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    .groups = "drop"
  )

## --- Helper function ---
run_overlap_temporal_analysis <- function(seasonal_data, annual_data, county_vec,
                                          heat_var, heat_label, out_stub) {
  
  if (length(county_vec) == 0) {
    cat("\n====================================================\n")
    cat("Skipping", out_stub, "- no counties entered.\n")
    return(NULL)
  }
  
  overlap_seasonal <- seasonal_data %>%
    filter(County %in% county_vec) %>%
    mutate(
      heat_value = .data[[heat_var]],
      Year_c = Year - 2018
    )
  
  overlap_annual_lung <- annual_data %>%
    filter(County %in% county_vec)
  
  ## Mean seasonal trend across overlap counties
  seasonal_mean <- overlap_seasonal %>%
    group_by(Year, season) %>%
    summarise(
      mean_heat = mean(heat_value, na.rm = TRUE),
      .groups = "drop"
    )
  
  p1 <- ggplot(
    seasonal_mean,
    aes(x = Year, y = mean_heat, color = season, group = season)
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_x_continuous(breaks = 2018:2023) +
    labs(
      title = paste0("Seasonal Trends: Overlap Counties (", heat_label, ")"),
      x = "Year",
      y = paste0("Mean seasonal ", heat_label),
      color = "Season"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  ggsave(
    file.path(output_dir, paste0(out_stub, "_seasonal_mean_trend.png")),
    p1, width = 10, height = 7, dpi = 300, bg = "white"
  )
  
  ## County-level seasonal trend
  p2 <- ggplot(
    overlap_seasonal,
    aes(x = Year, y = heat_value, color = season, group = season)
  ) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    facet_wrap(~ County) +
    scale_x_continuous(breaks = 2018:2023) +
    labs(
      title = paste0("County-Level Seasonal Trends: ", heat_label),
      x = "Year",
      y = paste0("Mean seasonal ", heat_label),
      color = "Season"
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  ggsave(
    file.path(output_dir, paste0(out_stub, "_county_seasonal_trends.png")),
    p2, width = 12, height = 8, dpi = 300, bg = "white"
  )
  
  ## Seasonal anomalies:
  ## anomaly = county-season-year value minus that county's 2018-2023 mean for that same season
  anomaly_df <- overlap_seasonal %>%
    group_by(County, season) %>%
    mutate(
      seasonal_baseline = mean(heat_value, na.rm = TRUE),
      anomaly = heat_value - seasonal_baseline
    ) %>%
    ungroup()
  
  anomaly_mean <- anomaly_df %>%
    group_by(Year, season) %>%
    summarise(
      mean_anomaly = mean(anomaly, na.rm = TRUE),
      .groups = "drop"
    )
  
  p3 <- ggplot(anomaly_df, aes(x = Year, y = anomaly, group = County)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    geom_line(alpha = 0.30, color = "gray65") +
    geom_point(alpha = 0.30, color = "gray65") +
    geom_line(
      data = anomaly_mean,
      aes(x = Year, y = mean_anomaly, group = 1),
      color = "#b2182b",
      linewidth = 1.2,
      inherit.aes = FALSE
    ) +
    geom_point(
      data = anomaly_mean,
      aes(x = Year, y = mean_anomaly),
      color = "#b2182b",
      size = 2.2,
      inherit.aes = FALSE
    ) +
    facet_wrap(~ season, nrow = 1) +
    scale_x_continuous(breaks = 2018:2023) +
    labs(
      title = paste0("Seasonal Anomalies: Overlap Counties (", heat_label, ")"),
      subtitle = "Gray lines = individual counties; Red line = mean",
      x = "Year",
      y = paste0(heat_label, " anomaly")
    ) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      strip.text = element_text(face = "bold")
    )
  
  ggsave(
    file.path(output_dir, paste0(out_stub, "_seasonal_anomalies.png")),
    p3, width = 12, height = 4.5, dpi = 300, bg = "white"
  )
  
  ## Annual lung cancer trend
  lung_mean <- overlap_annual_lung %>%
    group_by(Year) %>%
    summarise(
      Lung_Cancer_Rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
      County = "Mean",
      .groups = "drop"
    )
  
  lung_plot_df <- bind_rows(
    overlap_annual_lung %>%
      select(Year, Lung_Cancer_Rate, County),
    lung_mean
  )
  
  p4 <- ggplot(
    lung_plot_df,
    aes(x = Year, y = Lung_Cancer_Rate, color = County, group = County)
  ) +
    geom_line(
      aes(linewidth = County),
      alpha = 0.95
    ) +
    geom_point(size = 2.6, alpha = 0.95) +
    scale_color_manual(
      values = c(
        "Orange County, CA" = "#1b9e77",
        "Riverside County, CA" = "#d95f02",
        "San Bernardino County, CA" = "#7570b3",
        "San Diego County, CA" = "#e7298a",
        "Mean" = "black"
      )
    ) +
    scale_linewidth_manual(
      values = c(
        "Orange County, CA" = 1.1,
        "Riverside County, CA" = 1.1,
        "San Bernardino County, CA" = 1.1,
        "San Diego County, CA" = 1.1,
        "Mean" = 1.8
      )
    ) +
    scale_x_continuous(breaks = 2018:2023) +
    labs(
      title = paste0("Annual Lung Cancer Trends: ", heat_label, " Overlap Counties"),
      x = "Year",
      y = "Lung cancer mortality rate (per 100,000)",
      color = "Series",
      linewidth = "Series"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  ggsave(
    file.path(output_dir, paste0(out_stub, "_annual_lung_cancer_trend.png")),
    p4, width = 10, height = 7, dpi = 300, bg = "white"
  )
  
  ## Simple fixed-effects trend model:
  ## tests whether slope differs by season after accounting for county baseline differences
  trend_model <- lm(
    heat_value ~ Year_c * season + County,
    data = overlap_seasonal
  )
  
  coef_table <- as.data.frame(summary(trend_model)$coefficients)
  coef_table$term <- rownames(coef_table)
  rownames(coef_table) <- NULL
  coef_table <- coef_table %>%
    select(term, everything())
  
  write_csv(
    overlap_seasonal,
    file.path(output_dir, paste0(out_stub, "_seasonal_values.csv"))
  )
  
  write_csv(
    anomaly_df,
    file.path(output_dir, paste0(out_stub, "_seasonal_anomalies.csv"))
  )
  
  write_csv(
    overlap_annual_lung,
    file.path(output_dir, paste0(out_stub, "_annual_lung_cancer_values.csv"))
  )
  
  write_csv(
    coef_table,
    file.path(output_dir, paste0(out_stub, "_trend_model_coefficients.csv"))
  )
  
  cat("\n====================================================\n")
  cat(out_stub, "\n")
  cat("Overlap counties:\n")
  print(county_vec)
  cat("\nTrend model summary:\n")
  print(summary(trend_model))
  
  return(list(
    seasonal_values = overlap_seasonal,
    anomalies = anomaly_df,
    annual_lung = overlap_annual_lung,
    model = trend_model
  ))
}

## --- Run analyses ---
tmin_results <- run_overlap_temporal_analysis(
  seasonal_data = seasonal_df,
  annual_data = annual_lung_df,
  county_vec = tmin_overlap_counties,
  heat_var = "Tmin",
  heat_label = "Tmin",
  out_stub = "Tmin_overlap"
)

wbgt_results <- run_overlap_temporal_analysis(
  seasonal_data = seasonal_df,
  annual_data = annual_lung_df,
  county_vec = wbgt_overlap_counties,
  heat_var = "WBGT",
  heat_label = "WBGT",
  out_stub = "WBGT_overlap"
)