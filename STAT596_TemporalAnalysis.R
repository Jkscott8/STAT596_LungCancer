### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 18, 2026
### This code creates a simple temporal analysis by grouping California counties into:
### 1. High heat + Low Income
### 2. High heat + High Income
### 3. Low heat + Low Income
### 4. Low heat + High Income
###
### It does this for BOTH:
### - WBGT
### - Tmin
###
### Then it plots the annual mean lung cancer death rate for each group.


rm(list = ls())

## --- Packages ---
library(dplyr)
library(readr)
library(ggplot2)
library(lmerTest)
library(emmeans)

## --- Paths and Settings ---
setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/temporal_group_trends"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## --- Read and Clean Data ---
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income),
    Tmin = as.numeric(Tmin),
    WBGT = as.numeric(WBGT)
  )

## --- Build County-Year Data ---
county_year <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    Tmin = mean(Tmin, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lung_Cancer_Rate = ifelse(is.nan(Lung_Cancer_Rate), NA, Lung_Cancer_Rate),
    Median_Household_Income = ifelse(is.nan(Median_Household_Income), NA, Median_Household_Income),
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

test_temporal_significance <- function(grouped_data, group_var, label_text) {
  
  model_df <- grouped_data %>%
    filter(
      !is.na(Lung_Cancer_Rate),
      !is.na(.data[[group_var]])
    ) %>%
    mutate(
      Year_c = Year - min(Year),   # center year so intercept is easier
      group = factor(.data[[group_var]])
    )
  
  ## Mixed model: counties repeated over time
  m <- lmer(
    log1p(Lung_Cancer_Rate) ~ Year_c * group + (1 | GEOID),
    data = model_df
  )
  
  cat("\n====================================================\n")
  cat("TEMPORAL SIGNIFICANCE TEST:", label_text, "\n")
  cat("====================================================\n")
  
  cat("\nOmnibus tests:\n")
  print(anova(m, type = 3))
  
  cat("\nEstimated mean differences among groups (averaged across years):\n")
  print(emmeans(m, pairwise ~ group, adjust = "tukey"))
  
  cat("\nEstimated slope differences among groups:\n")
  print(emtrends(m, pairwise ~ group, var = "Year_c", adjust = "tukey"))
  
  return(m)
}

## --- Function to Run Grouped Temporal Analysis ---
run_temporal_groups <- function(data, heat_var, heat_label, output_dir) {
  
  county_mean <- data %>%
    group_by(GEOID, County) %>%
    summarise(
      mean_lung = mean(Lung_Cancer_Rate, na.rm = TRUE),
      mean_income = mean(Median_Household_Income, na.rm = TRUE),
      mean_heat = mean(.data[[heat_var]], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      mean_lung = ifelse(is.nan(mean_lung), NA, mean_lung),
      mean_income = ifelse(is.nan(mean_income), NA, mean_income),
      mean_heat = ifelse(is.nan(mean_heat), NA, mean_heat)
    )
  
  heat_cutoff <- median(county_mean$mean_heat, na.rm = TRUE)
  income_cutoff <- median(county_mean$mean_income, na.rm = TRUE)
  
  county_groups <- county_mean %>%
    mutate(
      heat_group = ifelse(mean_heat >= heat_cutoff,
                          paste("High", heat_label),
                          paste("Low", heat_label)),
      income_group = ifelse(mean_income < income_cutoff,
                            "Low Income",
                            "High Income"),
      heat_income_group = case_when(
        heat_group == paste("High", heat_label) & income_group == "Low Income"  ~ paste("High", heat_label, "+ Low Income"),
        heat_group == paste("High", heat_label) & income_group == "High Income" ~ paste("High", heat_label, "+ High Income"),
        heat_group == paste("Low", heat_label)  & income_group == "Low Income"  ~ paste("Low", heat_label, "+ Low Income"),
        heat_group == paste("Low", heat_label)  & income_group == "High Income" ~ paste("Low", heat_label, "+ High Income")
      )
    )
  
  grouped_data <- data %>%
    left_join(
      county_groups %>%
        select(GEOID, County, heat_group, income_group, heat_income_group),
      by = c("GEOID", "County")
    )
  
  sig_model <- test_temporal_significance(
    grouped_data = grouped_data,
    group_var = "heat_income_group",
    label_text = heat_label
  )
  
  group_levels <- c(
    paste("High", heat_label, "+ Low Income"),
    paste("High", heat_label, "+ High Income"),
    paste("Low", heat_label, "+ Low Income"),
    paste("Low", heat_label, "+ High Income")
  )
  
  annual_group_trends <- grouped_data %>%
    filter(!is.na(Lung_Cancer_Rate), !is.na(heat_income_group)) %>%
    group_by(Year, heat_income_group) %>%
    summarise(
      mean_lung_cancer_rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
      n_counties = n(),
      .groups = "drop"
    ) %>%
    mutate(
      heat_income_group = factor(heat_income_group, levels = group_levels)
    )
  
  group_colors <- c("#b2182b", "#ef8a62", "#67a9cf", "#2166ac")
  names(group_colors) <- group_levels
  
  p <- ggplot(
    annual_group_trends,
    aes(
      x = Year,
      y = mean_lung_cancer_rate,
      color = heat_income_group,
      group = heat_income_group
    )
  ) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2.5) +
    scale_color_manual(
      values = group_colors,
      name = "County Group"
    ) +
    labs(
      title = paste("Annual Mean Lung Cancer Rate by", heat_label, "+ Income Group"),
      x = "Year",
      y = "Mean Lung Cancer Death Rate"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5),
      legend.title = element_text(face = "bold")
    )
  
  ggsave(
    filename = file.path(output_dir, paste0("annual_lung_cancer_trends_by_", heat_var, "_income_group.png")),
    plot = p,
    width = 10,
    height = 7,
    dpi = 300,
    bg = "white"
  )
  
  return(list(
    county_groups = county_groups,
    annual_group_trends = annual_group_trends,
    significance_model = sig_model,
    plot = p
  ))
}

## --- Run for WBGT ---
wbgt_results <- run_temporal_groups(
  data = county_year,
  heat_var = "WBGT",
  heat_label = "WBGT",
  output_dir = output_dir
)

## --- Run for Tmin ---
tmin_results <- run_temporal_groups(
  data = county_year,
  heat_var = "Tmin",
  heat_label = "Tmin",
  output_dir = output_dir
)