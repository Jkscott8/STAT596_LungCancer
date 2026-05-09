### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 21, 2026
### This code creates a temporal analysis using ONLY the four counties identified
### as high-risk in the bivariate LISA maps:
### 1. San Bernardino County
### 2. Riverside County
### 3. Orange County
### 4. San Diego County
###
### This code:
### - treats each county separately
### - ranks the four counties by mean household income
### - plots annual lung cancer death rate by county
### - includes income rank in the legend

rm(list = ls())

## --- Packages ---
library(dplyr)
library(readr)
library(ggplot2)
library(lmerTest)
library(emmeans)
library(stringr)

## --- Paths and Settings ---
setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/temporal_four_highrisk_counties"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

## --- Read and Clean Data ---
df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income)
  )

## --- Build County-Year Data ---
county_year <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    Lung_Cancer_Rate = ifelse(is.nan(Lung_Cancer_Rate), NA, Lung_Cancer_Rate),
    Median_Household_Income = ifelse(is.nan(Median_Household_Income), NA, Median_Household_Income)
  )

## --- Keep only the four high-risk counties ---
target_counties <- c(
  "San Bernardino",
  "Riverside",
  "Orange",
  "San Diego"
)

county_year_4 <- county_year %>%
  mutate(
    County_clean = County %>%
      str_remove(", CA$") %>%
      str_remove(" County$") %>%
      str_trim()
  ) %>%
  filter(County_clean %in% target_counties)

cat("\nCounties retained for analysis:\n")
print(unique(county_year_4$County))

## --- Rank Counties by Mean Income ---
income_ranking <- county_year_4 %>%
  group_by(GEOID, County, County_clean) %>%
  summarise(
    mean_income = mean(Median_Household_Income, na.rm = TRUE),
    mean_lung_rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    mean_income = ifelse(is.nan(mean_income), NA, mean_income),
    mean_lung_rate = ifelse(is.nan(mean_lung_rate), NA, mean_lung_rate)
  ) %>%
  arrange(desc(mean_income)) %>%
  mutate(
    income_rank = row_number(),
    income_rank_label = paste0("#", income_rank, " income")
  )

cat("\n====================================================\n")
cat("INCOME RANKING OF THE FOUR HIGH-RISK COUNTIES\n")
cat("====================================================\n")
print(income_ranking)

write_csv(
  income_ranking,
  file.path(output_dir, "income_ranking_four_highrisk_counties.csv")
)

## --- Merge Income Ranking Back Into Annual Data ---
analysis_df <- county_year_4 %>%
  left_join(
    income_ranking %>%
      select(GEOID, County, mean_income, income_rank, income_rank_label),
    by = c("GEOID", "County")
  ) %>%
  mutate(
    County_legend = paste0(County_clean, " (Income Rank ", income_rank, ")")
  )

## --- Temporal Significance Test by County ---
test_temporal_significance <- function(data) {
  
  model_df <- data %>%
    filter(!is.na(Lung_Cancer_Rate)) %>%
    mutate(
      Year_c = Year - min(Year),
      County_legend = factor(County_legend)
    )
  
  m <- lm(
    log1p(Lung_Cancer_Rate) ~ Year_c * County_legend,
    data = model_df
  )
  
  cat("\n====================================================\n")
  cat("TEMPORAL SIGNIFICANCE TEST: FOUR HIGH-RISK COUNTIES\n")
  cat("====================================================\n")
  
  cat("\nOmnibus tests:\n")
  print(anova(m))
  
  cat("\nEstimated mean differences among counties (averaged across years):\n")
  print(emmeans(m, pairwise ~ County_legend, adjust = "tukey"))
  
  cat("\nEstimated slope differences among counties:\n")
  print(emtrends(m, pairwise ~ County_legend, var = "Year_c", adjust = "tukey"))
  
  return(m)
}

sig_model <- test_temporal_significance(analysis_df)

## --- Annual County Trends ---
annual_county_trends <- analysis_df %>%
  filter(!is.na(Lung_Cancer_Rate)) %>%
  group_by(Year, County_clean, County_legend, income_rank) %>%
  summarise(
    mean_lung_cancer_rate = mean(Lung_Cancer_Rate, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(income_rank)

write_csv(
  annual_county_trends,
  file.path(output_dir, "annual_lung_cancer_trends_four_highrisk_counties.csv")
)

## --- Set Legend Order by Income Rank ---
legend_order <- income_ranking %>%
  arrange(income_rank) %>%
  mutate(
    County_legend = paste0(County_clean, " (Income Rank ", income_rank, ")")
  ) %>%
  pull(County_legend)

annual_county_trends <- annual_county_trends %>%
  mutate(
    County_legend = factor(County_legend, levels = legend_order)
  )

## --- Plot ---
p <- ggplot(
  annual_county_trends,
  aes(
    x = Year,
    y = mean_lung_cancer_rate,
    color = County_legend,
    group = County_legend
  )
) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.8) +
  labs(
    title = "Annual Lung Cancer Death Rate in High-Risk Counties",
    x = "Year",
    y = "Lung Cancer Death Rate (per 100,000 people)",
    color = "County"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5),
    legend.title = element_text(face = "bold")
  )

ggsave(
  filename = file.path(output_dir, "annual_lung_cancer_trends_four_highrisk_counties.png"),
  plot = p,
  width = 10,
  height = 7,
  dpi = 300,
  bg = "white"
)

print(p)
