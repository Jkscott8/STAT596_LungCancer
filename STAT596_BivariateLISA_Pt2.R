### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 18, 2026
### This code identifies which counties fall into the same-county
### WBGT + Income and Tmin + Income groups shown in the annual figures.
###
### Group definitions:
### - High heat = county value >= statewide median for that year
### - Low heat  = county value < statewide median for that year
### - Low income = county income < statewide median for that year
### - High income = county income >= statewide median for that year

rm(list = ls())

## --- Packages ---
library(sf)
library(dplyr)
library(readr)
library(ggplot2)

## --- Paths and Settings ---
setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/heat_income_groups"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "maps_wbgt_income"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(output_dir, "maps_tmin_income"), showWarnings = FALSE, recursive = TRUE)

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

## --- Helper for safe means ---
safe_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  } else {
    return(mean(x, na.rm = TRUE))
  }
}

## --- Build County-Year Data ---
county_year <- df %>%
  group_by(GEOID, County, Year) %>%
  summarise(
    Lung_Cancer_Rate = first(na.omit(Lung_Cancer_Rate), default = NA_real_),
    Median_Household_Income = first(na.omit(Median_Household_Income), default = NA_real_),
    Tmin = safe_mean(Tmin),
    Tmax = safe_mean(Tmax),
    HI = safe_mean(HI),
    WBGT = safe_mean(WBGT),
    .groups = "drop"
  )

## --- County Means for 2018-2023 ---
county_mean <- county_year %>%
  group_by(GEOID, County) %>%
  summarise(
    Lung_Cancer_Rate = safe_mean(Lung_Cancer_Rate),
    Median_Household_Income = safe_mean(Median_Household_Income),
    Tmin = safe_mean(Tmin),
    Tmax = safe_mean(Tmax),
    HI = safe_mean(HI),
    WBGT = safe_mean(WBGT),
    .groups = "drop"
  )

## --- Group labels and colors ---
wbgt_levels <- c(
  "High WBGT + Low Income",
  "High WBGT + High Income",
  "Low WBGT + Low Income",
  "Low WBGT + High Income"
)

tmin_levels <- c(
  "High Tmin + Low Income",
  "High Tmin + High Income",
  "Low Tmin + Low Income",
  "Low Tmin + High Income"
)

group_colors <- c(
  "High WBGT + Low Income"  = "#b2182b",
  "High WBGT + High Income" = "#ef8a62",
  "Low WBGT + Low Income"   = "#67a9cf",
  "Low WBGT + High Income"  = "#2166ac",
  "High Tmin + Low Income"  = "#b2182b",
  "High Tmin + High Income" = "#ef8a62",
  "Low Tmin + Low Income"   = "#67a9cf",
  "Low Tmin + High Income"  = "#2166ac"
)

## --- Annual county-year group assignment: WBGT + Income ---
wbgt_income_year <- county_year %>%
  group_by(Year) %>%
  mutate(
    wbgt_year_median = median(WBGT, na.rm = TRUE),
    income_year_median = median(Median_Household_Income, na.rm = TRUE),
    wbgt_level = case_when(
      is.na(WBGT) ~ NA_character_,
      WBGT >= wbgt_year_median ~ "High",
      WBGT < wbgt_year_median ~ "Low"
    ),
    income_level = case_when(
      is.na(Median_Household_Income) ~ NA_character_,
      Median_Household_Income < income_year_median ~ "Low",
      Median_Household_Income >= income_year_median ~ "High"
    ),
    wbgt_income_group = case_when(
      wbgt_level == "High" & income_level == "Low" ~ "High WBGT + Low Income",
      wbgt_level == "High" & income_level == "High" ~ "High WBGT + High Income",
      wbgt_level == "Low"  & income_level == "Low" ~ "Low WBGT + Low Income",
      wbgt_level == "Low"  & income_level == "High" ~ "Low WBGT + High Income",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    wbgt_income_group = factor(wbgt_income_group, levels = wbgt_levels)
  )

## --- Annual county-year group assignment: Tmin + Income ---
tmin_income_year <- county_year %>%
  group_by(Year) %>%
  mutate(
    tmin_year_median = median(Tmin, na.rm = TRUE),
    income_year_median = median(Median_Household_Income, na.rm = TRUE),
    tmin_level = case_when(
      is.na(Tmin) ~ NA_character_,
      Tmin >= tmin_year_median ~ "High",
      Tmin < tmin_year_median ~ "Low"
    ),
    income_level = case_when(
      is.na(Median_Household_Income) ~ NA_character_,
      Median_Household_Income < income_year_median ~ "Low",
      Median_Household_Income >= income_year_median ~ "High"
    ),
    tmin_income_group = case_when(
      tmin_level == "High" & income_level == "Low" ~ "High Tmin + Low Income",
      tmin_level == "High" & income_level == "High" ~ "High Tmin + High Income",
      tmin_level == "Low"  & income_level == "Low" ~ "Low Tmin + Low Income",
      tmin_level == "Low"  & income_level == "High" ~ "Low Tmin + High Income",
      TRUE ~ NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    tmin_income_group = factor(tmin_income_group, levels = tmin_levels)
  )

## --- 2018-2023 overall mean group assignment ---
overall_income_median <- median(county_mean$Median_Household_Income, na.rm = TRUE)
overall_wbgt_median   <- median(county_mean$WBGT, na.rm = TRUE)
overall_tmin_median   <- median(county_mean$Tmin, na.rm = TRUE)

county_mean_groups <- county_mean %>%
  mutate(
    wbgt_level_mean = case_when(
      is.na(WBGT) ~ NA_character_,
      WBGT >= overall_wbgt_median ~ "High",
      WBGT < overall_wbgt_median ~ "Low"
    ),
    tmin_level_mean = case_when(
      is.na(Tmin) ~ NA_character_,
      Tmin >= overall_tmin_median ~ "High",
      Tmin < overall_tmin_median ~ "Low"
    ),
    income_level_mean = case_when(
      is.na(Median_Household_Income) ~ NA_character_,
      Median_Household_Income < overall_income_median ~ "Low",
      Median_Household_Income >= overall_income_median ~ "High"
    ),
    wbgt_income_group_2018_2023 = case_when(
      wbgt_level_mean == "High" & income_level_mean == "Low" ~ "High WBGT + Low Income",
      wbgt_level_mean == "High" & income_level_mean == "High" ~ "High WBGT + High Income",
      wbgt_level_mean == "Low"  & income_level_mean == "Low" ~ "Low WBGT + Low Income",
      wbgt_level_mean == "Low"  & income_level_mean == "High" ~ "Low WBGT + High Income",
      TRUE ~ NA_character_
    ),
    tmin_income_group_2018_2023 = case_when(
      tmin_level_mean == "High" & income_level_mean == "Low" ~ "High Tmin + Low Income",
      tmin_level_mean == "High" & income_level_mean == "High" ~ "High Tmin + High Income",
      tmin_level_mean == "Low"  & income_level_mean == "Low" ~ "Low Tmin + Low Income",
      tmin_level_mean == "Low"  & income_level_mean == "High" ~ "Low Tmin + High Income",
      TRUE ~ NA_character_
    )
  ) %>%
  mutate(
    wbgt_income_group_2018_2023 = factor(wbgt_income_group_2018_2023, levels = wbgt_levels),
    tmin_income_group_2018_2023 = factor(tmin_income_group_2018_2023, levels = tmin_levels)
  )

## --- Annual summaries matching the line figures ---
annual_wbgt_summary <- wbgt_income_year %>%
  filter(!is.na(wbgt_income_group)) %>%
  group_by(Year, wbgt_income_group) %>%
  summarise(
    n_counties = n(),
    mean_lung_cancer_rate = safe_mean(Lung_Cancer_Rate),
    mean_wbgt = safe_mean(WBGT),
    mean_income = safe_mean(Median_Household_Income),
    .groups = "drop"
  )

annual_tmin_summary <- tmin_income_year %>%
  filter(!is.na(tmin_income_group)) %>%
  group_by(Year, tmin_income_group) %>%
  summarise(
    n_counties = n(),
    mean_lung_cancer_rate = safe_mean(Lung_Cancer_Rate),
    mean_tmin = safe_mean(Tmin),
    mean_income = safe_mean(Median_Household_Income),
    .groups = "drop"
  )

## --- County counts by group ---
wbgt_group_counts <- wbgt_income_year %>%
  filter(!is.na(wbgt_income_group)) %>%
  count(Year, wbgt_income_group, name = "n_counties")

tmin_group_counts <- tmin_income_year %>%
  filter(!is.na(tmin_income_group)) %>%
  count(Year, tmin_income_group, name = "n_counties")

## --- Join geometry for annual maps ---
wbgt_map_sf <- counties_sf %>%
  left_join(
    wbgt_income_year %>%
      select(
        GEOID, County, Year, Lung_Cancer_Rate, Median_Household_Income,
        WBGT, wbgt_year_median, income_year_median, wbgt_income_group
      ),
    by = c("GEOID", "County")
  )

tmin_map_sf <- counties_sf %>%
  left_join(
    tmin_income_year %>%
      select(
        GEOID, County, Year, Lung_Cancer_Rate, Median_Household_Income,
        Tmin, tmin_year_median, income_year_median, tmin_income_group
      ),
    by = c("GEOID", "County")
  )

county_mean_map_sf <- counties_sf %>%
  left_join(county_mean_groups, by = c("GEOID", "County"))

## --- Helper plotting function ---
plot_group_map <- function(sf_obj, group_var, title_text, subtitle_text,
                           legend_title, out_file, color_values) {
  
  p <- ggplot(sf_obj) +
    geom_sf(aes(fill = .data[[group_var]]), color = "white", linewidth = 0.15) +
    scale_fill_manual(
      values = color_values,
      drop = FALSE,
      na.value = "grey90",
      name = legend_title
    ) +
    labs(
      title = title_text,
      subtitle = subtitle_text,
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
    filename = out_file,
    plot = p,
    width = 10,
    height = 8,
    dpi = 300,
    bg = "white"
  )
}

## --- Annual WBGT + Income maps ---
years <- sort(unique(county_year$Year))

for (yr in years) {
  sf_year <- wbgt_map_sf %>%
    filter(Year == yr)
  
  plot_group_map(
    sf_obj = sf_year,
    group_var = "wbgt_income_group",
    title_text = paste("County Groups:", yr, "WBGT + Income"),
    subtitle_text = "High/Low based on yearly statewide medians",
    legend_title = "County Group",
    out_file = file.path(
      output_dir, "maps_wbgt_income",
      paste0("WBGT_income_groups_", yr, ".png")
    ),
    color_values = group_colors[wbgt_levels]
  )
}

## --- Annual Tmin + Income maps ---
for (yr in years) {
  sf_year <- tmin_map_sf %>%
    filter(Year == yr)
  
  plot_group_map(
    sf_obj = sf_year,
    group_var = "tmin_income_group",
    title_text = paste("County Groups:", yr, "Tmin + Income"),
    subtitle_text = "High/Low based on yearly statewide medians",
    legend_title = "County Group",
    out_file = file.path(
      output_dir, "maps_tmin_income",
      paste0("Tmin_income_groups_", yr, ".png")
    ),
    color_values = group_colors[tmin_levels]
  )
}

## --- Overall 2018-2023 mean maps ---
plot_group_map(
  sf_obj = county_mean_map_sf,
  group_var = "wbgt_income_group_2018_2023",
  title_text = "County Groups: 2018-2023 Mean WBGT + Income",
  subtitle_text = "High/Low based on 2018-2023 county mean medians",
  legend_title = "County Group",
  out_file = file.path(output_dir, "WBGT_income_groups_2018_2023_mean.png"),
  color_values = group_colors[wbgt_levels]
)

plot_group_map(
  sf_obj = county_mean_map_sf,
  group_var = "tmin_income_group_2018_2023",
  title_text = "County Groups: 2018-2023 Mean Tmin + Income",
  subtitle_text = "High/Low based on 2018-2023 county mean medians",
  legend_title = "County Group",
  out_file = file.path(output_dir, "Tmin_income_groups_2018_2023_mean.png"),
  color_values = group_colors[tmin_levels]
)

## --- Export CSVs: county-year membership tables ---
write_csv(
  wbgt_income_year %>%
    select(
      GEOID, County, Year,
      Lung_Cancer_Rate,
      Median_Household_Income, income_year_median,
      WBGT, wbgt_year_median,
      wbgt_level, income_level, wbgt_income_group
    ) %>%
    arrange(Year, wbgt_income_group, County),
  file.path(output_dir, "county_year_WBGT_income_group_membership.csv")
)

write_csv(
  tmin_income_year %>%
    select(
      GEOID, County, Year,
      Lung_Cancer_Rate,
      Median_Household_Income, income_year_median,
      Tmin, tmin_year_median,
      tmin_level, income_level, tmin_income_group
    ) %>%
    arrange(Year, tmin_income_group, County),
  file.path(output_dir, "county_year_Tmin_income_group_membership.csv")
)

## --- Export CSVs: overall 2018-2023 mean county groups ---
write_csv(
  county_mean_groups %>%
    select(
      GEOID, County,
      Lung_Cancer_Rate,
      Median_Household_Income,
      Tmin,
      WBGT,
      wbgt_income_group_2018_2023,
      tmin_income_group_2018_2023
    ) %>%
    arrange(wbgt_income_group_2018_2023, tmin_income_group_2018_2023, County),
  file.path(output_dir, "county_mean_group_membership_2018_2023.csv")
)

## --- Export CSVs: annual summaries used for the figures ---
write_csv(
  annual_wbgt_summary,
  file.path(output_dir, "annual_WBGT_income_group_summary.csv")
)

write_csv(
  annual_tmin_summary,
  file.path(output_dir, "annual_Tmin_income_group_summary.csv")
)

write_csv(
  wbgt_group_counts,
  file.path(output_dir, "annual_WBGT_income_group_counts.csv")
)

write_csv(
  tmin_group_counts,
  file.path(output_dir, "annual_Tmin_income_group_counts.csv")
)

write_csv(
  county_year,
  file.path(output_dir, "county_year_values_2018_2023.csv")
)

write_csv(
  county_mean,
  file.path(output_dir, "county_mean_values_2018_2023.csv")
)

## --- Print county membership to console ---
for (yr in years) {
  cat("\n====================================================\n")
  cat("WBGT + Income Groups | Year:", yr, "\n")
  cat("====================================================\n")
  
  print(
    wbgt_income_year %>%
      filter(Year == yr) %>%
      select(
        County, GEOID, WBGT, wbgt_year_median,
        Median_Household_Income, income_year_median,
        wbgt_income_group
      ) %>%
      arrange(wbgt_income_group, County)
  )
}

for (yr in years) {
  cat("\n====================================================\n")
  cat("Tmin + Income Groups | Year:", yr, "\n")
  cat("====================================================\n")
  
  print(
    tmin_income_year %>%
      filter(Year == yr) %>%
      select(
        County, GEOID, Tmin, tmin_year_median,
        Median_Household_Income, income_year_median,
        tmin_income_group
      ) %>%
      arrange(tmin_income_group, County)
  )
}

cat("\n====================================================\n")
cat("Done. Outputs saved to:\n")
cat(output_dir, "\n")
cat("====================================================\n")