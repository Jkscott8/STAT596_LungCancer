### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 17, 2026
### This code conducts a hotspot analysis in R to visualize any statistically 
### significant hotspots where lung cancer death rates are the highest, median 
### household income is the lowest, and which counties are the hottest in terms 
### of each CHIRTS-ERA5 variable of interest (TOTAL = 6 maps).


rm(list = ls())

## ---Packages---

library(sf)
library(dplyr)
library(readr)
library(ggplot2)
library(spdep)

## ---Paths and Settings---

setwd("E:/STAT596_SP26")

input_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"
output_dir <- "E:/STAT596_SP26/hotspots"

## ---Read and Clean Data---

df <- read_csv(input_file, show_col_types = FALSE) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Lung_Cancer_Rate = ifelse(Lung_Cancer_Rate == "Suppressed", NA, Lung_Cancer_Rate),
    Lung_Cancer_Rate = as.numeric(Lung_Cancer_Rate),
    Median_Household_Income = as.numeric(Median_Household_Income),
    Tmin = as.numeric(Tmin),
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
    HI = mean(HI, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
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
    HI = mean(HI, na.rm = TRUE),
    WBGT = mean(WBGT, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Lung_Cancer_Rate = ifelse(is.nan(Lung_Cancer_Rate), NA, Lung_Cancer_Rate),
    Median_Household_Income = ifelse(is.nan(Median_Household_Income), NA, Median_Household_Income),
    Tmin = ifelse(is.nan(Tmin), NA, Tmin),
    HI = ifelse(is.nan(HI), NA, HI),
    WBGT = ifelse(is.nan(WBGT), NA, WBGT)
  )

county_sf <- counties_sf %>%
  left_join(county_mean, by = c("GEOID", "County"))

## ---Variables---

vars <- c(
  "Lung_Cancer_Rate",
  "Median_Household_Income",
  "Tmin",
  "HI",
  "WBGT"
)

hotspot_tables <- list()
topbottom_tables <- list()

## ---Hotspot Analysis & Maps

for (v in vars) {
  
  sf_v <- county_sf %>%
    filter(!is.na(.data[[v]]))
  
  nb <- poly2nb(sf_v, queen = TRUE)
  lw <- nb2listw(include.self(nb), style = "B", zero.policy = TRUE)
  
  sf_v <- sf_v %>%
    mutate(
      gi_z = as.numeric(localG(.data[[v]], lw, zero.policy = TRUE)),
      p_value = 2 * pnorm(abs(gi_z), lower.tail = FALSE),
      hotspot_class = case_when(
        gi_z >=  2.58 & p_value <= 0.01 ~ "Hotspot (99%)",
        gi_z >=  1.96 & p_value <= 0.05 ~ "Hotspot (95%)",
        gi_z <= -2.58 & p_value <= 0.01 ~ "Coldspot (99%)",
        gi_z <= -1.96 & p_value <= 0.05 ~ "Coldspot (95%)",
        TRUE ~ "Not significant"
      )
    )
  
  hotspot_tables[[v]] <- sf_v %>%
    st_drop_geometry() %>%
    select(GEOID, County, all_of(v), gi_z, p_value, hotspot_class) %>%
    rename(value = all_of(v)) %>%
    mutate(variable = v) %>%
    select(variable, everything())
  
  topbottom_tables[[v]] <- bind_rows(
    sf_v %>%
      st_drop_geometry() %>%
      arrange(desc(.data[[v]])) %>%
      slice_head(n = 10) %>%
      select(GEOID, County, all_of(v), gi_z, p_value, hotspot_class) %>%
      rename(value = all_of(v)) %>%
      mutate(variable = v, extreme = "Highest"),
    
    sf_v %>%
      st_drop_geometry() %>%
      arrange(.data[[v]]) %>%
      slice_head(n = 10) %>%
      select(GEOID, County, all_of(v), gi_z, p_value, hotspot_class) %>%
      rename(value = all_of(v)) %>%
      mutate(variable = v, extreme = "Lowest")
  ) %>%
    select(variable, extreme, everything())
  
  p <- ggplot(sf_v) +
    geom_sf(aes(fill = hotspot_class), color = "white", linewidth = 0.2) +
    scale_fill_manual(
      values = c(
        "Coldspot (99%)"  = "#2166ac",
        "Coldspot (95%)"  = "#67a9cf",
        "Not significant" = "grey90",
        "Hotspot (95%)"   = "#ef8a62",
        "Hotspot (99%)"   = "#b2182b"
      ),
      drop = FALSE,
      name = "Gi* result"
    ) +
    labs(
      title = paste("Hotspots and Coldspots of", v, "in California"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.title = element_text(face = "bold")
    )
  
  print(p)
  
  ggsave(
    filename = file.path(output_dir, paste0(v, "_hotspots.png")),
    plot = p,
    width = 12,
    height = 9,
    dpi = 300,
    bg = "white"
  )
}

## ---Export CSVs---

hotspot_summary <- bind_rows(hotspot_tables)
topbottom_summary <- bind_rows(topbottom_tables)

write_csv(county_mean, file.path(output_dir, "county_mean_values_2018_2023.csv"))
write_csv(hotspot_summary, file.path(output_dir, "county_hotspot_summary_2018_2023.csv"))
write_csv(topbottom_summary, file.path(output_dir, "county_top_bottom_summary_2018_2023.csv"))

## ---Print Top/Bot Counties---

for (v in vars) {
  cat("\n====================================================\n")
  cat(v, "\n")
  
  cat("\nHighest counties:\n")
  print(
    topbottom_summary %>%
      filter(variable == v, extreme == "Highest") %>%
      select(County, GEOID, value, hotspot_class, gi_z, p_value)
  )
  
  cat("\nLowest counties:\n")
  print(
    topbottom_summary %>%
      filter(variable == v, extreme == "Lowest") %>%
      select(County, GEOID, value, hotspot_class, gi_z, p_value)
  )
}