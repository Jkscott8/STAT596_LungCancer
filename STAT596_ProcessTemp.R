### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 16, 2026
### This code processes the Tmin and Tmax data from CHIRTS-ERA5 for California
### counties and merges them with censusCDC_merged.csv

rm(list = ls())

## ---Packages---

library(terra)
library(sf)
library(dplyr)
library(tidyr)
library(ggplot2)
library(readr)
library(exactextractr)
library(purrr)
library(gganimate)
library(gifski)

## ---Settings & Paths---

terraOptions(tempdir = tempdir())

setwd("E:/STAT596_SP26")

tmin_path <- "E:/STAT596_SP26/tmin/cropped_ca"
tmax_path <- "E:/STAT596_SP26/tmax/cropped_ca"

tmin_file <- function(yr, mo) file.path(tmin_path, sprintf("CHIRTS-ERA5.monthly_Tmin.%04d.%02d.tif", yr, mo))
tmax_file <- function(yr, mo) file.path(tmax_path, sprintf("CHIRTS-ERA5.monthly_Tmax.%04d.%02d.tif", yr, mo))

census_file <- "E:/STAT596_SP26/censusCDC_merged.csv"

## ---Read and clean the census data---

census_df <- read_csv(census_file, show_col_types = FALSE)

census_df_clean <- census_df %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),   # 6001 -> "06001"
    Year = as.integer(Year),
    Month_num = match(trimws(Month), month.name)  # January -> 1, etc.
  )

## ---Build county geometry object---

county_geom <- census_df %>%
  select(GEOID, geometry) %>%
  distinct()

# Convert WKT text to sf geometry
counties_sf <- st_as_sf(county_geom, wkt = "geometry", crs = 3857)
counties_sf <- st_make_valid(counties_sf)

# Keep GEOID as 5-character county FIPS
counties_sf <- counties_sf %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))

## ---Function to clean temperature data---

# -9999.00000 values indicate water pixels that are NOT needed

clean_temp_raster <- function(r, fill_value = NULL) {
  
  if (!is.null(fill_value)) {
    NAflag(r) <- fill_value
  }
  
  r
}

## ---Function to extract county means per month---

extract_county_Tmean <- function(yr, mo, counties_sf) {
  
  message("Processing: ", yr, "-", sprintf("%02d", mo))
  
  # Read and clean rasters
  r_tmin <- rast(tmin_file(yr, mo))
  r_tmax <- rast(tmax_file(yr, mo))
  
  r_tmin <- clean_temp_raster(r_tmin, fill_value = -9999)
  r_tmax <- clean_temp_raster(r_tmax, fill_value = -9999)
  
  # Reproject county polygons to match raster CRS
  counties_use <- st_transform(counties_sf, crs = crs(r_tmin))
  
  # Extract county means
  tmin_mean <- exact_extract(r_tmin, counties_use, "mean")
  tmax_mean <- exact_extract(r_tmax, counties_use, "mean")
  
  tibble(
    GEOID = counties_use$GEOID,
    Year  = yr,
    Month = mo,
    Tmin  = tmin_mean,
    Tmax  = tmax_mean
  )
}

## ---Loop over all months 2018-2023---

years  <- 2018:2023
months <- 1:12

county_climate_monthly <- map_dfr(years, function(yr) {
  map_dfr(months, function(mo) {
    extract_county_Tmean(yr, mo, counties_sf)
  })
})

## ---Check result---

head(county_climate_monthly)

## ---Prepare county climate table---

county_climate_monthly_clean <- county_climate_monthly %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Month = as.integer(Month)
  )

## ---Merge back into census/CDC CSV---

census_merged_out <- census_df_clean %>%
  left_join(
    county_climate_monthly_clean,
    by = c("GEOID", "Year", "Month_num" = "Month")
  )

## ---Write output---

write_csv(census_merged_out, "censusCDC_CHIRTS_merged.csv")

## ---Plotting---

# 1. Make a county-month sf object for all months, 2018-2023
tmin_anim_sf <- county_climate_monthly %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    frame_date = as.Date(sprintf("%04d-%02d-01", Year, Month)),
    frame_lab  = format(frame_date, "%B %Y"),
    Tmin_bin = cut(
      Tmin,
      breaks = c(-Inf, 0, 5, 10, 15, 20, 25, 30, Inf),
      labels = c("< 0", "0 to 5", "5 to 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "> 30"),
      right = FALSE
    )
  ) %>%
  left_join(
    counties_sf %>%
      mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
      st_transform(4326),
    by = "GEOID"
  ) %>%
  st_as_sf()

# Keep frames in chronological order
tmin_anim_sf$frame_lab <- factor(
  tmin_anim_sf$frame_lab,
  levels = format(sort(unique(tmin_anim_sf$frame_date)), "%B %Y")
)

# 2. Build the animated plot
p_anim <- ggplot(tmin_anim_sf) +
  geom_sf(aes(fill = Tmin_bin), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "< 0"      = "#2c7bb6",
      "0 to 5"   = "#74add1",
      "5 to 10"  = "#abd9e9",
      "10 to 15" = "#fee090",
      "15 to 20" = "#fdae61",
      "20 to 25" = "#f46d43",
      "25 to 30" = "#d73027",
      "> 30"     = "#a50026"
    ),
    name = "Tmin (°C)",
    na.value = "grey90",
    drop = FALSE
  ) +
  labs(
    title = "County Mean Tmin in California",
    subtitle = "{current_frame}",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "cm")
  ) +
  transition_manual(frame_lab)

# 3. Render GIF
anim <- animate(
  p_anim,
  nframes = length(levels(tmin_anim_sf$frame_lab)),
  fps = 2,
  width = 900,
  height = 700,
  units = "px",
  end_pause = 8,
  renderer = gifski_renderer()
)

anim

# 4. Save GIF
anim_save(
  filename = "E:/STAT596_SP26/california_tmin_2018_2023.gif",
  animation = anim
)

## ---Plotting August Months Only---

# 1. Make a county-month sf object for August only
tmin_anim_sf_aug <- county_climate_monthly %>%
  filter(Month == 8) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    frame_date = as.Date(sprintf("%04d-%02d-01", Year, Month)),
    frame_lab  = format(frame_date, "%B %Y"),
    Tmin_bin = cut(
      Tmin,
      breaks = c(-Inf, 0, 5, 10, 15, 20, 25, 30, Inf),
      labels = c("< 0", "0 to 5", "5 to 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "> 30"),
      right = FALSE
    )
  ) %>%
  left_join(
    counties_sf %>%
      mutate(GEOID = sprintf("%05d", as.integer(GEOID))) %>%
      st_transform(4326),
    by = "GEOID"
  ) %>%
  st_as_sf()

# Keep frames in chronological order
tmin_anim_sf_aug$frame_lab <- factor(
  tmin_anim_sf_aug$frame_lab,
  levels = format(sort(unique(tmin_anim_sf_aug$frame_date)), "%B %Y")
)

# 2. Build the animated plot
p_anim_aug <- ggplot(tmin_anim_sf_aug) +
  geom_sf(aes(fill = Tmin_bin), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "< 0"      = "#2c7bb6",
      "0 to 5"   = "#74add1",
      "5 to 10"  = "#abd9e9",
      "10 to 15" = "#fee090",
      "15 to 20" = "#fdae61",
      "20 to 25" = "#f46d43",
      "25 to 30" = "#d73027",
      "> 30"     = "#a50026"
    ),
    name = "Tmin (°C)",
    na.value = "grey90",
    drop = FALSE
  ) +
  labs(
    title = "County Mean Tmin in California",
    subtitle = "{current_frame}",
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold"),
    plot.subtitle = element_text(hjust = 0.5, size = 14),
    legend.title = element_text(size = 16, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.1, "cm")
  ) +
  transition_manual(frame_lab)

# 3. Render GIF
anim_aug <- animate(
  p_anim_aug,
  nframes = length(levels(tmin_anim_sf_aug$frame_lab)) + 8,
  fps = 1,
  width = 900,
  height = 700,
  units = "px",
  end_pause = 8,
  renderer = gifski_renderer()
)

anim_aug

# 4. Save GIF
anim_save(
  filename = "E:/STAT596_SP26/california_tmin_August_2018_2023.gif",
  animation = anim_aug
)