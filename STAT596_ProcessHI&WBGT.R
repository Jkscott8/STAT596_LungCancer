### STAT 596: Spatiotemporal Data Analysis
### Spring 2026
### Written by: Paolo Urbano
### Last Updated: April 16, 2026
### This code processes the HI (heat index) and WBGT data from CHIRTS-ERA5 for 
### California counties and merges them with censusCDC_CHIRTS_merged.csv

rm(list = ls())

## ---Packages---

library(terra)
library(sf)
library(dplyr)
library(readr)
library(exactextractr)
library(purrr)

## ---Settings & Paths---

terraOptions(tempdir = tempdir())

setwd("E:/STAT596_SP26")

# Existing already-merged file from your Tmin/Tmax workflow
census_file <- "E:/STAT596_SP26/censusCDC_CHIRTS_merged.csv"

# Daily raster folders
hi_path   <- "E:/STAT596_SP26/HI/cropped_ca"
wbgt_path <- "E:/STAT596_SP26/WBGT/cropped_ca"

# Daily file naming functions
hi_file   <- function(yr, mo, dy) file.path(hi_path,   sprintf("HI.%04d.%02d.%02d.tif", yr, mo, dy))
wbgt_file <- function(yr, mo, dy) file.path(wbgt_path, sprintf("WBGT.%04d.%02d.%02d.tif", yr, mo, dy))

## ---Read existing merged census/CDC/CHIRTS file---

census_df <- read_csv(census_file, show_col_types = FALSE)

census_df_clean <- census_df %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Month_num = match(trimws(Month), month.name)
  )

## ---Build county geometry object---

county_geom <- census_df %>%
  select(GEOID, geometry) %>%
  distinct()

# Convert WKT text to sf geometry
counties_sf <- st_as_sf(county_geom, wkt = "geometry", crs = 3857)
counties_sf <- st_make_valid(counties_sf)

counties_sf <- counties_sf %>%
  mutate(GEOID = sprintf("%05d", as.integer(GEOID)))

## ---Helper functions---

# Set fill value to NA
clean_raster <- function(r, fill_value = -9999) {
  NAflag(r) <- fill_value
  r
}

# Get all day numbers in a month
get_days_in_month <- function(yr, mo) {
  start_date <- as.Date(sprintf("%04d-%02d-01", yr, mo))
  next_month <- seq(start_date, by = "month", length.out = 2)[2]
  end_date <- next_month - 1
  as.integer(format(seq(start_date, end_date, by = "day"), "%d"))
}

# Build one monthly mean raster from daily rasters
make_monthly_mean_raster <- function(yr, mo, file_fun, fill_value = -9999) {
  
  days <- get_days_in_month(yr, mo)
  
  daily_files <- map_chr(days, ~file_fun(yr, mo, .x))
  daily_files <- daily_files[file.exists(daily_files)]
  
  
  r_daily <- rast(daily_files)
  r_daily <- clean_raster(r_daily, fill_value = fill_value)
  
  # Average across all daily layers
  r_monthly <- app(r_daily, mean, na.rm = TRUE)
  
  return(r_monthly)
}

# Extract county monthly mean from daily raster set
extract_county_monthly_from_daily <- function(yr, mo, counties_sf, file_fun, var_name,
                                              fill_value = -9999) {
  
  message("Processing ", var_name, ": ", yr, "-", sprintf("%02d", mo))
  
  r_monthly <- make_monthly_mean_raster(
    yr = yr,
    mo = mo,
    file_fun = file_fun,
    fill_value = fill_value
  )
  
  counties_use <- st_transform(counties_sf, crs = crs(r_monthly))
  
  var_mean <- exact_extract(r_monthly, counties_use, "mean")
  
  out <- tibble(
    GEOID = counties_use$GEOID,
    Year = yr,
    Month = mo,
    value = var_mean
  )
  
  names(out)[4] <- var_name
  return(out)
}

## ---Years and months to process---

years  <- 2018:2023
months <- 1:12

## ---Extract monthly county HI---

county_hi_monthly <- map_dfr(years, function(yr) {
  map_dfr(months, function(mo) {
    extract_county_monthly_from_daily(
      yr = yr,
      mo = mo,
      counties_sf = counties_sf,
      file_fun = hi_file,
      var_name = "HI",
      fill_value = -9999
    )
  })
})

## ---Extract monthly county WBGT---

county_wbgt_monthly <- map_dfr(years, function(yr) {
  map_dfr(months, function(mo) {
    extract_county_monthly_from_daily(
      yr = yr,
      mo = mo,
      counties_sf = counties_sf,
      file_fun = wbgt_file,
      var_name = "WBGT",
      fill_value = -9999
    )
  })
})

## ---Combine HI and WBGT tables---

county_heat_monthly <- county_hi_monthly %>%
  left_join(county_wbgt_monthly, by = c("GEOID", "Year", "Month")) %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    Year = as.integer(Year),
    Month = as.integer(Month)
  )

## ---Join onto existing census file---

census_final_out <- census_df_clean %>%
  left_join(
    county_heat_monthly,
    by = c("GEOID", "Year", "Month_num" = "Month")
  )

## ---Write output---

write_csv(census_final_out, "censusCDC_CHIRTS_merged.csv")

## ---Plotting HI---

# 1. Make a county-month sf object for all months, 2018-2023
hi_anim_sf <- county_heat_monthly %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    frame_date = as.Date(sprintf("%04d-%02d-01", Year, Month)),
    frame_lab  = format(frame_date, "%B %Y"),
    HI_bin = cut(
      HI,
      breaks = c(-Inf, 10, 15, 20, 25, 30, 35, 40, Inf),
      labels = c("< 10", "10 to 15", "15 to 20", "20 to 25", "25 to 30", "30 to 35", "35 to 40", "> 40"),
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
hi_anim_sf$frame_lab <- factor(
  hi_anim_sf$frame_lab,
  levels = format(sort(unique(hi_anim_sf$frame_date)), "%B %Y")
)

# 2. Build the animated plot
p_hi_anim <- ggplot(hi_anim_sf) +
  geom_sf(aes(fill = HI_bin), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "< 10"      = "#2c7bb6",
      "10 to 15"  = "#74add1",
      "15 to 20"  = "#abd9e9",
      "20 to 25"  = "#fee090",
      "25 to 30"  = "#fdae61",
      "30 to 35"  = "#f46d43",
      "35 to 40"  = "#d73027",
      "> 40"      = "#a50026"
    ),
    name = "HI (°C)",
    na.value = "grey90",
    drop = FALSE
  ) +
  labs(
    title = "County Mean Heat Index in California",
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
anim_hi <- animate(
  p_hi_anim,
  nframes = length(levels(hi_anim_sf$frame_lab)),
  fps = 2,
  width = 900,
  height = 700,
  units = "px",
  end_pause = 8,
  renderer = gifski_renderer()
)

anim_hi

# 4. Save GIF
anim_save(
  filename = "E:/STAT596_SP26/california_HI_2018_2023.gif",
  animation = anim_hi
)

## ---Plotting WBGT---

# 1. Make a county-month sf object for all months, 2018-2023
wbgt_anim_sf <- county_heat_monthly %>%
  mutate(
    GEOID = sprintf("%05d", as.integer(GEOID)),
    frame_date = as.Date(sprintf("%04d-%02d-01", Year, Month)),
    frame_lab  = format(frame_date, "%B %Y"),
    WBGT_bin = cut(
      WBGT,
      breaks = c(-Inf, 10, 15, 20, 25, 28, 30, 32, Inf),
      labels = c("< 10", "10 to 15", "15 to 20", "20 to 25", "25 to 28", "28 to 30", "30 to 32", "> 32"),
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
wbgt_anim_sf$frame_lab <- factor(
  wbgt_anim_sf$frame_lab,
  levels = format(sort(unique(wbgt_anim_sf$frame_date)), "%B %Y")
)

# 2. Build the animated plot
p_wbgt_anim <- ggplot(wbgt_anim_sf) +
  geom_sf(aes(fill = WBGT_bin), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "< 10"      = "#2c7bb6",
      "10 to 15"  = "#74add1",
      "15 to 20"  = "#abd9e9",
      "20 to 25"  = "#fee090",
      "25 to 28"  = "#fdae61",
      "28 to 30"  = "#f46d43",
      "30 to 32"  = "#d73027",
      "> 32"      = "#a50026"
    ),
    name = "WBGT (°C)",
    na.value = "grey90",
    drop = FALSE
  ) +
  labs(
    title = "County Mean WBGT in California",
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
anim_wbgt <- animate(
  p_wbgt_anim,
  nframes = length(levels(wbgt_anim_sf$frame_lab)),
  fps = 2,
  width = 900,
  height = 700,
  units = "px",
  end_pause = 8,
  renderer = gifski_renderer()
)

anim_wbgt

# 4. Save GIF
anim_save(
  filename = "E:/STAT596_SP26/california_WBGT_2018_2023.gif",
  animation = anim_wbgt
)
