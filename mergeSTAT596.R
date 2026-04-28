
# Read in the data
census = read_csv("~/Downloads/ca_census.csv")
lung_cancer = read_csv("~/Downloads/Lung Cancer 2018-2023.csv")


# Reshape the data frame
lung_cancer_long <- lung_cancer |>
  pivot_longer(cols = -c("County", "County Code"), names_to = "Month_Year", values_to = "Lung_Cancer_Rate") |>
  separate(col = Month_Year, into = c("Month", "Year"), sep = " ") |>
  mutate(Year = as.integer(Year)) |> 
  mutate(GEOID = `County Code` ) |> 
  select(-`County Code`)

# Get names for columns to know what to merge
names(lung_cancer_long)
names(census)

# Merge by GEOID & Year 
merged = left_join(lung_cancer_long, census, by = c('GEOID' = 'GEOID','Year' = 'Data_Year'))
write_csv(merged, 'census_cdc_merged.csv')


# Re Add Polygons since there was an issue earlier

library(sf)
library(spnaf)

censusCDC_geom = left_join(merged |> select(-geometry), CA_polygon |> mutate(id = as.integer(id)), by = c('COUNTYFP'='id'))
str(censusCDC_geom)



