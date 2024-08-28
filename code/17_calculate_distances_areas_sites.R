library(here)
library(tidyverse)
library(sf)
library(leaflet)

setwd(here::here("data"))

load("neon_mammal_box_trapping_v01.RData")

d <- readr::read_csv("disease_with_biodiversity_metrics_v01.csv")

coords <- boxtrap[[5]] |> 
  tibble::as_tibble() |> 
  dplyr::select(siteID, plotID,
                y = decimalLatitude, x = decimalLongitude) |> 
  dplyr::distinct() |> 
  dplyr::filter(siteID %in% unique(d$siteID)) |> 
  dplyr::group_by(siteID, plotID) |> 
  dplyr::summarise(x = mean(x), 
                   y = mean(y)) |> 
  sf::st_as_sf(coords = c("x", "y"),
               crs = 4326) 

# Function to calculate pairwise distances within a single site
calculate_pairwise_distances <- function(site_data) {
  # Compute pairwise distances
  distances <- st_distance(site_data, site_data) |>
    as.data.frame() |>
    tibble::rownames_to_column("point1") |>
    tidyr::pivot_longer(-point1, names_to = "point2", values_to = "distance") |>
    dplyr::filter(point1 != point2) |>
    dplyr::mutate(siteID = unique(site_data$siteID)) |>
    dplyr::arrange(siteID, point1, point2)
  
  distances
}

# Group data by siteID and compute distances for each group
all_distances <- coords |>
  dplyr::group_by(siteID) |>
  dplyr::group_split() |>
  purrr::map_dfr(calculate_pairwise_distances)

# View the results
all_distances |> 
  dplyr::mutate(distance = as.numeric(distance) / 1000 ) |> 
  filter(!distance == 0) |> 
  dplyr::select(siteID, distance) |> 
  dplyr::distinct() |> 
  group_by(siteID) |> 
  
  dplyr::summarise(mean_dist = mean(distance), 
                   sd_dist = sd(distance)) 

calculate_mcp_area <- function(site_data) {
  # Compute Minimum Convex Polygon (Convex Hull)
  mcp <- sf::st_convex_hull(sf::st_union(site_data))
  
  # Calculate the area of the MCP
  area <- sf::st_area(mcp)
  
  # Return the siteID and area
  tibble::tibble(siteID = unique(site_data$siteID), mcp_area = as.numeric(area))
}

mcp_areas <- coords |>
  dplyr::group_by(siteID) |>
  dplyr::group_split() |>
  purrr::map_dfr(calculate_mcp_area)

# View the results
#GRSM is huge
mcp_areas |> 
  mutate(area_km2 = mcp_area/1e6) |> 
  ggplot(aes(x = area_km2)) + 
  geom_histogram()

# plot for a sanity check
site_data <- coords |> 
  filter(siteID == "GRSM")

map <- leaflet() |>
  addProviderTiles(providers$OpenStreetMap) |>  # Add base map tiles
  addCircleMarkers(data = site_data,            # Add point data
                   color = "red",               # Color of the markers
                   radius = 5,                  # Radius of the markers
                   fillOpacity = 0.8,           # Fill opacity
                   popup = ~plotID)             # Popup with plotID on click

# Print the map to view it
print(map)

site_distance_area <- all_distances |>
  dplyr::mutate(distance = as.numeric(distance)/1000) |> 
  dplyr::filter(!distance == 0) |> 
  group_by(siteID) |> 
  dplyr::summarise(
    min = min(distance),
    mean = mean(distance),
    med = median(distance),
    max = max(distance),
    sd = sd(distance)) |> 
  dplyr::left_join(mcp_areas) |> 
  dplyr::mutate(area_km2 = mcp_area/1e6) |> 
  dplyr::select(-mcp_area)

setwd(here::here("data"))
write_csv(site_distance_area, "site_distance_area.csv")