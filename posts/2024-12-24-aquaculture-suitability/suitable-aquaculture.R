# Author: Leilanie Rubinstein
# Date: December 3, 2024

# Custom function to source into `hwk4.qmd`
# Given a marine species and its max and min sea surface temperature and depth
# this function creates a map of EEZ regions colored by amount of suitable area


# `min_temp` Minimum temperature threshold in Celsius
# `max_temp` Maximum temperature threshold in Celsius
# `min_depth` Minimum depth threshold in meters (negative values)
# `max_depth` Maximum depth threshold in meters (negative values)
# `species` Scientific name of species
# retuns a tmap plot showing suitable areas by EEZ region

suitable_aquaculture <- function(min_temp, max_temp, min_depth, max_depth, species) {
  # Read and prepare data
  files <- list.files(here::here("data"), pattern = "average*", full.names = TRUE)
  sst <- terra::rast(files)
  names(sst) <- c(2008, 2009, 2010, 2011, 2012)
  
  # Import bathymetry data
  depth <- terra::rast(here::here("data", "depth.tif"))
  
  # Import EEZ and state boundaries
  wc_eez <- st_read(here::here("data", "wc_regions_clean.shp"), quiet = TRUE) %>%
    st_transform(., crs = crs(sst))
  
  states <- st_read(here::here("data", "states", "cb_2023_us_state_20m.shp"), quiet = TRUE) %>%
    filter(NAME %in% c("California", "Oregon", "Washington", "Nevada")) %>%
    st_transform(., crs = crs(sst)) %>%
    st_crop(., st_bbox(wc_eez))
  
  # Process SST data
  sst_mean <- sst %>%
    terra::mean() %>%
    - 273.15  # Convert from Kelvin to Celsius
  
  # Process depth data
  depth <- project(depth, y = crs(sst))
  depth_cropped <- crop(depth, sst_mean)
  depth_resample <- resample(depth_cropped, sst_mean, method = "near")
  
  # Create reclassification matrices
  rcl_depth <- matrix(c(-Inf, -max_depth, 0,
                        -max_depth, min_depth, 1,
                        min_depth, Inf, 0),
                      ncol = 3, byrow = TRUE)
  rcl_sst <- matrix(c(-Inf, min_temp, 0,
                      min_temp, max_temp, 1,
                      max_temp, Inf, 0),
                    ncol = 3, byrow = TRUE)
  
  # Apply matrices to depth and SST rasters
  depth_rcl <- terra::classify(depth_resample, rcl = rcl_depth)
  sst_rcl <- terra::classify(sst_mean, rcl = rcl_sst)
  
  # Find locations that satisfy both SST and depth conditions
  suitablility <- function(depth_rcl, sst_rcl) {
    depth_rcl*sst_rcl
  }
  suitable <- lapp(c(depth_rcl, sst_rcl), fun = suitablility)
  
  # Set unsuitable locations to NA
  suitable[suitable == 0] <- NA
  
  # Find the total suitable area within each EEZ
  eez_area <- terra::cellSize(suitable, mask = T, unit = "km")
  suitable_eez_area <- exactextractr::exact_extract(eez_area, wc_eez, 
                                                    fun = "sum", 
                                                    append_cols=c("rgn", "area_km2"), 
                                                    progress = FALSE)
  
  # Join to original data to obtain geometries for visualization
  suitable_eez_join <- left_join(wc_eez, suitable_eez_area, by = join_by(rgn, area_km2))
  
  # Visualize suitable area
  tm_shape(suitable_eez_join) +
    tm_fill(col = "sum",
            style = "pretty",
            palette = "Blues",
            title = "Suitable Area (km\u00b2)") +
    tm_text(text = "rgn",
            size = 0.7) +
    tm_shape(states) +
    tm_polygons(col = "#e3d3b8",
                border.col = "#402618") +
    tm_shape(wc_eez) +
    tm_borders(col = "#402618") +
    tm_layout(main.title = paste("Suitable Area for", 
                                 species,
                                 "\nFisheries in West Coast EEZ"),
              bg.color = "#e8ebea",
              legend.bg.color = "#e3d3b8",
              legend.position = c(0.61, 0.85)) +
    tm_scale_bar(position = c("left", "bottom")) +
    tm_compass(position = c("right", "bottom"),
               size = 2)
}
