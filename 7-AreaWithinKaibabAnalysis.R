#' Tracking solutions to a persistent threat: 
#' Spatial movement patterns reflect lead exposure in critically endangered 
#' California Condors

#' Script by: Varalika Jain

#' The final KML files for this Area within Kaibab analysis are available upon request

####----LIBRARIES----####
library(sf)
library(rgdal)
library(dplyr)

#### (1) Prepare working space----
rm(list = ls())
dev.off()

#### BACKGROUND----
#### (1) Read in KML files ----
# List all KML files in the directory
background <- list.files(path = "./monthly KDE/1 - background/", pattern = "\\.kml$", full.names = TRUE)

# Load all KML files into a list of sf objects
condor_background <- lapply(background, function(file) {
  st_read(file) # Replace "layer_name" with the actual layer name if needed
})

# Combine all the ranges into a single sf object
combined_background <- do.call(rbind, condor_background)
combined_background$Name <- seq.int(nrow(combined_background))

#### (2) Read in Kaibab polygon file ----
# Load the specific polygon (replace with your actual file path)
polygon <- st_read("polygon x.shp")


#### (3) Match crs ----
# Reproject condor ranges to match the CRS of the polygon
combined_background <- st_transform(combined_background, st_crs(polygon))


#### (4) Intersection ----
# Calculate intersection between condor ranges and polygon
intersection_background <- st_intersection(combined_background, polygon)


#### (5) Area proportion ----
# Calculate area of the original condor ranges
background_areas <- st_area(combined_background)

# Calculate area of the intersection
intersection_background_areas <- st_area(intersection_background)

# Calculate the proportion of each range inside the polygon
proportion_background <- intersection_background_areas / background_areas

# Summarize or explore the proportions
summary(proportion_background)

# If you want the overall proportion across all ranges
overall_proportion_background <- sum(intersection_background_areas) / sum(background_areas)

library(ggplot2)

ggplot() +
  geom_sf(data = combined_background, fill = "lightblue", alpha = 0.05, lwd = 0) +
  geom_sf(data = polygon, color = "red", alpha = 0) +
  theme_minimal() +
  labs(title = "Condor Ranges and Specific Polygon Intersection")

#### (6) Area within polygon ----
combined_background$range_area <- st_area(combined_background)

intersection_background$intersected_area <- st_area(intersection_background)

intersection_background_no_geom <- intersection_background %>% 
  st_drop_geometry() %>% 
  select(Name, intersected_area)

# Join the intersected area back to the original data using bird ID or another unique identifier
# (assuming 'bird_id' is a unique identifier in your data)
combined_background <- combined_background %>%
  left_join(intersection_background_no_geom %>% select(Name, intersected_area), by = "Name")

# Calculate the Area-Weighted Concentration Index (AWCI)
combined_background$awci <- combined_background$intersected_area / combined_background$range_area

# Aggregate the AWCI by lead exposure level
mean(combined_background$awci, na.rm = TRUE)



#### EXPOSURE----
#### (1) Read in KML files ----
# List all KML files in the directory
exposure <- list.files(path = "./monthly KDE/2 - exposure/", pattern = "\\.kml$", full.names = TRUE)

# Load all KML files into a list of sf objects
condor_exposure <- lapply(exposure, function(file) {
  st_read(file) # Replace "layer_name" with the actual layer name if needed
})

# Combine all the ranges into a single sf object
combined_exposure <- do.call(rbind, condor_exposure)
combined_exposure$Name <- seq.int(nrow(combined_exposure))

#### (2) Read in Kaibab polygon file ----
# Load the specific polygon (replace with your actual file path)
polygon <- st_read("./Condor overlap map/polygon x.shp")


#### (3) Match crs ----
# Reproject condor ranges to match the CRS of the polygon
combined_exposure <- st_transform(combined_exposure, st_crs(polygon))


#### (4) Intersection ----
# Calculate intersection between condor ranges and polygon
intersection_exposure <- st_intersection(combined_exposure, polygon)


#### (5) Area proportion ----
# Calculate area of the original condor ranges
exposure_areas <- st_area(combined_exposure)

# Calculate area of the intersection
intersection_exposure_areas <- st_area(intersection_exposure)

# Calculate the proportion of each range inside the polygon
proportion_exposure <- intersection_exposure_areas / exposure_areas

# Summarize or explore the proportions
summary(proportion_exposure)

# If you want the overall proportion across all ranges
overall_proportion_exposure <- sum(intersection_exposure_areas) / sum(exposure_areas)

library(ggplot2)

ggplot() +
  geom_sf(data = combined_exposure, fill = "lightblue", alpha = 0.05, lwd = 0) +
  geom_sf(data = polygon, color = "red", alpha = 0) +
  theme_minimal() +
  labs(title = "Condor Ranges and Specific Polygon Intersection")

#### (6) Area within polygon----
combined_exposure$range_area <- st_area(combined_exposure)

intersection_exposure$intersected_area <- st_area(intersection_exposure)

intersection_exposure_no_geom <- intersection_exposure %>% 
  st_drop_geometry() %>% 
  select(Name, intersected_area)

# Join the intersected area back to the original data using bird ID or another unique identifier
# (assuming 'bird_id' is a unique identifier in your data)
combined_exposure <- combined_exposure %>%
  left_join(intersection_exposure_no_geom %>% select(Name, intersected_area), by = "Name")

# Calculate the Area-Weighted Concentration Index (AWCI)
combined_exposure$awci <- combined_exposure$intersected_area / combined_exposure$range_area

# Aggregate the AWCI by lead exposure level
mean(combined_exposure$awci, na.rm = TRUE)


#### HIGH EXPOSURE----
#### (1) Read in KML files ----
# List all KML files in the directory
h.exposure <- list.files(path = "./monthly KDE/3 - high exposure/", pattern = "\\.kml$", full.names = TRUE)

# Load all KML files into a list of sf objects
condor_h.exposure <- lapply(h.exposure, function(file) {
  st_read(file) # Replace "layer_name" with the actual layer name if needed
})

# Combine all the ranges into a single sf object
combined_h.exposure <- do.call(rbind, condor_h.exposure)
combined_h.exposure$Name <- seq.int(nrow(combined_h.exposure))

#### (2) Read in Kaibab polygon file ----
# Load the specific polygon (replace with your actual file path)
polygon <- st_read("./Condor overlap map/polygon x.shp")


#### (3) Match crs ----
# Reproject condor ranges to match the CRS of the polygon
combined_h.exposure <- st_transform(combined_h.exposure, st_crs(polygon))
combined_h.exposure <- st_make_valid(combined_h.exposure)


#### (4) Intersection ----
# Calculate intersection between condor ranges and polygon
intersection_h.exposure <- st_intersection(combined_h.exposure, polygon)


#### (5) Area proportion ----
# Calculate area of the original condor ranges
h.exposure_areas <- st_area(combined_h.exposure)

# Calculate area of the intersection
intersection_h.exposure_areas <- st_area(intersection_h.exposure)

# Calculate the proportion of each range inside the polygon
proportion_h.exposure <- intersection_h.exposure_areas / h.exposure_areas

# Summarize or explore the proportions
summary(proportion_h.exposure)

# If you want the overall proportion across all ranges
overall_proportion_h.exposure <- sum(intersection_h.exposure_areas) / sum(h.exposure_areas)

library(ggplot2)

ggplot() +
  geom_sf(data = combined_h.exposure, fill = "lightblue", alpha = 0.05, lwd = 0) +
  geom_sf(data = polygon, color = "red", alpha = 0) +
  theme_minimal() +
  labs(title = "Condor Ranges and Specific Polygon Intersection")

#### (6) Area within polygon----
combined_h.exposure$range_area <- st_area(combined_h.exposure)

intersection_h.exposure$intersected_area <- st_area(intersection_h.exposure)

intersection_h.exposure_no_geom <- intersection_h.exposure %>% 
  st_drop_geometry() %>% 
  select(Name, intersected_area)

# Join the intersected area back to the original data using bird ID or another unique identifier
# (assuming 'bird_id' is a unique identifier in your data)
combined_h.exposure <- combined_h.exposure %>%
  left_join(intersection_h.exposure_no_geom %>% select(Name, intersected_area), by = "Name")

# Calculate the Area-Weighted Concentration Index (AWCI)
combined_h.exposure$awci <- combined_h.exposure$intersected_area / combined_h.exposure$range_area

# Aggregate the AWCI by lead h.exposure level
mean(combined_h.exposure$awci, na.rm = TRUE)



#### CLINICALLY AFFECTED----
#### (1) Read in KML files ----
# List all KML files in the directory
cl.aff <- list.files(path = "./monthly KDE/4 - clinically affected/", pattern = "\\.kml$", full.names = TRUE)

# Load all KML files into a list of sf objects
condor_cl.aff <- lapply(cl.aff, function(file) {
  st_read(file) # Replace "layer_name" with the actual layer name if needed
})

# Combine all the ranges into a single sf object
combined_cl.aff <- do.call(rbind, condor_cl.aff)
combined_cl.aff$Name <- seq.int(nrow(combined_cl.aff))

#### (2) Read in Kaibab polygon file ----
# Load the specific polygon (replace with your actual file path)
polygon <- st_read("./Condor overlap map/polygon x.shp")


#### (3) Match crs ----
# Reproject condor ranges to match the CRS of the polygon
combined_cl.aff <- st_transform(combined_cl.aff, st_crs(polygon))


#### (4) Intersection ----
# Calculate intersection between condor ranges and polygon
intersection_cl.aff <- st_intersection(combined_cl.aff, polygon)


#### (5) Area proportion ----
# Calculate area of the original condor ranges
cl.aff_areas <- st_area(combined_cl.aff)

# Calculate area of the intersection
intersection_cl.aff_areas <- st_area(intersection_cl.aff)

# Calculate the proportion of each range inside the polygon
proportion_cl.aff <- intersection_cl.aff_areas / cl.aff_areas

# Summarize or explore the proportions
summary(proportion_cl.aff)

# If you want the overall proportion across all ranges
overall_proportion_cl.aff <- sum(intersection_cl.aff_areas) / sum(cl.aff_areas)

library(ggplot2)

ggplot() +
  geom_sf(data = combined_cl.aff, fill = "lightblue", alpha = 0.05, lwd = 0) +
  geom_sf(data = polygon, color = "red", alpha = 0) +
  theme_minimal() +
  labs(title = "Condor Ranges and Specific Polygon Intersection")

#### (6) Area within polygon----
combined_cl.aff$range_area <- st_area(combined_cl.aff)

intersection_cl.aff$intersected_area <- st_area(intersection_cl.aff)

intersection_cl.aff_no_geom <- intersection_cl.aff %>% 
  st_drop_geometry() %>% 
  select(Name, intersected_area)

# Join the intersected area back to the original data using bird ID or another unique identifier
# (assuming 'bird_id' is a unique identifier in your data)
combined_cl.aff <- combined_cl.aff %>%
  left_join(intersection_cl.aff_no_geom %>% select(Name, intersected_area), by = "Name")

# Calculate the Area-Weighted Concentration Index (AWCI)
combined_cl.aff$awci <- combined_cl.aff$intersected_area / combined_cl.aff$range_area

# Aggregate the AWCI by lead cl.aff level
mean(combined_cl.aff$awci, na.rm = TRUE)



#### AREA WITHIN POLYGON CONCENTRATION FINAL ----
# Assuming 'lead_exposure_level' is a column in your data

acwi <- c(mean(combined_background$awci, na.rm = TRUE),
          mean(combined_exposure$awci, na.rm = TRUE),
          mean(combined_h.exposure$awci, na.rm = TRUE),
          mean(combined_cl.aff$awci, na.rm = TRUE))


c(range(combined_background$awci, na.rm = TRUE),
          range(combined_exposure$awci, na.rm = TRUE),
          range(combined_h.exposure$awci, na.rm = TRUE),
          range(combined_cl.aff$awci, na.rm = TRUE))

library(units)
ggplot(combined_background, aes(x =as.numeric(awci))) +  
  xlim(c(0,1))+ylim(c(0,15))+  
  geom_histogram(color = "#000000", fill = "grey", bins = 35)+
  theme_minimal()


ggplot(combined_exposure, aes(x =as.numeric(awci))) + 
  xlim(c(0,1))+ylim(c(0,15))+  
  geom_histogram(color = "#000000", fill = "grey", bins = 35)+
  theme_minimal()

ggplot(combined_h.exposure, aes(x =as.numeric(awci))) + 
  xlim(c(0,1))+ylim(c(0,15))+  
  geom_histogram(color = "#000000", fill = "grey", bins = 35)+
  theme_minimal()

ggplot(combined_cl.aff, aes(x =as.numeric(awci))) + 
  xlim(c(0,1))+ylim(c(0,15))+  
  geom_histogram(color = "#000000", fill = "grey", bins = 35)+
  theme_minimal()
