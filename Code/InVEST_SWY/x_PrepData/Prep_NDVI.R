# /*===========================================*/
#'= NDVI Data is downloaded from  =
# /*===========================================*/



library(data.table)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)

base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

ls_file_names <- list.files(
  file.path(base_dt_dir, "Raw/NDVI/MOD13C2_061-20250708_040832"), 
  full.names = FALSE
)

ls_ndvi_raw <- 
  lapply(
    ls_file_names, 
    function(x) {
      # x = ls_file_names[1]

      file_path <- file.path(base_dt_dir, "Raw/NDVI/MOD13C2_061-20250708_040832", x)

      layer_raw <- rast(file_path, raw = TRUE) 

      layer <- layer_raw[['\"CMG 0.05 Deg Monthly NDVI\"']]
      
      layer <- layer * 0.0001

      plot(layer)
      
      # Layer names
      day <- as.integer(sub(".*A\\d{4}(\\d{3})\\..*", "\\1", x))
      month <- month(as.Date(day, origin = paste0(substr(x, 10, 13), "-01-01")))
      names(layer) <- paste0("ndvi_m", month)
      
      return(layer)
    }
  ) %>%
  rast()

# plot(ls_ndvi_raw[[5]])

writeRaster(
  ls_ndvi_raw, 
  filename = file.path(base_dt_dir, "Raw/NDVI/ndvi_m_2020.tif"), 
  overwrite = TRUE
)

