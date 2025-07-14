library(data.table)
library(dplyr)
library(terra)
library(sf)
library(tmap)
library(ggplot2)


gep_gwr_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"


# /*===========================================*/
#'= World Data Grids (from SWAT+) =
# /*===========================================*/

# /*===== Country Boundary Data =====*/
world_bd <- 
  st_read(file.path(gep_gwr_dt_dir, "Raw_all/World_Data_Grids/country.shp"))
st_crs(world_bd) = 4326

tm_shape(world_bd) +
  tm_borders()

# /*===== DEM SRTM =====*/
dem_srtm <- st_read(file.path(gep_gwr_dt_dir, "Raw/World_Data_Grids/dem_srtm_grid.shp"))

tm_shape(dem_srtm) +
  tm_borders()


# /*===========================================*/
#'=  WHYMAP =
# /*===========================================*/
# /*===== MAJOR GROUNDWATER BASINS =====*/
aquifer_bd <- 
  st_read(file.path(gep_gwr_dt_dir, "Raw/WHYMAP_GWR/shp/whymap_GW_aquifers_v1_poly.shp")) %>%
  setnames(., tolower(names(.))) %>%
  # Remove Antarctica
  filter(continent != "99") %>%
  filter(hygeo2 %in% 12:15)
  
ggplot() +
  geom_sf(data = world_bd) +
  geom_sf(
    data = aquifer_bd,
    aes(fill = factor(hygeo2))
  )
         
# aquifer_bd_v2 <- 
#   st_read(file.path(gep_gwr_dt_dir, "Raw/WHYMAP_GWR/whymap_gw_aquifers_v2_poly/whymap_gw_aquifers_v2_poly.shp")) %>%
#   setnames(., tolower(names(.))) %>%
#   # Remove Antarctica
#   filter(continent != "99")
  
# sort(unique(aquifer_bd_v2$hygeo2))

# ggplot() +
#   geom_sf(data = world_bd) +
#   geom_sf(
#     # data = aquifer_bd,
#     data = filter(aquifer_bd_v2, hygeo2 %in% c(11:15)),
#     aes(fill = factor(hygeo2))
#   )

#/*--------------------------------*/
#' ## Jasechko_et_al_2024
#/*--------------------------------*/
aquifer_jsechko <- 
  st_read(file.path(gep_gwr_dt_dir, "Raw_all/AquiferBoundary/Jasechko_et_al_2024/sourcedatafig2.shp"))

ggplot() +
  geom_sf(data = world_bd) +
  geom_sf(
    data = aquifer_jsechko,
    fill = "blue"
  )


#/*--------------------------------*/
#' ## Global Geo-processed Data of Aquifer Properties by 0.5° Grid, Country and Water Basins
#/*--------------------------------*/
df_in <- fread(file.path(gep_gwr_dt_dir, "Raw_all/9sany-rht26/inputs/aquifer_properties.csv")) 

aq_bd <- 
  st_read(file.path(gep_gwr_dt_dir, "Raw_all/9sany-rht26/shapefiles/All_merged.shp"))

ggplot() +
  geom_sf(data = st_transform(world_bd, st_crs(aq_bd))) +
  geom_sf(
    data = filter(aq_bd, WHYClass == 10),
    fill = "blue"
  )


inputs_sf <- 
  st_read(file.path(gep_gwr_dt_dir, "Raw_all/9sany-rht26/shapefiles/inputs.shp")) %>%
  st_make_valid()
