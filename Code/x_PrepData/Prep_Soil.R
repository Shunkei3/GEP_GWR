library(data.table)
library(terra)
library(sf)
library(dplyr)
library(ggplot2)
library(tmap)

base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

#/*--------------------------------*/
#' ## Example Soil Group Data
#/*--------------------------------*/
example_soil_group <- rast('/Users/shunkeikakimoto/Dropbox/ResearchProject/GEP/GWR/Data/z_InVest_SampleData/Seasonal_Water_Yield/soil_group_gura.tif')

tm_shape(example_soil_group) +
  tm_raster()

unique(values(example_soil_group))


#/*--------------------------------*/
#' ## TODO:
#' + Create four soil groups (A, B, C, D)
#' + HYSOGs250m contains values: 1,2,3,4,11,12,13,14
#' 
#' To create four soil groups, we need to reclassify the values (1,2,3,4,11,12,13,14) into four groups (A, B, C, D, or 1, 2, 3, 4). Group A is the least runoff potential, while group D is the most runoff potential.
#' + A : 1, Soils having high infiltration rates even when thoroughly wetted, consisting chiefly of deep, well to excessively drained sand and/or gravel. These soils have a high rate of water transmission and would result in a low runoff potential.
#' + B : 2, Soils having moderate infiltration rates when thoroughly wetted, consisting chiefly of moderately deep or deep, moderately well or well drained soils with moderately fine to moderately coarse textures. These soils have a moderate rate of water transmission.
#' + C : 3, Soils having slow infiltration rates when thoroughly wetted, consisting chiefly of (1) soils with a layer that impedes the downward movement of water, or (2) soils with moderately fine or fine textures and slow infiltration rate. These soils have a slow rate of water transmission.
#' + D : 4, Soils having very slow infiltration rates when thoroughly wetted, consisting chiefly of (1) clayey soils with high swelling capacity or potential, (2) soils with a high permanent water table, (3) soils with claypan or clay layer at or near the surface, and (4) shallow soils over nearly impervious materials. These soils have a very slow rate of water transmission.
#/*--------------------------------*/


# /*===========================================*/
#'=  HYSOGs250m =
# /*===========================================*/
hygeo <- rast(file.path(base_dt_dir, "Raw/Soil/Global_Hydrologic_Soil_Group_1566/data/HYSOGs250m.tif"))

tm_shape(hygeo) +
  tm_raster()

# unique(hygeo) # 1  2  3  4 11 12 13 14

reclass_matrix <- 
  matrix(
    c(
      1, 1,
      2, 2,
      3, 3,
      4, 4,
      11, 1,
      12, 2,
      13, 3,
      14, 4
    ), ncol = 2, byrow = TRUE
  )

# Apply reclassification
hygeo_reclassed <- classify(hygeo, rcl = reclass_matrix)

unique(hygeo_reclassed) # 1 2 3 4

writeRaster(
  hygeo_reclassed, 
  filename = file.path(base_dt_dir, "Raw/Soil/hygeo_reclassed.tif"), 
  overwrite = TRUE
)
