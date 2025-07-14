# /*===========================================*/
#'=  This code is unnecessary because I create the base biophysical table by myself  =
#' See Raw/Biophysical_Tables/biophysical_table_base.csv
# /*===========================================*/


library(data.table)
library(dplyr)

library(terra)

base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

# lulc <- rast(file.path(base_dt_dir, "Raw/LULC/lulc_esa_2020.tif"))

# unique_lulc_codes <- unique(lulc)
# c(unique_lulc_codes)

#/*--------------------------------*/
#' ## Biophysical Table for SWY
#' A table mapping each LULC code to biophysical properties of the corresponding LULC class.

#' "All LULC classes in the LULC raster MUST have corresponding values in this table."

#'  Each row is a land use/land cover class and columns must be named and defined as follows:

#' lucode: 
#' + LULC codes from the LULC raster. Each code must be a unique integer.
#' 
#' cn_[SOIL_GROUP]:
#' + Soil group code: A, B, C, D,
#' + Curve number values for each combination of soil group and LULC class.
#' 
#' kc_[MONTH]:
#' + Crop/vegetation coefficient (Kc) values for this LULC class in each month. 
#/*--------------------------------*/


#/*--------------------------------*/
#' ## Reference:
#' + Biophysical_Colombia_plus_Carbon_and_SWY_condensed_10Aug2016_SW.xlsx
#' + Kc is "a coefficient expressing the difference in evapotranspiration between the cropped and reference grass surface" (FAO guide).
#/*--------------------------------*/

#/*--------------------------------*/
#' ## What is Curve Number?
#' + The CN method helps estimate the volume of runoff produced by a rainfall event  
#/*--------------------------------*/

#/*--------------------------------*/
#' ## What is crop/vegetation coefficient (Kc)?
#' + The Kc coefficient adjusts reference evapotranspiration (ET₀) to estimate actual crop evapotranspiration (ETc):
#' + Kc values vary by:
#' - Land use/land cover (LULC) class
#' - Crop growth stage
#' - Climate and management conditions
#/*--------------------------------*/


# /*===========================================*/
#'= Build Biophysical Tables = (One time)
# /*===========================================*/
# === This is a base table === #
# lulc_codes <- 
#   fread(file.path(base_dt_dir, "Raw/Biophysical_Tables/esa_and_modis_biophysical_table.csv")) %>%
#   .[!is.na(esacci_lucode)] %>%
#   .[,.(
#     desc = LULC_desc, 
#     lucode = as.integer(esacci_lucode) 
#   )]

# # save as a csv file
# write.csv(
#   lulc_codes, 
#   file.path(base_dt_dir, "Raw/Biophysical_Tables/biophysical_table_base.csv"), 
#   row.names = FALSE
# )






