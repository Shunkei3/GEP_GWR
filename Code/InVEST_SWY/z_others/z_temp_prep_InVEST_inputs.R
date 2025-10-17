# /*===========================================*/
#'= Objective  =
#' + Prepare input data for each InVEST model
# /*===========================================*/

# --- for data wrangling --- #
library(here)
library(data.table)
library(tidyverse)
# --- for spatial operation --- #
library(sf)
library(terra)
# --- visualization --- #
library(tmap)
httpgd::hgd()
httpgd::hgd_browse()

project_path <- here("seals/projects/project_slv/intermediate")

# /*===== Load boundary of SLV =====*/
bd_slv <- 
	st_read(
		file.path(project_path, "project_aoi/aoi_SLV.gpkg"),
		quiet = TRUE
	)


# /*===== Define functions =====*/
crop_mask_raster_for_aoi <- function(rast_object, ...){	
	crop(
		rast_object, 
		vect(st_transform(bd_slv, crs = crs(rast_object))),
		mask = TRUE
	) %>%
	project(., "epsg:32616")
	# project(., "ESRI:54030")
}

map_raster_aoi <- function(rast_object, ...){
	tm_shape(rast_object)+
	tm_raster(
		style = "cont", 
		palette = "seq"
	)+
	tm_shape(bd_slv)+
	tm_borders(col="blue")
}



#/*--------------------------------*/
#' ## Precipitation
#/*--------------------------------*/
ppt <- rast(here("base_data/mesh/worldclim/baseline/5min/baseline_bio12_Annual_Precipitation.tif"))

ppt_slv <- crop_mask_raster_for_aoi(ppt)

map_raster_aoi(ppt_slv)

writeRaster(
	ppt_slv,
	here("base_data/mesh/worldclim/ppt_slv.tif"),
	overwrite=TRUE
	)


#/*--------------------------------*/
#' ## Evapotranspiration
#/*--------------------------------*/
pet <- rast(here("base_data/mesh/cgiar_csi/pet.tif"))

pet_slv <- crop_mask_raster_for_aoi(pet)

map_raster_aoi(pet_slv)

writeRaster(
	pet_slv,
	here("base_data/mesh/cgiar_csi/pet_slv.tif"),
	overwrite=TRUE
	)

#/*--------------------------------*/
#' ## Watershed
#/*--------------------------------*/
# test <- 
# 	st_read(here("base_data/invest_sample_data/Annual_Water_Yield/watershed_gura.shp"))

hydrobasin_lev06 <- 
	st_read(here("base_data/mesh/hydrosheds/hydrobasins/hybas_na_lev01-06_v1c/hybas_na_lev06_v1c.shp"))

# tm_shape(hydrobasin_lev)+
# 	tm_borders() + 
# tm_shape(bd_slv)+
# 	tm_polygons(col="blue", alpha=0.5)

sf_use_s2(FALSE)
hydrobasin_slv_lev06 <- 
	st_crop(hydrobasin_lev06, bd_slv, crop = TRUE) %>%
	st_transform(., crs = 32616) %>%
	# st_transform(., crs = st_crs("ESRI:54030")) %>%
	# --- rename some columns --- #
	rename(ws_id = "HYBAS_ID")
	
# tm_shape(hydrobasin_slv_lev06)+
# 	tm_borders() + 
# tm_shape(bd_slv)+
# 	tm_polygons(col="blue", alpha=0.5)
st_write(
	hydrobasin_slv_lev06, 
	here("base_data/mesh/hydrosheds/hydrobasin_slv_lev06.shp"),
	delete_layer = TRUE
)

#/*--------------------------------*/
#' ## Root Restricting Layer Depth
#/*--------------------------------*/
depth_to_root <- 
	rast(here("base_data/mesh/isric/depth_to_root_restricting_layer.tif"))

depth_to_root_slv <- crop_mask_raster_for_aoi(depth_to_root)

map_raster_aoi(depth_to_root_slv)

writeRaster(
	depth_to_root_slv,
	here("base_data/mesh/isric/depth_to_root_restricting_layer_slv.tif"),
	overwrite=TRUE
	)


#/*--------------------------------*/
#' ## Plant Available Water Content
#/*--------------------------------*/
pawc <- rast(here("base_data/mesh/soil/pawc_30s.tif"))

pawc_slv <- crop_mask_raster_for_aoi(pawc)

map_raster_aoi(pawc_slv)


writeRaster(
	pawc_slv,
	here("base_data/mesh/soil/pawc_slv.tif"),
	overwrite=TRUE
	)

#/*--------------------------------*/
#' ## Biophysical Table
#/*--------------------------------*/
# /*===== Load tables =====*/
esa_seals7_tb <-
	fread(here("base_data/seals/default_inputs/esa_seals7_correspondence.csv")) %>%
	setnames("src_label", "LULC_desc") %>%
	.[!duplicated(LULC_desc)]
# any(duplicated(esa_seals7_tb$src_label))

bio_tb <- 
	fread(here("base_data/mesh/esa_and_modis_biophysical_table.csv")) %>%
	.[!is.na(esacci_lucode)] %>%
	.[,names(.)[!duplicated(names(.))], with = FALSE]

# bio_tb2 <- 
# 	readxl::read_excel(here("base_data/mesh/Koh_et_al_2016_pnas.xlsx")) %>%
# 	data.table()

bio_pollination <- 
	fread(here("base_data/mesh/landcover_biophysical_table.csv")) %>%
	setnames("src_label", "LULC_desc")


bio_seals7_dt <- 
	left_join(esa_seals7_tb,bio_tb, by = "LULC_desc") %>%
	# esa_seals7_tb[bio_tb, on=.(LULC_desc)] %>%
	# bio_tb[esa_seals7_tb, on=.(LULC_desc)] %>%
	.[, `:=`(
		dst_index = NULL,
		native_veg = NULL,
		LULC_desc = NULL,
		modis_lucode = NULL
	)] %>%
	.[, lapply(.SD, mean, na.rm=TRUE), by = .(dst_label, dst_id)] %>%
	.[order(dst_id)] %>%
	setnames("dst_id", "lucode") %>%
	left_join(., bio_pollination, by = "lucode")

#' lulc_veg (integer, required): 
#' Code indicating whether the the LULC class is vegetated 
#' for the purpose of AET. Enter 1 for all vegetated classes 
#' except wetlands, and 0 for all other classes, including 
#' wetlands, urban areas, water bodies, etc.

# --- update LULC_veg column --- #
ls_veg_lucode <- c(2, 3, 4)
bio_seals7_dt[,
	LULC_veg := ifelse(lucode %in% ls_veg_lucode, 1, 0)
]

bio_seals7_dt[is.na(bio_seals7_dt)] <- 0


write.csv(
	bio_seals7_dt, 
	here("base_data/mesh/biophysical_table.csv")
)


#/*--------------------------------*/
#' ## Elevation 
#/*--------------------------------*/
elev <- rast(here("base_data/seals/static_regressors/alt_m.tif"))

slv_elev <- crop_mask_raster_for_aoi(elev)

map_raster_aoi(slv_elev)

writeRaster(
	slv_elev,
	here("base_data/mesh/DEM/slv_elev.tif"),
	overwrite=TRUE
	)


plot(slv_elev)
tm_shape(slv_elev)+tm_raster()

#/*--------------------------------*/
#' ## Rainfall Erosivity
#/*--------------------------------*/
slv_erosivity <- 
	rast(here("base_data/mesh/IEEM-ES-SLV-main/R_factor_rainfall_erosivity/slv_R_factor_rainfall_erosivity.tif"))


#/*--------------------------------*/
#' ## Soil Erodibility
#/*--------------------------------*/
slv_soil_erodibility <- 
	rast(here("base_data/mesh/IEEM-ES-SLV-main/K_factor_soil_erodibility/slv_K_factor_soil_erodibility.tif"))















