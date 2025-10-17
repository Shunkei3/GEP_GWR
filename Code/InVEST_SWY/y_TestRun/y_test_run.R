# /*===========================================*/
#'=  Run by Watersheds in a continent =
# /*===========================================*/

library(here)
library(data.table)
library(dplyr)
library(terra)
library(exactextractr)
library(sf)
library(tmap)
library(ggplot2)

#/*--------------------------------*/
#' ## Preparation
#/*--------------------------------*/
base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

# === Load Functions === #
source(here("Scr/Code/InVEST_SWY/x_Functions.R"))

# === Load Watersheds === #
case_hybas_tbl <- readRDS(here("Data/Intermediate/case_hybas_tbl.rds"))

test <- tidyr::unnest(case_hybas_tbl, tg_hybas_sf)


tmp_region <- "na"
tmp_case_hybas_sf <- 
  filter(case_hybas_tbl, region_abbr == tmp_region) %>%
  .[["tg_hybas_sf"]] %>%
  .[[1]]

# === Select a watershed === #
i = 1
tmp_hybas_i <- tmp_case_hybas_sf[1,]

tmp_hybas_j <- tmp_case_hybas_sf[2,]


#/*--------------------------------*/
#' ## Test Run
#/*--------------------------------*/
# Prepare Workspace directory
workspace_dir <- file.path(out_dir, "x_output/test_run")
if(!dir.exists(workspace_dir)){
  dir.create(workspace_dir, recursive = TRUE)
}

tmp_region <- tmp_hybas_i$continent
tmp_file_name <- tmp_hybas_i$id

# === AOI === #
# Create .gpkg file in "AOI" folder
fn_save_aoi(aoi = tmp_hybas_i)
test_aoi <- 
  st_read(
    file.path(out_dir, "AOI", tmp_region, paste0(tmp_file_name, ".gpkg"))
  )
   
# === LULC === #
# Create . tif file in "LULC" folder
fn_get_lulc_aoi(aoi = tmp_hybas_i)

test_lulc <- 
  rast(
    file.path(out_dir, "LULC", tmp_region, paste0(tmp_file_name, ".tif"))
  )

# unique(test_lulc, na.rm = FALSE)


# === DEM === #
# Create .tif file in "DEM" folder
fn_get_dem_aoi(aoi = tmp_hybas_i) 

test_dem <- 
  rast(
    file.path(out_dir, "DEM", tmp_region, paste0(tmp_file_name, ".tif"))
  )

# === Flow Accumulation === #
# Get a number for Threshold Flow Accumulation
flow_acc_dt_r <- 
  readRDS(file.path(out_dir, "FlowAcc", paste0("flow_acc_tbl_", tmp_hybas_i$continent, ".rds")))

print(round(flow_acc_dt_r[id == tmp_hybas_i$id, fa_q25]))

# === Monthly Precipitation and Monthly Alpha Table === #
#' Two outputs are created:
#' - A folder containing monthly precipitation (twelve .tif files) in "Precip" folder
#' - Monthly Alpha Table (.csv)
fn_get_month_ppt(aoi = tmp_hybas_i, type = "mean")

test_m_ppt <- 
  rast(
    list.files(
      file.path(out_dir, "Precip", tmp_region, tmp_file_name),
      full.names = TRUE
    )
  ) 

test_m_alpha_tbl <- 
  fread(
    file.path(out_dir, "Alpha_m", tmp_region, paste0(tmp_file_name, ".csv"))
  )

# === Climate Zone and Rain Events Table === #
# Two ouputs are created:
# - A raster file for climate zone (.tif) in "ClimateZone" folder
# - A table for rain events by climate zone (.csv) in "RainEvents" folder
get_rain_events_tbl_cz(aoi = tmp_hybas_i)

test_cz <-
  rast(
    file.path(out_dir, "ClimateZone", tmp_region, paste0(tmp_file_name, ".tif"))
  )

test_rain_events_tbl_cz <-
  fread(
    file.path(out_dir, "RainEvents", tmp_region, paste0(tmp_file_name, ".csv"))
  )

# === Monthly ET0 === #
# Create a folder containing monthly ET0 (twelve .tif files) in "ET0" folder
fn_get_month_et0(aoi = tmp_hybas_i, type = "mean")

test_m_et0 <- 
  rast(
    list.files(
      file.path(out_dir, "ET0", tmp_region, tmp_file_name),
      full.names = TRUE
    )
  )

# === Soil Group === #
# Create a raster file for soil group (.tif) in "SoilGroup" folder
fn_get_soil_group_aoi(aoi = tmp_hybas_i)

test_sg <- 
  rast(
    file.path(out_dir, "SoilGroup", tmp_region, paste0(tmp_file_name, ".tif"))
  )


# === Beta_i parameter === #
Beta_i <- fn_get_beta_twi(aoi = tmp_hybas_i)
Beta_i
# === Gamma parameter === #
Gamma <- fn_get_gamma(aoi = tmp_hybas_i)
Gamma

# /*===== Biophysical Table  =====*/
# Create a biophysical table (.csv) in "Biophysical_Table" folder
fn_gen_biophysical_table(aoi = tmp_hybas_i)

test_biophysical_tbl <- 
  fread(
    file.path(out_dir, "Biophysical_Tables", tmp_region, paste0(tmp_file_name, ".csv")),
    na.strings = c("", "NA", "NaN")
  )





# /*===========================================*/
#'= Sample Data  =
# /*===========================================*/
sample_dt_dir <- '/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/z_InVest_SampleData'

sample_lulc <- 
  rast(file.path(sample_dt_dir, "Seasonal_Water_Yield/land_use_gura.tif"))

sample_biophysical_tbl <- 
  fread(file.path(sample_dt_dir, "Seasonal_Water_Yield/biophysical_table_gura_SWY.csv"))

sample_dem <- 
  rast(file.path(sample_dt_dir, "Seasonal_Water_Yield/dem_gura.tif"))

sample_aoi <- 
  st_read(file.path(sample_dt_dir, "Seasonal_Water_Yield/watershed_gura.shp"))
# plot(st_geometry(sample_aoi))

sample_et0 <-
  rast(
    list.files(
      file.path(sample_dt_dir, "Seasonal_Water_Yield/ET0_monthly"),
      full.names = TRUE
    )
  )
   
sample_ppt <- 
  rast(
    list.files(
      file.path(sample_dt_dir, "Seasonal_Water_Yield/Precipitation_monthly"),
      full.names = TRUE
    )
  )

sample_sg <-
  rast(file.path(sample_dt_dir, "Seasonal_Water_Yield/soil_group_gura.tif"))


sample_rain_events_tbl <- 
  fread(
    file.path(sample_dt_dir, "Seasonal_Water_Yield/rain_events_gura.csv")
  )

# Can I use the wide format?
sample_rain_events_tbl_wide <-
  sample_rain_events_tbl %>%
  dcast(. ~ month, value.var = "events") %>%
  .[,c(as.character(1:12))]

write.csv(
  sample_rain_events_tbl_wide, 
  file.path(sample_dt_dir, "Seasonal_Water_Yield/rain_events_gura_wide.csv"), 
  row.names = FALSE
)


file.exists("/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/z_InVest_SampleData/Seasonal_Water_Yield/Precipitation_monthly/precip_gura_1.tif")

file.exists("/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/z_InVest_SampleData/Seasonal_Water_Yield/Precipitation_monthly/._precip_gura_1.tif")


junk_files <- list.files("Precipitation_monthly", pattern = "^\\._", full.names = TRUE)
file.remove(junk_files)