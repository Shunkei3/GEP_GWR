library(data.table)
library(dplyr)
library(terra)
library(sf)
library(tmap)
library(ggplot2)


base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"


#' I need to get the daily precipitation data to calculate monthly rain events (>= 0.1mm) 
#' The daily precipitation data is aggregated to monthly values,
#' Also, Monthly ET0 data is derived based on the same precipitation data. Use the modified Hargreaves method.
#' 



# CPC Global Unified Gauge-Based Analysis of Daily Precipitation (https://psl.noaa.gov/data/gridded/data.cpc.globalprecip.html)



#' Global high-resolution precipitation: MSWEP
#' https://www.gloh2o.org/mswep/
#' 0.1 degree resolution
#' "Precipitation estimates from MSWEP are thus likely more accurate than those from MSWX-Past in (a) densely gauged regions and (b) convection-dominated regions due to the satellite data."
#' I downloaded the daily precipitation data for 2020.


#' Multi-Source Weather (MSWX) 
#' https://www.gloh2o.org/mswx/
#' See: https://journals.ametsoc.org/view/journals/bams/103/3/BAMS-D-21-0145.1.xml
#' 0.1 degree resolution
#' Tmin (degrees Celsius)
#' Tmax (degrees Celsius)
#' SWd: Shortwave Downward Radiation (W m-2)
#' LWd: Longwave Downward Radiation (W m-2)
# !Radiation data is not available!
# How can I get the monthly 


#' (Currently in use) AgERA5, daily
#' 0.1 degree resolution
#' https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=overview
#' P (Preciptation_Flux, mm/day): Total precipitation
#' ET0 (mm/day): Penman-Monteith reference evapotranspiration according to the FAO56 approach


test <- rast("/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Raw/Weather/Agromet_Temp_d_2020_Jan_Mar/Temperature-Air-2m_Max-24h_C3S-glob-agric_AgERA5_20200501_final-v2.0.0.nc")


# /*===========================================*/
#'=  Daily Precipitation Data (MSWEP) =
# /*===========================================*/
# === Daily Precipitation (CPU) === #
# ppt_d_cpu_raw <- rast('/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Raw/Precip/precip.2019.nc')

# === Daily Precipitation (MSWEP) === #
ls_mswp_2020 <-
  list.files(
    file.path(base_dt_dir, "Raw/Weather/MSWEP_2020"),
    full.names = FALSE
  )

ppt_d_mswep_raw <- 
  lapply(
    ls_mswp_2020, 
    function(x) {
      # x = ls_mswp_2020[1]
      layer <- rast(file.path(base_dt_dir, "Raw/Weather/MSWEP_2020", x))
      
      names(layer) <-  paste0("d", as.numeric(gsub("\\.nc$", "", x)) - 2020000)

      return(layer)
    }
  ) %>%
  rast()

writeRaster(
  ppt_d_mswep_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/ppt_d_mswep_raw_2020.tif"), 
  overwrite = TRUE
)

# /*===========================================*/
#'= Daily Precipitation (AgERA5) =
# /*===========================================*/
ls_files_p_d_Agromet_2020 <- 
  list.files(
    file.path(base_dt_dir, "Raw/Weather/Agromet_P_d_2020"),
    full.names = FALSE
  )

p_d_Agromet_raw <-
  lapply(
    ls_files_p_d_Agromet_2020, 
    function(x) {
      # x = ls_files_p_d_Agromet_2020[10]
      layer <- rast(file.path(base_dt_dir, "Raw/Weather/Agromet_P_d_2020", x))

      # plot(layer)
      
     #  Layer names
      day <- gsub(".*?(\\d{8}).*", "\\1", x)
      names(layer) <- paste0("d", as.integer(as.Date(day, format = "%Y%m%d") - as.Date(paste0(substr(day, 1, 4), "-01-01")) + 1))
      
      return(layer)
    }
  ) %>%
  rast()

writeRaster(
  p_d_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/p_d_Argomet_raw_2020.tif"), 
  overwrite = TRUE
)

#/*--------------------------------*/
#' ## Monthly Average and Total Precipitation
#/*--------------------------------*/
ls_dates <- gsub(".*?(\\d{8}).*", "\\1", ls_files_p_d_Agromet_2020)
ls_dates_months <- month(as.Date(ls_dates, format = "%Y%m%d"))


ls_p_m_mean_Agromet_raw <- list(rep(NA, length(unique(ls_dates_months))))

ls_p_m_sum_Agromet_raw <- list(rep(NA, length(unique(ls_dates_months))))


for (m in unique(ls_dates_months)) {

  # m = 1
  tg_file_names <- ls_files_p_d_Agromet_2020[ls_dates_months == m]

  layers_d_m <- rast(file.path(base_dt_dir, "Raw/Weather/Agromet_P_d_2020", tg_file_names))

  # Aggregate daily layers to monthly mean
  monthly_mean <- mean(layers_d_m, na.rm = TRUE)
  names(monthly_mean) <- paste0("precip_", m)
  ls_p_m_mean_Agromet_raw[[m]] <- monthly_mean

  # Aggregate daily layers to monthly total
  monthly_total <- sum(layers_d_m, na.rm = TRUE)
  names(monthly_total) <- paste0("precip_", m)
  ls_p_m_sum_Agromet_raw[[m]] <- monthly_total
}

p_m_mean_Agromet_raw <- rast(ls_p_m_mean_Agromet_raw)
p_m_sum_Agromet_raw <- rast(ls_p_m_sum_Agromet_raw)

writeRaster(
  p_m_mean_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/p_m_mean_Agromet_raw_2020.tif"), 
  overwrite = TRUE
)

writeRaster(
  p_m_sum_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/p_m_sum_Agromet_raw_2020.tif"), 
  overwrite = TRUE
)


# /*===========================================*/
#'=  Daily ET0 (AgERA5) =
# /*===========================================*/
ls_files_et0_d_Agromet_2020 <- 
  list.files(
    file.path(base_dt_dir, "Raw/Weather/Agromet_ET0_d_2020"),
    full.names = FALSE
  )

et0_d_Agromet_raw <-
  lapply(
    ls_files_et0_d_Agromet_2020, 
    function(x) {
      # x = ls_files_et0_d_Agromet_2020[10]
      layer <- rast(file.path(base_dt_dir, "Raw/Weather/Agromet_ET0_d_2020", x))

      # plot(layer)
      
     #  Layer names
      day <- gsub(".*?(\\d{8}).*", "\\1", x)
      names(layer) <- paste0("d", as.integer(as.Date(day, format = "%Y%m%d") - as.Date(paste0(substr(day, 1, 4), "-01-01")) + 1))
      
      return(layer)
    }
  ) %>%
  rast()

writeRaster(
  et0_d_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/et0_d_Agromet_raw_2020.tif"), 
  overwrite = TRUE
)


#/*--------------------------------*/
#' ## Monthly Average and Total ET0
#/*--------------------------------*/
ls_dates <- gsub(".*?(\\d{8}).*", "\\1", ls_files_et0_d_Agromet_2020)
ls_dates_months <- month(as.Date(ls_dates, format = "%Y%m%d"))


ls_et0_m_mean_Agromet_raw <- list(rep(NA, length(unique(ls_dates_months))))

ls_et0_m_sum_Agromet_raw <- list(rep(NA, length(unique(ls_dates_months))))


for (m in unique(ls_dates_months)) {

  # m = 1
  tg_file_names <- ls_files_et0_d_Agromet_2020[ls_dates_months == m]

  layers_d_m <- rast(file.path(base_dt_dir, "Raw/Weather/Agromet_ET0_d_2020", tg_file_names))

  # Aggregate daily layers to monthly mean
  monthly_mean <- mean(layers_d_m, na.rm = TRUE)
  names(monthly_mean) <- paste0("ET0_", m)
  ls_et0_m_mean_Agromet_raw[[m]] <- monthly_mean

  # Aggregate daily layers to monthly total
  monthly_total <- sum(layers_d_m, na.rm = TRUE)
  names(monthly_total) <- paste0("ET0_", m)
  ls_et0_m_sum_Agromet_raw[[m]] <- monthly_total
}

et0_m_mean_Agromet_raw <- rast(ls_et0_m_mean_Agromet_raw)
et0_m_sum_Agromet_raw <- rast(ls_et0_m_sum_Agromet_raw)

writeRaster(
  et0_m_mean_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/et0_m_mean_Agromet_raw_2020.tif"), 
  overwrite = TRUE
)

writeRaster(
  et0_m_sum_Agromet_raw, 
  filename = file.path(base_dt_dir, "Raw/Weather/et0_m_sum_Agromet_raw_2020.tif"), 
  overwrite = TRUE
)
