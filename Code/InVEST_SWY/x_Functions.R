if(!exists("base_dt_dir")){
  stop("base_dt_dir is not set. Please set the base data directory.")
} else {
  message("base_dt_dir is set to: ", base_dt_dir)
}

if(!exists("out_dir")){
  out_dir <- here("Data/Intermediate")
  message("out_dir is not set. Using default: ", out_dir)
} else {
  message("out_dir is set to: ", out_dir)
}


espg_code <- 6933 # WGS 84 / World Cylindrical Equal Area

year <- 2020 # Default year for the data

# /*===========================================*/
#'= Generic Functions =
# /*===========================================*/
crop_mask_raster_for_aoi <- function(rast_object, bd_sf){	
	crop(
		rast_object, 
		vect(st_transform(bd_sf, crs = crs(rast_object, proj = TRUE))),
		mask = TRUE
  )
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


# /*===========================================*/
#'= Save AOI as .gpkg file =
# /*===========================================*/
fn_save_aoi <- function(aoi){
  
  aoi_dir_path <- file.path(out_dir, "AOI", aoi$continent)
  
  if(!dir.exists(aoi_dir_path)){
    dir.create(aoi_dir_path, recursive = TRUE)
  }

  if(st_crs(aoi)$epsg != espg_code){
    aoi <- st_transform(aoi, crs = espg_code)
  }

  st_write(
    aoi, 
    file.path(aoi_dir_path, paste0(aoi$id, ".gpkg")),
    delete_dsn = TRUE,
    quiet = TRUE
  )
}

# /*===========================================*/
#'=  DEM =
# /*===========================================*/
fn_get_dem_aoi <- function(aoi) {
  # aoi = tmp_hybas_i

  region_abbr = aoi$continent

  # === Load DEM === #  
  dem <- 
    rast(file.path(base_dt_dir, "Raw/DEM", paste0("hyd_", region_abbr, "_dem_15s.tif"))) %>%
    crop_mask_raster_for_aoi(.,  aoi) %>%
    project(., paste0("epsg:", espg_code)) #  WGS 84 / World Cylindrical Equal Area 
    # project(., "ESRI:54030") #World Robinson projection
  
  # plot(st_geometry(aoi))
  # plot(dem)

  # === Save === #
  save_path <- file.path(here(out_dir, "DEM", region_abbr))
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  filename <- paste0(aoi$id, ".tif")

  writeRaster(
  	dem,
	  file.path(save_path, filename),
	  overwrite=TRUE
	)
}


# /*===========================================*/
#'= Land Use/Land Cover =
# /*===========================================*/
fn_get_lulc_aoi <- function(aoi) {
  # aoi = tmp_hybas_i
  
  region_abbr = aoi$continent

  # === Load LULC === #
  lulc_aoi <- 
    rast(file.path(base_dt_dir, "Raw/LULC/lulc_esa_2020.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi) %>%
    project(., paste0("epsg:", espg_code), method = "near") 
  
  lulc_aoi[is.na(lulc_aoi)] <- 999
    
  # unique(lulc_aoi, na.rm = FALSE)
  # plot(lulc_aoi)
  # all_values <- values(lulc_aoi) 
  # unique(all_values)

  # === Save === #
  save_path <- file.path(out_dir, "LULC", region_abbr)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  filename <- paste0(aoi$id, ".tif")

  writeRaster(
  	lulc_aoi,
	  file.path(save_path, filename),
	  overwrite=TRUE
  )

  # test <- rast(file.path(save_path, filename))
  # unique(test, na.rm = FALSE)
}


# /*===========================================*/
#'=  Precipitation =
# /*===========================================*/
#/*--------------------------------*/
#' ## Monthly Total Precipitation (mm/month, .tif)
#/*--------------------------------*/
# Twelve files, one for each month. File names must end with the month number (1-12). For example, the filenames ‘precip_1.tif’ and ‘precip1.tif’ are both valid names for the month of January.

# example_ppt <- 
#   rast(here("Data/z_InVest_SampleData/Seasonal_Water_Yield/Precipitation_monthly/precip_gura_1.tif"))

fn_get_month_ppt <- function(aoi, type = "mean"){
  # aoi = tmp_hybas_i
  # type = "mean" # or "total"
  
  region_abbr = aoi$continent

  # === Load Monthly Total Precipitation === #  
  which_p <- paste0("p_m_", type, "_Agromet_raw_2020.tif")
  
  p_monthly_aoi <- 
    rast(file.path(base_dt_dir, "Raw/Weather", which_p)) %>%
    crop_mask_raster_for_aoi(.,  aoi) %>%
    project(., paste0("epsg:", espg_code))

  # plot(p_monthly_aoi)

  # === Save Rasters === #
  save_path <- file.path(out_dir, "Precip", region_abbr, aoi$id)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  output_filenames <- paste0(names(p_monthly_aoi), ".tif")
  
  writeRaster(
    p_monthly_aoi, 
    filename = file.path(save_path, output_filenames), 
    overwrite = TRUE
  )

  # === Monthly Alpha Table === #
  # Alpha: In practice, we suggest that for highly seasonal climates, alpha should be set to the antecedent monthly precipitation values, relative to the total precipitation: Pm-1/Pannual

  alpha_tbl <- 
    global(p_monthly_aoi, "mean", na.rm = TRUE) %>%
    as.data.table(., keep.rownames = "month_id") %>%
    .[, month := as.integer(gsub("precip_", "", month_id))] %>%
    .[, alpha := 
      ifelse(
        month == 1, 1 / 12,
        data.table::shift(mean, type = "lag") / sum(alpha_tbl$mean, na.rm = TRUE)
    )]
  
  save_path <- file.path(out_dir, "Alpha_m", region_abbr)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }
  write.csv(
    alpha_tbl[, .(month, alpha)], 
    file.path(save_path, paste0(aoi$id, ".csv")), 
    row.names = FALSE
  )
}

#/*--------------------------------*/
#' ## Rain Events Table (.csv)
#/*--------------------------------*/
# the rain events table has a count of the number of rain events per month, where each rain event exceeds 0.1mm. 

# example_rain_tbl <- 
#   fread(here("Data/z_InVest_SampleData/Seasonal_Water_Yield/rain_events_gura.csv")) 

# example_rain_tbl_climate <-
#   fread(here("Data/z_InVest_SampleData/Seasonal_Water_Yield/climate_zone_table_gura.csv"))


# https://community.naturalcapitalproject.org/t/rain-events-table/5488/2
# cz <- 
#   rast(file.path(base_dt_dir, "Raw/ClimateZone/koppen_geiger_climatezones_1991_2020_1km.tif")) %>%
#   crop_mask_raster_for_aoi(.,  aoi)

# example_cz <-
#   rast(here("Data/z_InVest_SampleData/Seasonal_Water_Yield/climate_zones_gura.tif"))
# plot(example_cz)

get_rain_events_tbl_cz <- function(aoi){
  # aoi = tmp_hybas_i
  
  region_abbr = aoi$continent

  # === Load Climate Zone === #
  cz <- 
    rast(file.path(base_dt_dir, "Raw/ClimateZone/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi)

  # plot(cz)
  
  # --- Save cz --- #
  cz_save_path <- file.path(out_dir, "ClimateZone", region_abbr)
  if(!dir.exists(cz_save_path)){
    dir.create(cz_save_path, recursive = TRUE)
  }

  writeRaster(
    project(cz, paste0("epsg:", espg_code), method = "near"), 
    file.path(cz_save_path, paste0(aoi$id, ".tif")), 
    overwrite = TRUE
  )

  # === Load Daily Precipitation === #
  d_rain <- 
    rast(file.path(base_dt_dir, "Raw/Weather/p_d_Argomet_raw_2020.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi)

  # Create a logical raster for rain events
  d_rain_event <- d_rain > 0.1

  # Define months for each layer
  dates <- seq.Date(as.Date(paste0(year, "-01-01")), by = "day", length.out = nlyr(d_rain))  # adjust year
  months <- tolower(format(dates, "%b"))

  # Get unique climate zones within the AOI
  cz_vals <- unique(na.omit(values(cz)))

  res <- 
    data.table(
      cz_id = integer(), 
      month = character(), 
      rain_days = integer()
    )

  # plot(cz)
  # plot(mask(cz, cz != 3, maskvalues = 1) )
    
  for (cz_i in cz_vals) {
    # cz_i <- cz_vals[1]
    # Create a mask for the current climate zone
    cz_mask <- mask(cz, cz != cz_i, maskvalues = 1) 
    # plot(cz_mask)

    for(m in unique(months)){
      # m <- unique(months)[1]
      # Get layers for the current month
      rain_month <- d_rain_event[[which(months == m)]] 

      # Step 1: Count rainy days at each pixel
      rain_pixel_count <- app(rain_month, fun = sum, na.rm = TRUE)

      # Step 2: Mask to current climate zone
      cz_mask_aligned <- resample(cz_mask, rain_pixel_count, method = "near")
      rain_pixel_zone <- mask(rain_pixel_count, cz_mask_aligned)
      # plot(rain_pixel_zone)

      # Step 3: Aggregate all rainy days across pixels in zone
      rain_agg <- 
        global(rain_pixel_zone, "mean", na.rm = TRUE)[1, 1]
      
      # Store result
      res <- 
        rbind(
          res, 
          data.table(
            cz_id = cz_i,
            month = m, 
            rain_days = as.integer(rain_agg))
          )
    }
  }

  res_wide <- 
    # copy(res)[, month := match(month, tolower(month.abb))] %>% # No
    dcast(res, cz_id ~ month, value.var = "rain_days") %>%
    .[, c("cz_id", tolower(month.abb)), with = FALSE]

  # === Save === #
  save_path <- file.path(out_dir, "RainEvents", region_abbr)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  write.csv(
    res_wide, 
    file.path(save_path, paste0(aoi$id, ".csv")), 
    row.names = FALSE
  )
}

# Loop through each climate zone

# /*===========================================*/
#'=  ET0 =
# /*===========================================*/
# NOTE: the “depth” of water that evapotranspirates from a given region; it is not “per” pixel, square meter, or any other area unit.
# NOTE: ET data is based on the same precipitation data that is being used as a model input.

# NOTE: You could use the ‘modified Hargreaves’ method.


#/*--------------------------------*/
#' ## ET0 Monthly, (mm/month, .tif)
#/*--------------------------------*/
# example_et0 <- 
#   rast(here("Data/z_InVest_SampleData/Seasonal_Water_Yield/ET0_monthly/ET0_gura_1.tif"))

fn_get_month_et0 <- function(aoi, type = "mean"){
  # aoi = tmp_hybas_i
  # type = "mean" 
  region_abbr = aoi$continent

  # === Load Monthly ET0 for aoi === #
  which_et0 <- paste0("et0_m_", type, "_Agromet_raw_2020.tif")
  et0_aoi <- 
    rast(file.path(base_dt_dir, "Raw/Weather", which_et0)) %>%
    crop_mask_raster_for_aoi(.,  aoi) %>%
    project(., paste0("epsg:", espg_code))

  # plot(et0_aoi)
  # === Save === #
  save_path <- file.path(out_dir, "ET0", region_abbr, aoi$id)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  output_filenames <- paste0(names(et0_aoi), ".tif")
  
  writeRaster(
    et0_aoi, 
    filename = file.path(save_path, output_filenames), 
    overwrite = TRUE
  )
}

# /*===========================================*/
#'=  Soil =
# /*===========================================*/
# For details on the creation of `hygeo_reclassed.tif`, see `Prep_Soil.R`.


fn_get_soil_group_aoi <- function(aoi) {
  # aoi ƒ= tmp_hybas_i
  region_abbr = aoi$continent

  # === Load Soil Group === #
  soi_aoi_group <- 
    rast(file.path(base_dt_dir, "Raw/Soil/hygeo_reclassed.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi) %>%
    project(., paste0("epsg:", espg_code), method = "near") 
  
  # plot(soi_aoi_group)

  # === Save === #
  save_path <- file.path(out_dir, "SoilGroup", region_abbr)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  writeRaster(
    soi_aoi_group,
    file.path(save_path, paste0(aoi$id, ".tif")),
    overwrite = TRUE
  )
}

# /*===========================================*/
#'=  Flow Accumulation =
# /*===========================================*/
#' ...
#' See "Prep_FlowAcc.R" for details on how to prepare the flow accumulation data.





# /*===========================================*/
#'=  Derive beta parameter as TWI =
# /*===========================================*/
#' TWI = ln(A/tan(slope))
#' A: flow accumulation x cell area (in m2)
#' slope: slope in radians


#' Objective: Given a watershed, calculate the beta parameter as TWI (Topographic Wetness Index).


library(dplyr)
library(data.table)
library(terra)


get_beta_twi <- function(aoi){
  region_abbr = aoi$continent
  file_path <- file.path(save_path, paste0(aoi$id, ".tif"))

  if(file.exists(file_path)){

    twi_norm <- rast(file_path)

  } else {
    # /*===== Load Data =====*/
    dem <- 
      rast(file.path(base_dt_dir, "Raw/DEM", paste0("hyd_", region_abbr, "_dem_15s.tif"))) %>%
      crop_mask_raster_for_aoi(., aoi) %>%
      project(., paste0("epsg:", espg_code))
    
    flow_acc <-
      rast(file.path(base_dt_dir, "Raw/FlowAccumulation", paste0("hyd_", region_abbr, "_acc_15s.tif"))) %>%
      crop_mask_raster_for_aoi(., aoi) %>%
      project(., paste0("epsg:", espg_code))

    # /*===== Slope Calculation =====*/  
    slope_rad <- terrain(dem, v = "slope", unit = "radians")
    
    # Estimate cell area in m2
    res_m <- res(slope_rad)[1]  # assuming x and y resolution are equal
    cell_area <- res_m^2

    # === Compute contributing area A (in m2) === #
    A <- flow_acc * cell_area

    # Add small epsilon to avoid division by zero
    epsilon <- 1e-6
    twi <- log((A + epsilon) / (tan(slope_rad) + epsilon))
    # Normalize TWI to [0, 1]
    twi_norm <- 
      (twi - global(twi, "min", na.rm = TRUE)[[1]]) / (global(twi, "max", na.rm = TRUE)[[1]] - global(twi, "min", na.rm = TRUE)[[1]])

    # plot(twi_norm)

    # === Save === #
    save_path <- file.path(out_dir, "Beta_i", region_abbr)
    if(!dir.exists(save_path)){
      dir.create(save_path, recursive = TRUE)
    }

    writeRaster(
      twi_norm,
      file_path,
      overwrite = TRUE
    )
  }

  # === Return === #
  mean_beta_i <- global(twi_norm, "mean", na.rm = TRUE)[[1]]

  return(mean_beta_i)
}


# /*===========================================*/
#'=  Derive Gamma parameter based on soil type =
# /*===========================================*/
# 1: high infiltration rate (A)
# 2: moderate infiltration rate (B)
# 3: slow infiltration rate (C)
# 4: very slow infiltration rate (D)

fn_get_gamma <- function(aoi){

  region_abbr = aoi$continent

  file_path <-  file.path(out_dir, "SoilGroup", region_abbr, paste0(aoi$id, ".tif"))

  # === Load soil hydrologic group data === #
  if(file.exists(file_path)){
    soil_aoi_group <- rast(file_path)
  } else {
    fn_get_soil_group_aoi(aoi)
    soil_aoi_group <- rast(file_path)
  }

  # plot(soil_aoi_group)
  
  # Count area of each soil group
  tab <- freq(soil_aoi_group, digits = 0)
  # Compute area fractions
  tab$area_fraction <- tab$count / sum(tab$count)

  gamma_lookup <- 
    data.frame(
      value = c(1, 2, 3, 4),
      gamma = c(0.95, 0.85, 0.75, 0.55)
    )

  # Merge lookup and compute weighted average
  tab <- merge(tab, gamma_lookup, by = "value")
  gamma_weighted <- sum(tab$gamma * tab$area_fraction)

  return(gamma_weighted)
}


# /*===========================================*/
#'= Biophysical Table =
# /*===========================================*/
# A table mapping each LULC code to biophysical properties of the corresponding LULC class. All values in the LULC raster must have corresponding entries in this table.

# A .csv (Comma Separated Value) table containing model information corresponding to each of the land use classes in the LULC raster. All LULC classes in the LULC raster MUST have corresponding values in this table. 

# Each row is a land use/land cover class and columns must be named and defined as follows:

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
#' ## GCN250 for Curve Number (cn_[SOIL_GROUP])
#/*--------------------------------*/
fn_get_cn_tbls <- function(aoi, lulc_aoi, soil_aoi){
  # aoi = tmp_hybas_i

  gcn_aoi <- 
    rast(file.path(base_dt_dir, "Raw/Biophysical_Tables/GCN250/GCN250_ARCII.tif")) %>%
    crop_mask_raster_for_aoi(., aoi)

  # Resample lulc_aoi to match GCN resolution and extent
  lulc_resampled <- resample(lulc_aoi, gcn_aoi, method = "near")
  lulc_resampled[is.na(lulc_resampled)] <- 999
  
  # unique(lulc_aoi, na.rm = FALSE)
  # unique(lulc_resampled, na.rm = FALSE)

  # Align soil_aoi to match gcn_aoi (same extent/resolution)
  soil_resampled <- resample(soil_aoi, gcn_aoi, method = "near")

  # unique(soil_resampled, na.rm = FALSE)

  # Stack the rasters
  stacked <- 
    c(
      lulc_resampled, 
      soil_resampled, 
      gcn_aoi
    ) 
  names(stacked) <- c("lucode", "soil_group", "cn")

  # unique(stacked['lucode'], na.rm = FALSE)

  # Extract values and remove NA
  
  # values_mat <- values(stacked)

  # unique(na.omit(values_mat[, "lucode"]))
  
  # unique(vals$lucode)
  # unique(na.omit(values_dt$lucode))
  
  values_dt <- 
    data.table(values(stacked)) %>%
    # na.omit(., cols = "lucode") %>% # Do not remove NA values in lucode!
    .[, soil_class := fcase(
      soil_group == 1, "A",
      soil_group == 2, "B",
      soil_group == 3, "C",
      soil_group == 4, "D"
    )] %>%
    # Aggregate mean CN
    .[, .(
      cn = round(mean(cn, na.rm = TRUE), 1)
      ), by = .(lucode, soil_class)
    ]%>%
    # For lucode == 210, insert 100 for CN value 
    .[lucode == 210, cn := 100] %>%
    dcast(lucode ~ soil_class, value.var = "cn")  %>%
    # Rename columns to match InVEST requirements
    setnames(
      old = intersect(c("A", "B", "C", "D"), names(.)), 
      new = paste0("CN_", intersect(c("A", "B", "C", "D"), names(.)))
    ) %>%
    setnames(names(.), tolower(names(.)))
  
    required_cn_cols <- c("cn_a", "cn_b", "cn_c", "cn_d")
    default_cn <- list(cn_a = 60, cn_b = 70, cn_c = 80, cn_d = 90)

    for (col in required_cn_cols) {
      # col = required_cn_cols[1]
      if (!col %in% names(values_dt)) {
        values_dt[, (col) := default_cn[[col]]]
      }
    }

    # Replace all NA values with a default numeric value (e.g., 0)
    values_dt[is.na(cn_a), cn_a := 60]
    values_dt[is.na(cn_b), cn_b := 70]
    values_dt[is.na(cn_c), cn_c := 80]
    values_dt[is.na(cn_d), cn_d := 90]
    
    # Remove 
    values_dt <- values_dt[, !names(values_dt) %in% "na", with = FALSE]
  
    return(values_dt)
}

# test1 <- get_cn_tbls(tmp_hybas_i)
# test2 <- get_cn_tbls(tmp_hybas_j)

#/*--------------------------------*/
#' ## Kc values (kc_[MONTH])
#/*--------------------------------*/
#' It is assumed that the Kc values are constant for some LULC classes across all months. 
#' Those classes are evergreen forest, deciduous forest, shrublands, urban, bare ground, snow/ice, and open water.
# === Forests: FLUXNET shows evergreen Kc ~0.9 year-round === #
#'  + https://community.naturalcapitalproject.org/t/seasonal-water-yield-limitations/79/9
#'  + https://community.naturalcapitalproject.org/t/kc-crop-coefficient-for-urban-and-mining-land-cover-classes/1555
#'  + https://farmwest.com/climate/calculator-information/et/crop-coefficients/
#'  + https://ucanr.edu/site/center-landscape-urban-horticulture/turfgrass-crop-coefficients-kc
#'  + @Liu.etal2017: "Environmental controls on seasonal ecosystem evapotranspiration/potential evapotranspiration ratio as determined by the global eddy flux measurements"
#'  + FAO 
# === Grass/shrub: mean Kc ~0.7 monthly === #
# + @Liu.etal2017: "Environmental controls on seasonal ecosystem evapotranspiration/potential evapotranspiration ratio as determined by the global eddy flux measurements"
# === Open water: FAO recommends ~0.75 average === #
# === Urban: blended impervious + vegetation, ~0.475 average === #
# === Bare ground and snow/ice: minimal ET, so low Kc === #


fn_get_kc_tbls <- function(aoi, lulc_aoi){
  # aoi = tmp_hybas_i

  region_abbr = aoi$continent

  # === Load NDVI Values === #
  ndvi_m_aoi <- 
    rast(file.path(base_dt_dir, paste0("Raw/NDVI/ndvi_m_", year, ".tif"))) %>%
    crop_mask_raster_for_aoi(., aoi)

  # Reproject NDVI
  ndvi_m_aoi <- project(ndvi_m_aoi, lulc_aoi)
  # Resample the projected NDVI to match the resolution and lightment of LULC raster
  # 2. Resample the projected NDVI stack to match the resolution and alignment of LULC
  ndvi_m_aoi <- resample(ndvi_m_aoi, lulc_aoi, method = "bilinear")
  
  # unique(lulc_aoi, na.rm = FALSE)
  
  # lulc_aoi
  kc_m_aoi <- 1.37 * ndvi_m_aoi - 0.069
  kc_m_aoi <- clamp(kc_m_aoi, lower = 0, upper = 1.2)

  # Extract zonal mean Kc per LULC and month
  monthly_kc_dt <- 
    lapply(
      1:12, 
      function(i) {
        # i = 1
        kc_i <- kc_m_aoi[[i]]
        stats <- 
          zonal(kc_i, lulc_aoi, fun = "mean", digits = 3, na.rm = TRUE) %>%
          data.table()
        
        stats[, month := i] %>%
        setnames(names(.), c("lucode", "Kc", "month"))
      }
    ) %>%
    rbindlist(., use.names = TRUE) %>%
    dcast(lucode ~ month, value.var = "Kc") %>%
    setnames(
      as.character(1:12),
      paste0("Kc_", 1:12)
    )

  # For crops, Kc = 1.37 * NDVI - 0.069

  # /*===== Base Kc_table =====*/
  dt <- 
    data.table(
      lucode = c(50, 60, 70, 190, 200, 210, 220),
      class = c("Evergreen forest","Deciduous forest","Shrublands","Urban","Bare ground","Snow/Ice","Open water")
    )
  
  dt[, paste0("Kc_", 1:12) := {
      Kc <- 
        switch(
          class[1],
          "Evergreen forest" = rep(0.90, 12),
          "Deciduous forest" = rep(0.85, 12),  # treating deciduous similarly
          "Shrublands" = rep(0.70, 12),
          "Urban" = rep(0.475, 12), #Impervious surfaces=0, Vegerated patches = 0.7 ~ 0.9, so average is 0.475
          "Bare ground" = rep(0.20, 12),
          "Snow/Ice" = rep(0, 12),
          "Open water" = rep(1.00, 12)
        )
      as.list(Kc)
    }, by = lucode
  ] %>%
    .[, class := NULL] # Remove class column

  ls_lucode <- unique(monthly_kc_dt$lucode)
  
  out_kc_tbl <-
    rbind(
      monthly_kc_dt[!lucode %in% dt$lucode],
      dt[lucode %in% ls_lucode] # Remove class column
    )

  return(out_kc_tbl)
}



#/*--------------------------------*/
#' ## Generate Biophysical Table (.csv)
#/*--------------------------------*/

fn_gen_biophysical_table <- function(aoi){

  # aoi = tmp_hybas_i

  region_abbr = aoi$continent

  # === Load Data === #
  lulc_aoi <- 
    rast(file.path(base_dt_dir, "Raw/LULC/lulc_esa_2020.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi) 
  lulc_aoi[is.na(lulc_aoi)] <- 999
  # unique(lulc_aoi, na.rm = FALSE)

  soil_aoi <- 
    rast(file.path(base_dt_dir, "Raw/Soil/hygeo_reclassed.tif")) %>%
    crop_mask_raster_for_aoi(.,  aoi)
  
  # unique(soil_aoi)

  # === Get Curve Number Table === #
  cn_tbl <- 
    fn_get_cn_tbls(aoi, lulc_aoi, soil_aoi) %>%
    setnames(names(.), tolower(names(.))) 

  # === Get Kc values === #
  kc_tbl <- fn_get_kc_tbls(aoi, lulc_aoi)

  # === Combine Tables === #
  biophysical_tbl <- 
    merge(cn_tbl, kc_tbl, by = "lucode", all = TRUE)

  # === Save CSV === #
  save_path <- file.path(out_dir, "Biophysical_Tables", region_abbr)
  if(!dir.exists(save_path)){
    dir.create(save_path, recursive = TRUE)
  }

  write.csv(
    biophysical_tbl, 
    file.path(save_path, paste0(aoi$id, ".csv")), 
    row.names = FALSE,
    na = "NaN" 
  )
}





# I need to extract the curve number values for each LULC class in the LULC raster.


# /*===========================================*/
#'=  Comprehensive =
# /*===========================================*/
# fn_comp_prep_data <- function(aoi) {
#   # aoi = tmp_hybas_i
  
#   id = aoi$id
#   region_abbr = aoi$continent
#   espg <- 6933

#   # === Save AOI === #
#   aoi_dir_path <- file.path(out_dir, "AOI", region_abbr)
#   if(!dir.exists(aoi_dir_path)){
#     dir.create(aoi_dir_path, recursive = TRUE)
#   }  
  
#   st_write(
#     st_transform(aoi, crs = espg), 
#     file.path(aoi_dir_path, paste0(aoi$id, ".gpkg")),
#     delete_dsn = TRUE,
#     quiet = TRUE
#   )

#   # === DEM === #
#   fn_get_dem_aoi(aoi)


#   # === Flow accumulation === #  



#   # === LULC === #

# }
