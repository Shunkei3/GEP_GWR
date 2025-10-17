#' ---
#' title: Functions to Prepare Inputs for InVEST SWY Model
#' ---
#' 

#| label: load-packages
library(here)
i_am("Scr/Code/InVEST_SWY/x_SWY_Functions.qmd")

library(data.table)
library(dplyr)
library(tidyr)

library(terra)
library(exactextractr)
library(sf)

library(tmap)
library(ggplot2)

# === Base path configurations === #
source(here("Scr/Code/0_path_config.R"))

#| label: test-aoi
#| eval: false
# 
# # /*===========================================*/
# #'= Skip!!!: This code chunk is just for testing the functions =
# # /*===========================================*/
# # /*===== AOI =====*/
# # Load Watersheds
# case_hybas_tbl <- readRDS(file.path(int_dir, "case_hybas_tbl.rds"))
# 
# tmp_region <- "na"
# tmp_region_hybas_sf <- 
#   filter(case_hybas_tbl, region_abbr == tmp_region) %>%
#   .[["tg_hybas_sf"]] %>%
#   .[[1]]
# 
# # === Select a watershed === #
# tmp_aoi_i <- tmp_region_hybas_sf[1,]
# 
# # st_crs(tmp_aoi_i)$epsg == espg_code
# 
# 
# aoi_sf_var <- tmp_aoi_i

#| label: setup

espg_code <- 6933 # WGS 84 / NSIDC EASE-Grid 2.0 Global
year <- 2020

# /*===========================================*/
#'=  Generic Functions =
# /*===========================================*/
mod_na_cells <- function(x) {
  x[is.nan(x)] <- NA  # make sure NaN -> NA
  x
}

# crop_mask_raster_for_aoi <- function(rast_object, bd_sf){	
# 	crop(
# 		rast_object, 
# 		vect(st_transform(bd_sf, crs = crs(rast_object, proj = TRUE))),
# 		mask = TRUE
#   )
# }

#| label: path-setup

# Create directory to save SWY inputs
# `fn_gen_save_dir` is defined in 0_path_config.R
swy_dt_dir <- 
  fn_gen_save_dir(
    where = int_dir,
    new_relative_dir_path = "SWY_inputs",
    return_path = TRUE
  )

# if(!exists("base_p_dir")){
#   stop("base_p_dir is not set. Please set the base data directory.")
# } else {
#   message("base_p_dir is set to: ", base_p_dir)
# }

# if(!exists("int_dir")){
#   stop("int_dir is not set. Please set the base data directory.")
# } else {
#   message("int_dir is set to: ", int_dir)
# }

#| label: prepare-inputs

fn_comp_prep_inputs <- 
  function(
    which_region, 
    which_aoi, 
    save_dir_return = FALSE, 
    overwrite = TRUE, 
    # If you want to test with a custom AOI, use the following two arguments
    custom_aoi_provided = FALSE, 
    custom_save_dir = NULL
  ){
  
  if(custom_aoi_provided){
    save_dir <- custom_save_dir
  } else {
    # Test Run
    # which_aoi <- tmp_aoi_i
    # which_aoi <- case_dt_r$hybas[[x]]
    # which_region <- which_aoi$continent 

    # /*===== Create a directory =====*/
    save_dir <-
      fn_gen_save_dir(
        where = swy_dt_dir, 
        new_relative_dir_path = file.path(which_region, which_aoi$HYBAS_ID),
        return_path = TRUE
      )
  }

  # /*===== Save AOI as a .gpkg file =====*/
  st_write(
    which_aoi,
    file.path(save_dir, "aoi.gpkg"),
    delete_dsn = TRUE,
    quiet = TRUE
  )

  # /*===== LULC (categorical) =====*/
  fn_get_lulc_aoi(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir
  )

  # /*===== DEM (continuous)=====*/
  fn_get_dem_aoi(
    region_abbr_var = which_region,
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir
  )
  # NOTE: I know that aoi_sf_var has region_abbr column, but for sake of generality, I pass region_abbr_var separately. Also, by doing so, I can use this function for custom AOI that does not have region_abbr column.

  # /*===== Soil Hydrologic Group (categorical) =====*/
  fn_get_soil_group_aoi(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir
  )

  # /*===== Biophysical Table (csv)=====*/
  fn_get_biophysical_table(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir
  )

  # /*===== ET0 (continuous) =====*/
  fn_get_monthly_et0(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir,
    type = "mean"
  )

  # /*===== Precipitation (continuous) =====*/
  # Monthly Alpha Table (csv) will be created within this function.
  fn_get_monthly_ppt(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir,
    type = "mean"
  )

  # /*===== Climate Zone (categorical) =====*/
  fn_get_cz(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir
  )

  # /*===== Rain Events Table (csv) =====*/
  # Create the table based on daily precipitation data and climate zone data
  fn_get_rain_events_tbl(
    aoi_sf_var = which_aoi,
    save_dir_var = save_dir,
    by_climate_zone = TRUE 
  )

  # /*===== Parameters: Threshold Flow Accumulation, Beta_i Parameter, Gamma Prameter (csv) =====*/ 
  # These parameter will be generated separately and saved as csv files. 

  if (save_dir_return){
    return(save_dir)
  }
}


fn_get_lulc_aoi <- function(aoi_sf_var, save_dir_var){
    # aoi_sf_var <- which_aoi
    # save_dir_var <- save_dir
  
  lulc_aoi <- 
    rast(file.path(raw_dir, "LULC/lulc_esa_2020.tif")) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE))),
      mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "near") %>%
    app(., mod_na_cells)

  # unique(values(lulc_aoi))
  # plot(lulc_aoi)

  # === Save === #
  writeRaster(
    lulc_aoi,
    here(save_dir_var, "lulc.tif"),
    overwrite = TRUE,
    datatype  = "INT1U",  # fits 0..255
    NAflag    = 255
  )
}


# lulc_test <- rast(here(save_dir, "lulc.tif"))
# unique(values(lulc_test))
# plot(lulc_test)



fn_get_dem_aoi <- function(region_abbr_var, aoi_sf_var, save_dir_var){
  # aoi_sf_var <- which_aoi
  # save_dir_var <- save_dir
  # region_abbr_var <- which_region

  dem_aoi <- 
    rast(file.path(raw_dir, "DEM", paste0("hyd_", region_abbr_var, "_dem_15s.tif"))) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE))),
      mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)
  
  # plot(dem_aoi)

  # === Save === #
  writeRaster(
    dem_aoi,
    file.path(save_dir_var, "dem.tif"),
    overwrite = TRUE
  )
}


fn_get_soil_group_aoi <- function(aoi_sf_var, save_dir_var) {

  # === Load Soil Group === #
  soil_group_aoi <- 
    rast(file.path(raw_dir, "Soil/hygeo_reclassed.tif")) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
      # mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "near") %>%
    app(., mod_na_cells)
  
  # plot(soil_group_aoi)
  # unique(values(soil_group_aoi))

  # === Save === #
  writeRaster(
    soil_group_aoi,
    here(save_dir_var, "soil_group.tif"),
    overwrite = TRUE,
    datatype  = "INT1U",  # fits 0..255
    NAflag    = 255
  )
}


# /*===== test =====*/
# soil_test <- rast(here(save_dir, "soil_group.tif"))
# unique(values(soil_test))


fn_get_biophysical_table <- function(aoi_sf_var, save_dir_var){

  # === Load LULC data === #
  # LULC data will be used to derive both Curve Number and Kc values
  lulc_aoi <- 
    rast(file.path(raw_dir, "LULC/lulc_esa_2020.tif")) %>%
    crop(., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))) %>%
    project(., paste0("epsg:", espg_code), method = "near")

  #/*--------------------------------*/
  #' ## Curve Number (CN) for each soil group
  #/*--------------------------------*/
  # /*===== Data Preparation =====*/
  # Load GCN250 data for curve number (250m resolution)
  gcn_aoi <- 
    rast(file.path(raw_dir, "Biophysical_Tables/GCN250/GCN250_ARCII.tif")) %>%
    crop(., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)

  # Load LULC data for the same AOI
  lulc_aoi_gcn <-
    project(lulc_aoi, gcn_aoi, method = "near") %>%
    app(., mod_na_cells)

  # Load Soil Hydrologic Group data for the same AOI
  soil_aoi_gcn <- 
    rast(file.path(raw_dir, "Soil/hygeo_reclassed.tif")) %>%
    project(., gcn_aoi, method = "near") %>%
    app(., mod_na_cells)

  # Stack the rasters
  stacked <- c(lulc_aoi_gcn, soil_aoi_gcn, gcn_aoi)
  names(stacked) <- c("lucode", "soil_group", "cn")

  # Mask to the AOI
  stacked_mask <- mask(stacked, vect(st_transform(aoi_sf_var, crs = crs(stacked, proj = TRUE))))

  # /*===== Curve Number Calculation =====*/
  cn_tbl_raw <- 
    data.table(values(stacked)) %>% 
    # na.omit(., cols = "lucode") %>% # remove rows with NA in "lucode"
    .[, soil_class := fcase(
      soil_group == 1, "A",
      soil_group == 2, "B",
      soil_group == 3, "C",
      soil_group == 4, "D"
    )]

  cn_tbl_out <- 
    cn_tbl_raw %>%
    # Aggregate mean CN
    .[, .(
      cn = round(mean(cn, na.rm = TRUE), 1)
      ), by = .(lucode, soil_class)
    ] %>%
    # For lucode == 210, insert 100 for CN value 
    .[lucode == 210, cn := 100] %>%
    dcast(lucode ~ soil_class, value.var = "cn")  %>%
    # Rename columns to match InVEST requirements
    setnames(
      old = intersect(c("A", "B", "C", "D"), names(.)), 
      new = paste0("CN_", intersect(c("A", "B", "C", "D"), names(.)))
    ) %>%
    setnames(names(.), tolower(names(.)))

  # /*===== Final Tuning =====*/
  # These values will not be used in the model.
  # 1. Make sure all required columns are present
  required_cn_cols <- c("cn_a", "cn_b", "cn_c", "cn_d")
  default_cn <- list(cn_a = 60, cn_b = 70, cn_c = 80, cn_d = 90)

  for (col in required_cn_cols) {
    # col = required_cn_cols[1]
    if (!col %in% names(cn_tbl_out)) {
      cn_tbl_out[, (col) := default_cn[[col]]]
    }
  }

  # 2. Make sure no missing values
  # Row-wise mean across available CN columns
  cn_tbl_out[, row_mean := round(rowMeans(.SD, na.rm = TRUE), 1), .SDcols = required_cn_cols]

  # Fill only missing cells in each row with that row's mean
  for (col in required_cn_cols) {
    cn_tbl_out[
      is.na(get(col)) & !is.na(row_mean),
      (col) := row_mean
    ]
  }

  # 3. Remove unnecessary columns
  cn_tbl_out <- cn_tbl_out[, c("lucode", required_cn_cols), with = FALSE]
  # Remove rows with NA in lucode (if NA exists, the model will fail)
  cn_tbl_out <- cn_tbl_out[!is.na(lucode)]

  #/*--------------------------------*/
  #' ## Kc values (kc_[MONTH])
  #/*--------------------------------*/
  ndvi_m_aoi <- 
    rast(file.path(raw_dir, paste0("NDVI/ndvi_m_", 2020, ".tif"))) %>%
    crop(., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))) %>%
    project(., lulc_aoi, method = "bilinear") %>%
    app(., mod_na_cells)

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

  kc_tbl_out <-
    rbind(
      monthly_kc_dt[!lucode %in% dt$lucode],
      dt[lucode %in% ls_lucode] # Remove class column
    )
  kc_cols <- grep("^Kc_", names(kc_tbl_out), value = TRUE)
  kc_tbl_out[, (kc_cols) := lapply(.SD, round, 2), .SDcols = kc_cols]


  #/*--------------------------------*/
  #' ## Combine CN and Kc tables to create biophysical table
  #/*--------------------------------*/ 
  biophysical_tbl <- merge(cn_tbl_out, kc_tbl_out, by = "lucode", all = TRUE)

  write.csv(
    biophysical_tbl, 
    here(save_dir_var, "biophysical_table.csv"),
    row.names = FALSE,
    na = "NaN" 
  )
}


fn_get_monthly_et0 <- function(aoi_sf_var, save_dir_var, type = "mean") {
  # type = "mean"  

  et0_m_aoi <- 
    rast(file.path(raw_dir, "Weather", paste0("et0_m_", type, "_Agromet_raw_", year, ".tif"))) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
      # mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)

  # === Save === #
  # Save each layer as a separate .tif file
  out_filenames <- paste0(names(et0_m_aoi), ".tif")

  out_file_paths <- 
    fn_gen_save_dir(
      where = save_dir_var,
      new_relative_dir_path = "ET0_monthly",
      return_path = TRUE
    )

  writeRaster(
    et0_m_aoi, 
    filename = file.path(out_file_paths, out_filenames), 
    overwrite = TRUE
  )
}


fn_get_monthly_ppt <- function(aoi_sf_var, save_dir_var, type = "mean"){
  # type = "mean"
  #/*--------------------------------*/
  #' ## Create monthly precipitation raster
  #/*--------------------------------*/
  ppt_m_aoi <- 
    rast(file.path(raw_dir, "Weather", paste0("p_m_", type, "_Agromet_raw_", year, ".tif"))) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
      # mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)

  # === Save === #
  # Save each layer as a separate .tif file
  out_filenames <- paste0(names(ppt_m_aoi), ".tif")

  out_file_paths <- 
    fn_gen_save_dir(
      where = save_dir_var,
      new_relative_dir_path = "Precipitation_monthly",
      return_path = TRUE
    )

  writeRaster(
    ppt_m_aoi, 
    filename = file.path(out_file_paths, out_filenames), 
    overwrite = TRUE
  )

  #/*--------------------------------*/
  #' ## Monthly Alpha Table
  #/*--------------------------------*/
  # Alpha: In practice, we suggest that for highly seasonal climates, alpha should be set to the antecedent monthly precipitation values, relative to the total precipitation: Pm-1/Pannual
  ppt_m_aoi_mask <- mask(ppt_m_aoi, aoi_sf_var)

  alpha_tbl <- 
    global(ppt_m_aoi_mask, "mean", na.rm = TRUE) %>%
    as.data.table(., keep.rownames = "month_id") %>%
    .[, month := as.integer(gsub("precip_", "", month_id))] %>%
    .[, alpha := 
      ifelse(
        month == 1, 1 / 12,
        data.table::shift(mean, type = "lag") / sum(.$mean, na.rm = TRUE)
    )]
  
  write.csv(
    alpha_tbl[, .(month, alpha)], 
    file.path(save_dir_var, "m_alpha.csv"), 
    row.names = FALSE
  )
}


fn_get_cz <- function(aoi_sf_var, save_dir_var){

  cz_aoi <- 
    rast(file.path(raw_dir, "ClimateZone/koppen_geiger_tif/1991_2020/koppen_geiger_0p1.tif")) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
      # mask = TRUE 
    ) %>%
    project(., paste0("epsg:", espg_code), method = "near") %>%
    app(., mod_na_cells)

  # === Save === #
  writeRaster(
    cz_aoi,
    filename = file.path(save_dir_var, "climate_zone.tif"),
    overwrite = TRUE,
    datatype  = "INT1U",  # fits 0..255
    NAflag    = 255
  )
}

# test_cz <- rast(here(save_dir, "climate_zone.tif"))
# unique(values(test_cz))
# plot(test_cz)


fn_get_rain_events_tbl <- function(aoi_sf_var, save_dir_var, by_climate_zone = TRUE){

  # /*===== Count Rain Events =====*/
  ppt_d_event <-
    rast(file.path(raw_dir, "Weather", paste0("p_d_Argomet_raw_", year, ".tif"))) %>%
    crop(
      ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
      # mask = TRUE
    ) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)
  
  ppt_d_event <- ppt_d_event > 0.1 # Logical raster
  
  # Define months for each layer
  dates <- seq.Date(as.Date(paste0(year, "-01-01")), by = "day", length.out = nlyr(ppt_d_event))
  months <- tolower(format(dates, "%b"))

  if(by_climate_zone){
    #/*--------------------------------*/
    #' ## Count Rain Events by Climate Zones in the AOI
    #/*--------------------------------*/
    cz_aoi_path <- file.path(save_dir_var, "climate_zone.tif")

    # If climate zone raster does not exist, create it
    if(!file.exists(cz_aoi_path)){
      fn_gen_climate_zone(
        aoi_sf_var = aoi_sf_var,
        save_dir_var = save_dir_var
      )
    }

    # Load Climate Zone raster
    cz_aoi <- rast(cz_aoi_path)
    # Unique climate zones in the AOI
    ls_cz_ids <- unique(cz_aoi, na.rm=TRUE)[[1]]

    res <- 
      data.table(
        cz_id = integer(),
        month = character(),
        rain_days = integer()
      )

    # Loop over climate zones
    for(cz_id in ls_cz_ids){
      # cz_id = ls_cz_ids[1]
      cz_i <- mask(cz_aoi, cz_aoi != cz_id, maskvalue = 1)

      for (m in unique(months)){
        # m = unique(months)[2]

        # Step 0: Get layers for the current month
        rain_month <- ppt_d_event[[which(months == m)]]
        # Step 1: Count rainy days at each pixel
        rain_pixel_count <- app(rain_month, fun = sum, na.rm = TRUE)
        # Step 2: Mask to current climate zone
        cz_mask_aligned <- resample(cz_i, rain_pixel_count, method = "near")
        rain_pixel_zone <- mask(rain_pixel_count, cz_mask_aligned)
        # plot(rain_pixel_zone)

        # Step 3: Aggregate all rainy days across pixels in zone
        rain_agg <- 
          global(rain_pixel_zone, "mean", na.rm = TRUE)[1, 1] %>%
          as.integer()

        # Store the result
        res <- 
          rbind(
            res, 
            data.table(
            cz_id = cz_id,
            month = m, 
            rain_days = as.integer(rain_agg))
          )
      }
    }

    res_return <- 
      dcast(res, cz_id ~ month, value.var = "rain_days") %>%
      .[, c("cz_id", tolower(month.abb)), with = FALSE]
    
    # === Save === #
    write.csv(
      res_return, 
      file.path(save_dir_var, "rain_events_by_cz.csv"), 
      row.names = FALSE
    )
    
  } else {
    #/*--------------------------------*/
    #' ## Count Rain Events for the whole AOI
    #/*--------------------------------*/
    res_return <-
      data.table(
        month = character(),
        events = integer()
      )
    
    for (m in unique(months)){
        # m = unique(months)[1]

        # Step 0: Get layers for the current month
        rain_month <- ppt_d_event[[which(months == m)]]
        # Step 1: Count rainy days at each pixel
        rain_pixel_count <- app(rain_month, fun = sum, na.rm = TRUE)
        # Step 2: Aggregate all rainy days across pixels in zone
        rain_agg <- 
          global(rain_pixel_count, "mean", na.rm = TRUE)[1, 1] %>%
          as.integer()

        # Store the result
        res_return <-
          rbind(
            res_return, 
            data.table(
              month = m, 
              events = as.integer(rain_agg)
            )
          )
    }
    res_return[, month := as.integer(match(tolower(month), tolower(month.abb)))]
    # === Save === #
    write.csv(
      res_return, 
      file.path(save_dir_var, "rain_events.csv"), 
      row.names = FALSE
    )
  }
}


# Check
# rain_event_tbl_sample <- fread()


get_beta_twi <- function(region_abbr_var, aoi_sf_var, save_dir_var){
  
  # /*===== Load Data =====*/
  dem_aoi <- 
      rast(file.path(raw_dir, "DEM", paste0("hyd_", region_abbr_var, "_dem_15s.tif"))) %>%
      crop(
        ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
        # mask = TRUE
      ) %>%
      project(., paste0("epsg:", espg_code), method = "bilinear") %>%
      app(., mod_na_cells)

  flow_acc <-
    rast(file.path(raw_dir, "FlowAccumulation", paste0("hyd_", region_abbr_var, "_acc_15s.tif"))) %>%
    crop(
        ., vect(st_transform(aoi_sf_var, crs = crs(., proj = TRUE)))
        # mask = TRUE
      ) %>%
    project(., paste0("epsg:", espg_code), method = "bilinear") %>%
    app(., mod_na_cells)

  # /*===== Slope Calculation =====*/
  slope_rad <- terrain(dem_aoi, v = "slope", unit = "radians")
    
  # Estimate cell area in m2
  res_m <- res(slope_rad)[1]  # assuming x and y resolution are equal
  cell_area <- res_m^2

  # /*===== Compute contributing area A (in m2) =====*/
  A <- flow_acc * cell_area

  # Add small epsilon to avoid division by zero
  epsilon <- 1e-6
  twi <- log((A + epsilon) / (tan(slope_rad) + epsilon))
  # Normalize TWI to [0, 1]
  twi_norm <- 
    (twi - global(twi, "min", na.rm = TRUE)[[1]]) / (global(twi, "max", na.rm = TRUE)[[1]] - global(twi, "min", na.rm = TRUE)[[1]])

  # plot(twi_norm)
  mean_beta_i <- global(twi_norm, "mean", na.rm = TRUE)[[1]]

  return(mean_beta_i)
}


fn_get_gamma <- function(aoi_sf_var, save_dir_var){

  # === Load soil hydrologic group data === #
  path_soil_group <- file.path(save_dir_var, "soil_group.tif")

  if(!file.exists(path_soil_group)){
    fn_get_soil_group_aoi(
      aoi_sf_var = aoi_sf_var,
      save_dir_var = save_dir_var
    )
  }

  soil_aoi <- rast(path_soil_group)

  # Count area of each soil group
  tab <- freq(soil_aoi, digits = 0)
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
}

