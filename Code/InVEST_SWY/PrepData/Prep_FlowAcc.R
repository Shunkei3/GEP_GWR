# /*===========================================*/
#'= Create a Flow Accumulation Table =
# /*===========================================*/
# Threshold Flow Accumulation (number, units: number of pixels, required): The number of upslope pixels that must flow into a pixel before it is classified as a stream.


library(here)
library(data.table)
library(dplyr)
library(lubridate)
library(terra)
library(sf)
library(ggplot2)

library(foreach)

base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_interm_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"


if(!dir.exists(file.path(out_interm_dir, "FlowAcc"))){
  dir.create(file.path(out_interm_dir, "FlowAcc"), recursive = TRUE)
}

# world_bd <- 
#   st_read(
#     file.path(base_dt_dir, "Raw/World_Data_Grids/country.shp"),
#     quiet = TRUE
#   ) %>%
#   setnames(names(.), tolower(names(.)))
# st_crs(world_bd) <- 4326
# world_bd <- st_transform(world_bd, crs = "epsg:6933") # WGS 84 / World

# === Load Watersheds === #
case_hybas_tbl <- readRDS(file.path(out_interm_dir, "case_hybas_tbl.rds"))


# /*===========================================*/
#'=  Define Functions =
# /*===========================================*/
get_flow_acc_values <- function(aoi, river_raw, flow_acc_raw){
  # aoi = filter(case_hybas_r, id == "na_7060000010")
  # aoi = filter(case_hybas_r, id == "au_5060077570")
  # aoi = filter(case_hybas_r, id == "au_5060081550") #fail

  # plot(st_geometry(aoi))
  
  river_aoi <- try(
    st_intersection(
      river_raw, 
      aoi
    ) %>%
    st_buffer(., dist = 100),
    silent = TRUE
  )

  # ggplot() +
  #   geom_sf(data = aoi) +
  #   geom_sf(data = river_aoi, fill = "blue")
  
  if(inherits(river_aoi, "try-error")){
    res_values <- NA    
  } else {
    if(nrow(river_aoi) == 0){
    # If no rivers in the AOI, return NA values
    res_values <- NA    
    } else {
      # Extract values
      flow_acc_values <- 
        exactextractr::exact_extract(
          flow_acc_raw, river_aoi,
          progress = FALSE
        ) %>% 
        rbindlist(idcol = "id") %>%
        .[, .(value = sum(coverage_fraction * value, na.rm = TRUE) / sum(coverage_fraction, na.rm = TRUE)), by = id]

      res_values <- flow_acc_values$value
    }
  }

  res_tbl <- 
      data.table(
        id = unique(aoi$id),
        continent = aoi$continent,
        any_rivers = !all(is.na(res_values)),
        fa_q05 = quantile(res_values, 0.05, na.rm = TRUE)[[1]],
        fa_q10 = quantile(res_values, 0.10, na.rm = TRUE)[[1]],
        fa_q25 = quantile(res_values, 0.25, na.rm = TRUE)[[1]],
        fa_q50 = quantile(res_values, 0.50, na.rm = TRUE)[[1]]
      )
  return(res_tbl)
}



#/*--------------------------------*/
#' ## Run for each continent
#/*--------------------------------*/
library(future.apply)
plan(multicore, workers = 10)
library(progressr)
handlers("txtprogressbar")
handlers(global = TRUE)

ls_region <- unique(case_hybas_tbl$region_abbr)

flow_acc_tbl <- 
  foreach(r = ls_region, .combine = rbind) %do% {
    # r = "au"
    message(paste0("Processing region: ", r))
    
    # --- Subset Watershed data for region r --- #
    case_hybas_r <- 
      filter(case_hybas_tbl, region_abbr == r) %>%
      .[["tg_hybas_sf"]] %>%
      .[[1]]
    
    # --- Load River data for region r --- #
    river_raw <- 
      st_read(
        file.path(base_dt_dir, "Raw/Rivers", paste0("HydroRIVERS_v10_", r, "_shp"), paste0("HydroRIVERS_v10_", r, ".shp")),
        quiet = TRUE
      ) %>%
      st_transform(crs = st_crs(case_hybas_r))
    
    # --- Load Flow Accumulation data for region r --- #
    flow_acc_raw <- 
      rast(file.path(base_dt_dir, "Raw/FlowAccumulation", paste0("hyd_", r, "_acc_15s.tif"))) %>%
      project(., "epsg:6933") #  WGS 84 / World
    

    # === Create Flow Accumulation Table for region r === #
    xs <- 1:nrow(case_hybas_r)
    
    with_progress({
      p <- progressor(along = xs)

      flow_acc_tbl <- 
        lapply(
        # future_lapply(
          xs,
          function(x) {
            message(paste0("Processing region: ", r, ", watershed ID: ", case_hybas_r$id[x]))
            # p()
            
            get_flow_acc_values(
              aoi = case_hybas_r[x, ],
              river_raw = river_raw,
              flow_acc_raw = flow_acc_raw
            )
          }
          # future.seed = NULL
        ) %>%
        rbindlist()
    })

    saveRDS(
      flow_acc_tbl,
      file.path(out_interm_dir, "FlowAcc", paste0("flow_acc_tbl_", r, ".rds"))
    )

    return(flow_acc_tbl)  
}

# === Save === #
saveRDS(
  flow_acc_tbl,
  here(out_dir, "FlowAcc/flow_acc_tbl.rds")
)

