# /*===========================================*/
#'=  Objective =
#' + Identify the watersheds that are overlapping with the aquifer boundaries.
# /*===========================================*/

library(data.table)
library(dplyr)
library(terra)
library(sf)
library(tmap)
library(ggplot2)
library(here)

base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"
out_interm_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

world_bd <- 
 st_read(file.path(base_dt_dir, "Raw/World_Data_Grids/country.shp")) %>%
 setnames(names(.), tolower(names(.))) 
st_crs(world_bd) = 4326

# === Major GW Basins === #
gw_basins_sf <- 
    st_read(file.path(base_dt_dir, "Raw/h2wjr-te507/shapefiles/All_merged.shp")) %>%
    filter(WHYClass == 10)

# st_crs(gw_basins_sf)
# ggplot() +
#     geom_sf(data = st_transform(world_bd, st_crs(gw_basins_sf))) +
#     geom_sf(data = gw_basins_sf)

# /*===========================================*/
#'= Identifying Target Watersheds =
# /*===========================================*/
fn_identify_wsheds <- function(region_abbr) {
    # Test run
    # region_abbr = "na" 

    # === Load watersheds data === # 
    dir <-  "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Raw/HydroBASINS"
    folder_name <- paste0("hybas_", region_abbr, "_lev01-12_v1c")
    file_name <- paste0("hybas_", region_abbr, "_lev06_v1c.shp")

    hybas_raw <- 
        st_read(file.path(dir, folder_name, file_name))
    
    # ggplot() +
    #     geom_sf(data = world_bd) +
    #     geom_sf(data = hybas_raw, fill = "blue")

    tg_hybas <-
        st_transform(hybas_raw, crs = st_crs(gw_basins_sf)) %>%
        .[gw_basins_sf,] %>%
        st_transform(crs = 6933) %>% # WGS 84 / World Cylindrical Equal Area
        mutate(
            ws_id = HYBAS_ID,
            continent = region_abbr,
            id = paste0(region_abbr, "_", HYBAS_ID)
        )

    # ggplot() +
    #     geom_sf(data = filter(world_bd)) + 
    #     geom_sf(data = tg_hybas) +
    #     geom_sf(data = gw_basins_sf, fill = "blue")

    return(tg_hybas)
}

case_hybas_tbl <- 
     data.table(
        region = c("Africa", "Asia", "Australia", "Europe", "North America", "South America"),
        region_abbr = c("af", "as", "au", "eu", "na", "sa")
    ) %>%
    tibble() %>%
    mutate(
        tg_hybas_sf = lapply(
            1:nrow(.),
            \(x) fn_identify_wsheds(region_abbr = .$region_abbr[x])
        )
    )

# case_hybas_tbl %>%
#     filter(region_abbr == "na") %>%
#     .[["tg_hybas_sf"]] %>%
#     .[[1]]

# === Save === #
saveRDS(
    case_hybas_tbl, 
    file.path(out_interm_dir, "case_hybas_tbl.rds")
)




