# coding=UTF-8
# -----------------------------------------------
# Model: Seasonal Water Yield
# This is a sample python script generated from InVEST SWY model.
# `x_SWY_Runner.py`, which will be used to run SWY model for GWR project, is created based on this sample code.
# This script itself is not intended to be run directly, but just for reference.

import logging
import sys

import natcap.invest.seasonal_water_yield.seasonal_water_yield
import natcap.invest.utils

LOGGER = logging.getLogger(__name__)
root_logger = logging.getLogger()

handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    fmt=natcap.invest.utils.LOG_FMT,
    datefmt='%m/%d/%Y %H:%M:%S ')
handler.setFormatter(formatter)
logging.basicConfig(level=logging.INFO, handlers=[handler])


args = {
    'alpha_m': '',
    'aoi_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/aoi.gpkg',
    'beta_i': '0.28',
    'biophysical_table_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/biophysical_table.csv',
    'climate_zone_raster_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/climate_zone.tif',
    'climate_zone_table_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/rain_events_by_cz.csv',
    'dem_raster_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/dem.tif',
    'et0_raster_table': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/et0_path_tbl.csv',
    'flow_dir_algorithm': 'D8',
    'gamma': '0.84',
    'l_path': '',
    'lulc_raster_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/lulc.tif',
    'monthly_alpha': True,
    'monthly_alpha_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/m_alpha.csv',
    'precip_raster_table': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/Precipitation_path_tbl.csv',
    'rain_events_table_path': '',
    'results_suffix': '',
    'soil_group_path': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150/soil_group.tif',
    'threshold_flow_accumulation': '37',
    'user_defined_climate_zones': True,
    'user_defined_local_recharge': False,
    'workspace_dir': '/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Final/SWY/af/1060000150',
}

if __name__ == '__main__':
    natcap.invest.seasonal_water_yield.seasonal_water_yield.execute(args)
