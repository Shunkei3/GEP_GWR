import sys
import logging
from pathlib import Path

import natcap.invest.utils
from natcap.invest.seasonal_water_yield import seasonal_water_yield

LOGGER = logging.getLogger(__name__)
handler = logging.StreamHandler(sys.stdout)
formatter = logging.Formatter(
    fmt=natcap.invest.utils.LOG_FMT,
    datefmt='%m/%d/%Y %H:%M:%S ')
handler.setFormatter(formatter)
logging.basicConfig(level=logging.INFO, handlers=[handler])


# /*===== Helper function to find required files =====*/
# Returns absolute path to <base>/<name>, or raises if missing
def _find_required_file(base: Path, name: str) -> str:
    p = (base / name).resolve()
    if not p.exists():
        raise FileNotFoundError(f"Required file missing: {p}")
    return str(p)


# /*===== SWY function =====*/
def run_swy(
    workspace_dir: str,
    inputs_dir: str,
    *,
    threshold_flow_accumulation: int = 50,
    beta_i: float = 1.0,
    gamma: float = 1.0,
    quiet: bool = True,
) -> None:
    
    # /*---- Paths ----*/
    inputs = Path(inputs_dir).resolve()
    workspace = Path(workspace_dir).resolve()
    workspace.mkdir(parents=True, exist_ok=True)
    
    # Collect inputs
    aoi   = _find_required_file(inputs, "aoi.gpkg")
    dem   = _find_required_file(inputs, "dem.tif")
    lulc  = _find_required_file(inputs, "lulc.tif")
    soil  = _find_required_file(inputs, "soil_group.tif")
    bpt   = _find_required_file(inputs, "biophysical_table.csv")
    precip_tbl = _find_required_file(inputs, "Precipitation_path_tbl.csv")
    et0_tbl    = _find_required_file(inputs, "et0_path_tbl.csv")
    cz_ras = _find_required_file(inputs, "climate_zone.tif")
    cz_tbl = _find_required_file(inputs, "rain_events_by_cz.csv")
    m_alpha = _find_required_file(inputs, "m_alpha.csv")
    
    # /*---- SWY args ----*/
    args = {
        "workspace_dir": str(workspace),
        "results_suffix": "",

        "lulc_raster_path": lulc,
        "biophysical_table_path": bpt,
        
        "dem_raster_path": dem,
        "aoi_path": aoi,
        
        "flow_dir_algorithm": "D8",
        "threshold_flow_accumulation": int(threshold_flow_accumulation),
        "beta_i": float(beta_i),
        "gamma": float(gamma),
        
        "user_defined_local_recharge": False,
        "l_path": "",
        "precip_raster_table": precip_tbl,
        "et0_raster_table": et0_tbl,
        "soil_group_path": soil,
        
        "monthly_alpha": True,
        "alpha_m": "",
        "monthly_alpha_path": m_alpha,
        
        "user_defined_climate_zones": True,
        "rain_events_table_path": "",
        "climate_zone_table_path": cz_tbl,
        "climate_zone_raster_path": cz_ras,
    }
    
    # /*---- Run ----*/
    if quiet:
        # Disable all logging temporarily
        prev_threshold = logging.root.manager.disable
        logging.disable(logging.CRITICAL)
        try:
            seasonal_water_yield.execute(args)
        finally:
            logging.disable(prev_threshold)  # restore
    else:
        seasonal_water_yield.execute(args)


# --- Test Run --- #
# base = Path("/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000150")
# print(_find_required_file(base, "aoi.gpkg"))

# workspace_dir = "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Final/SWY/af/1060000160"
# inputs_dir = "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/1060000160"
# threshold_flow_accumulation = 37
# beta_i = 0.28
# gamma = 0.84

# --- Example usage ---
# tmp_id = 1060000160
# run_swy(
#     workspace_dir = "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Final/SWY/af/" + str(tmp_id),
#     inputs_dir="/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs/af/" + str(tmp_id),
#     threshold_flow_accumulation=37,
#     beta_i=0.28,
#     gamma=0.84,
# )
