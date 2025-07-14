# Directories for Raw Data
base_dt_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data"


# Directories for Intermediate Output
out_interm_dir <- "/Volumes/ext_crucial/NatCapTEEMs/GEP/GWR/Data/Intermediate"

if(!dir.exists(out_interm_dir)){
  dir.create(out_interm_dir, recursive = TRUE)
}

