# ---
# title: "Setting up the environment"
# ---

# /*===========================================*/
#'=  Path Setting =
# /*===========================================*/
# Reminder: I needed to do this because I was storing the data in an external drive, which is outside of the project directory. If you are storing the data withing the project directory, you may not need this part. Simply use `here::here()` function to set up the paths.


# Parent Directory (Change this part)
p_dir <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR"

# Data Directory
base_dir <- file.path(p_dir, "Data")

# /*---- Path to the Raw (input) data directory ----*/
raw_dir <- file.path(base_dir, "Raw")

# /*---- Create Intermediate output directories if they do not exist ----*/
int_dir <- file.path(base_dir, "Intermediate")

if (!dir.exists(int_dir)){
  dir.create(int_dir, recursive = TRUE)
}

# /*---- Create Final output directories if they do not exist ----*/
final_dir <- file.path(base_dir, "Final")

if (!dir.exists(final_dir)){
  dir.create(final_dir, recursive = TRUE)
}

# /*===========================================*/
#'=  Functions to create directory =
# /*===========================================*/
# I thought it would be useful to have a function to create directory if it does not exist.

fn_gen_save_dir <- function(where, new_relative_dir_path, return_path = FALSE){
  
  new_dir_path <- file.path(where, new_relative_dir_path)

  if(!dir.exists(new_dir_path)){
    dir.create(new_dir_path, recursive = TRUE)
  }

  if(return_path){
    return(new_dir_path)
  }
}

# Test
# fn_gen_save_dir(where = out_dt_p_dir, new_relative_dir_path = file.path(tmp_region, tmp_aoi_i$HYBAS_ID))