Sys.unsetenv("RETICULATE_PYTHON")

Sys.getenv("RETICULATE_PYTHON")
library(reticulate); py_config()

library(here)
readRenviron(here::here(".Renviron"))

library(reticulate)
py_config()
use_python("/Users/shunkeikakimoto/miniforge3/envs/gep_gwr/bin/python")





library(here)






use_python("/Users/shunkeikakimoto/miniforge3/envs/gep_gwr/bin/python", required = TRUE)

use_python("/Users/shunkeikakimoto/miniforge3/envs/gep_gwr/bin/python")
use_condaenv("gep_gwr")

library(reticulate)
library(here)

# Use your conda env (must actually contain Python):
use_condaenv(
  "/Users/shunkeikakimoto/miniforge3/envs/gep_gwr", 
  required = TRUE
)

# Inspect what R is using:
py_config()


file.exists(here("Scr/Code/InVEST_SWY/1_Run/test.py"))
source_python(here("Scr/Code/InVEST_SWY/1_Run/test.py"))


source_python(here("Scr/Code/InVEST_SWY/x_SWY_Runner.py"))

tmp_id = 1060000790
tmp_region = "af"
final_dir_swy_out <- "/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Final/SWY"
where_to_save <- file.path(final_dir_swy_out, tmp_region, tmp_id)
which_inputs_dir <- file.path("/Volumes/baseHD/NatCapTEEMs/GEP/GWR/Data/Intermediate/SWY_inputs", tmp_region, tmp_id)

dir.exists(which_inputs_dir)

run_swy(
    workspace_dir = where_to_save,
    inputs_dir = which_inputs_dir,
    threshold_flow_accumulation = 37,
    beta_i = 0.28,
    gamma = 0.84
)

