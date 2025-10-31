# Calculating GEP of Groundwater Recharge Service

This repository contains the code used to calculate the Global Ecosystem Product (GEP) of Groundwater Recharge Service. The majority of the code is written in R, with a small portion in Python for automating model execution.

## Path Configuration
`0_path_config.R` sets up the directory paths to the data folders used in this analysis. The variable `p_dir` in this script must be updated to point to the base (project) directory where the `Data` folder and `Scr` folder are located on your local machine. All the raw input data (available from the shared Google Drive folder) should be placed in the `Data/Raw` directory. This script also creates required intermediate and final output directories if they do not already exist.

## Data Preparation for InVEST Seasonal Water Yield Model
`x_PrepData` directory contains R scripts to preprocess raw input data for the InVEST Seasonal Water Yield (SWY) model. These data preparation steps include 

+ Identifying watersheds of interest, 
+ Aggregating global daily weather data to monthly data
+ Loading NDVI data and converting it to decimal form
+ Reclassifying soil datasets to generate a soil hydrologic group map

The prepared input data are then used to generate watershed-specific input folders for the SWY model.

The script `1_Gen_SWY_Inputs.qmd` creates the watershed-specific input data required to run the SWY model. It generates one input folder per watershed, stored in:

```
Data/Intermediate/SWY_inputs/<region>/<watershed_id>
```

Each watershed directory is automatically generated with all the necessary input data files. `x_SWY_Functions.qmd` contains all the functions used in the `1_Gen_SWY_Inputs.qmd` script.


## Running the SWY Model

`x_SWY_Runner.py` is written in Python and it serves as a wrapper script to automate the execution of the InVEST Seasonal Water Yield (SWY) model. For consistency of the computing environment, the Python functions written in this script are loaded in R using the `reticulate` package. Then, the SWY model is run for multiple watersheds in parallel.


## GEP Calculation and Result Preparation

`3_Local_Recharge.qmd` derive local groundwater recharge raster data from the SWY model outputs. The values in the raster data are aggregated to country level by calculating the area-weighted average of groundwater recharge rates within each country boundary in `4_Country_L_values.qmd`. 

Finally, `5_GEP_calculation.qmd` calculates the county-level GEP of groundwater recharge service by multiplying the country-level groundwater recharge rates and the country-level groundwater withdrawal amounts. The final GEP results are stored in `Data/Final/ground_water_gep.csv`.



## Code for Writing the Manuscript

The `Writing` directory contains Quarto files used to generate the manuscript associated with this analysis. Figures and Tables for the manuscript are produced in `x_PrepareResults.qmd`.


















