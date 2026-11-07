# Calculating GEP of Groundwater Recharge Service

This repository contains the code used to calculate the Gross Ecosystem Product (GEP) of groundwater recharge services, expressed in 2019 U.S. dollars. The majority of the code is written in R, with a small portion in Python for automating model execution. The analysis uses the InVEST Seasonal Water Yield (SWY) model, version 3.17.0.

## Path Configuration

`0_path_config.R` sets up the directory paths to the data folders used in this analysis. Update the variable `p_dir` to point to the directory containing your `Data` folder. All raw input data (available from the shared Google Drive folder) should be placed in `Data/Raw`. This script also creates the required intermediate and final output directories if they do not already exist.

The scripts use `Src/Code/` paths relative to the local project root. To match this layout, place this repository in a folder named `Src` within your local project directory; the repository's `Code` folder will then be located at `Src/Code`. If you use a different layout, update the code paths accordingly. `Src` stands for "source."

## Data Preparation for InVEST Seasonal Water Yield Model

The `x_PrepData` directory contains R scripts to preprocess raw input data for the InVEST Seasonal Water Yield (SWY) model. These data preparation steps include:

+ Identifying watersheds of interest
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

`x_SWY_Runner.py` is a Python wrapper script that automates execution of the InVEST Seasonal Water Yield model (version 3.17.0). Its Python functions are loaded in R using the `reticulate` package. The SWY model is then run for multiple watersheds in parallel. These executions are managed in `2_Run_SWY.qmd`.

## GEP Calculation and Result Preparation

`3_Local_Recharge.qmd` derives local groundwater recharge raster data from the SWY model outputs. `4_Country_L_values.qmd` aggregates these values to the country level by calculating the area-weighted average of groundwater recharge rates within each country boundary.

Finally, `5_GEP_calculation.qmd` calculates the country-level GEP of groundwater recharge services as the avoided pumping costs attributable to recharge. The calculation combines country-level groundwater recharge rates, groundwater withdrawal amounts, electricity prices, and an assumed pumping efficiency. Groundwater withdrawal data and electricity prices are selected for 2019. Electricity prices are converted from U.S. cents per kilowatt-hour to U.S. dollars per kilowatt-hour, and GEP values are expressed in 2019 U.S. dollars.

The final GEP results are stored in `Data/Final/ground_water_gep.csv`. The output includes country identifiers, the GEP value (`ground_water_gep`), and a `note` column explaining missing values. The note identifies missing recharge, groundwater withdrawal, or electricity price data, lists multiple reasons separated by semicolons, and is blank when all required inputs are available.
