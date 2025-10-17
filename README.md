# Estimate Groundwater Recharge using InVEST SWY Model

The InVEST SWY model is used to estimate the ecosystem production of groundwater recharge (GEP) across various regions. This document outlines the data requirements, mecalibrating the model parameters to replicate the analysis.


## Data Requirements

Below, I outline the data required for running the InVEST SWY model in the bold section. Also, I provide the source of each dataset and briefly describe how the input data is prepared for the InVEST SWY model. For more details on the data preparation, see the section of calibrating the model parameters.

+ **Watershed Boundary (AOI) Data**
  + I used major groundwater basin boundary data is obtained from the [@Niazi.etal2024] and the watershed boundary data from HydroBASINS v1 (Level 6) [@Lehner.Grill2013]. Specifically, I identified the watersheds that overlap with major groundwater basins and run the InVEST SWY model for each of those watersheds. In total, 7,722 watersheds are selected for the analysis. For each watershed, the input data is prepared and run the InVEST SWY model.

+ **Land Use/Land Cover Data**
  + LULC data is obtained from European Space Agency Climate Change Initiative land cover at 300m resolution

+ **Biophysical Table Data**
  + Curve number values from GCN250 dataset [@Jaafar.etal2019]
  + Global data of monthly crop evapotranspiration coefficient (Kc) is not available. I calculate it using the reflectance-based Kc (NDVI-Kc methods). The NDVI data is obtained from monthly data on MODIS/Terra Vegetation Indices with 0.05 degree resolution [@Didan2021]. See the section of calibrating the model parameters for details on the NDVI-Kc methods.

+ **Digital Elevation Model (DEM)**
  + Void-filled DEM data with 15 arc-seconds resolution is obtained from HydroSHEDS v1.1 [@Lehner.etal2008]. 

+ **Monthly precipitation and reference evapotranspiration data** 
  + Daily precipitation and reference evapotranspiration data is obtained AgERA5 dataset (version 2.0) from Copernicus Climate Change Service (C3S) Climate Data Store (CDS) [@CopernicusClimateChangeServic2020]. This is gridded weather datasets with 0.1 degree resolution. For each grid cell, I calculate the monthly mean precipitation and reference evapotranspiration values.
  
+ **Monthly Rain Events Table by Climate Zone**
  + A Rain event table is a table that contains the number of monthly rain events ($\ge 0.1 mm$) for each watershed. By default, the InVEST SWY model only requires rain event table for entire watershed. However, given the spatial variation of the rainfall events, I prepare the rain event table for each climate zone within each watershed.  Climate zone data is obtained from the Köppen-Geiger classification map for the period 1991-2020 downloaded from GloH2O [@Beck.etal2023]. For each watershed, I calculate the number of monthly rain events ($\ge 0.1 mm$) for each climate zone using the daily precipitation data.


+ **Monthly Alpha Table**
  + Following the description in the InVEST SWY model documentation, the monthly alpha table is prepared by calculating the ratio of the previous month's precipitation to the total precipitation for each watershed (i.e., $\alpha_{i,m} = P_{i,m-1} / \sum_{m=1}^{12} P_{i,m}$). 

+ **Soil Hydrologic Group**
  + HYSOGs250m dataset [@Ross.etal2018] is used to define the soil hydrologic groups (HSGs) for each watershed. This dataset provides the global hydrologic soil groups at 250m resolution.

+ **Threshold Flow Accumulation value**
  + I derived the threshold flow accumulation value for each watershed using the global layer of streams (HydroRIVERS v1) and a flow accumulation map obtained from HydroSHEDS dataset [@Lehner.Grill2013]. Specifically, I extracted the flow accumulation value along the streams in the watershed and used the 25th percentile of the flow accumulation values as the threshold value. 25th percentile is used to avoid the outliers in the flow accumulation values. 

+ **Beta_i Parameter**
  + For each watershed, Beta_i Parameter is calculated from the DEM data and the flow accumulation map of that watershed. See the section of calibrating the model parameters for details on how to calculate the Beta_i parameter.

+ **Gamma Parameter**
  + Using the soil hydrologic group data, I set the Gamma parameter for each watershed. See the section of calibrating the model parameters for details on how to calculate the Gamma parameter for each watershed.


## Calibrating the Model Parameters

In this section, I describe how I calibrate the model parameters of the SWY model to the local conditions of each watershed. The InVEST SWY model requires several parameters to be calibrated for each watershed, including the threshold flow accumulation value, the biophysical table data, Beta_i and Gamma parameters.  Although it is common to use the default values for these parameters, I calibrate them to better represent the local conditions of each watershed.

### Threshold Flow Accumulation Value

Using the flow accumulation map and the global layer of streams, I computed the threshold flow accumulation value for each watershed. The threshold flow accumulation value is the minimum number of upstream cells that contribute to a stream cell. I extracted the flow accumulation value along the streams in the watershed and used the 25th percentile of the flow accumulation values as the threshold value. 25th percentile is used to avoid the outliers in the flow accumulation values. 

### Biophysical Table Preparation
Biophysical table data for the InVEST SWY model requires soil-group specific curve numbers and monthly crop evapotranspiration coefficient (Kc) values for each land use/land cover (LULC) class. There is no single numbers for the curve number and Kc values for global scale. To account for the heterogeneity of the curve number and Kc values across regions, I calibrate the biophysical table data for each watershed. 

For watershed-specific curve numbers, I use the GCN250 dataset [@Jaafar.etal2019]. GCN250 dataset provides the gridded data of curve numbers for global scale at 250m resolution. Using the GCN250 dataset, the map of soil hydrologic groups (HYSOGs250m), and the LULC data for the watershed, I calculated the mean values of curve numbers for each soil group and LULC class.

Global Kc values are not available nor there is no consensus on the Kc values for each of the various land use/land cover classes. However, it is safe to assume that the Kc values for some land use/land cover classes are constant and similar across months and regions. For example, the Kc values for evergreen forests, grasslands, and wetlands are generally constant throughout the year. 


Therefore, I use constant Kc value for land use/land cover classes including including evergreen forest, deciduous forest, shrublands, urban, bare ground, snow/ice, and open water. These constant Kc values are borrowed from the previous studies [@Allen.etal1998; @Liu.etal2017]. Meanwhile, for the other land use/land cover classes such as rainfed or irrigated cropland and deciduous forests, the Kc values vary by month. For these classes, I calculate the monthly Kc values using the reflectance-based Kc (NDVI-Kc methods). This reflectance-based Kc method is widely used in the previous studies to estimate the Kc values for various land use/land cover classes (e.g., @Bausch1995; @Choudhury.etal1994; @Glenn.etal2010; @Kamble.etal2023; @Rayes-Gonzalez.etal2015). Following the previous studies, I calculate the monthly Kc values for each LULC class using the NDVI-Kc methods. The monthly NDVI data is obtained from the MODIS/Terra Vegetation Indices with 0.05 degree resolution [@Didan2021]. 

With NDVI data, Kc value is calculated using the following equation [@Bausch1995]:

$$K_c = 1.37 \cdot NDVI = 0.069$$


### Beta_i and Gamma Parameters

Another parameter that needs to be calibrated for each watershed is the Beta_i and Gamma parameters. The Beta_i parameter in the SWY model represents the fraction of upgradient water subsidy that is available for evapotranspiration at a given location. Following the description in the InVEST SWY model documentation, Beta_i parameter is calculated as the topographic wetness index (TWI):  

$$TWI = ln(\frac{A}{tan b})$$

, where $A$ is upstream contributing area (i.e. A = flow accumulation $\times$ cell area, in $m^2$) and $b$ is the slope in radians. $A$ is calculated from the flow accumulation map. The slope is calculated from the DEM data using the `terrain` function in the `terra` R package. After deriving TWI for each grid cell of the DEM data, I calculate the mean TWI value within each watershed. This mean TWI value is used as the watershed-specific Beta_i parameter for each watershed in the SWY model to represent the average topographic subsidy available to each watershed.

The Gamma parameter controls how much of the water recharge is made available to the downgradient pixels. Gamma = 0 means that no water moves downgradient from the pixel, while Gamma = 1 means that all local recharge becomes available to downgradient pixels. This parameter is set to 1 by default in the SWY model. However, I set the Gamma parameter for each watershed based on the soil hydrologic group raster data, which is one of the input data for the SWY model as described above. This raster dataset classifies each grid cell into one of four USDA-defined hydrologic soil groups: A (group 1), B (2), C (3), and D (4). As the group number increases, soils become more impermeable, with Group A having the highest infiltration rates and Group D having the lowest.  

Given this gradient in lateral subsurface connectively, I assigned following Gamma values to each group:  0.95 for Group A, 0.85 for Group B, 0.75 for Group C, and 0.55 for Group D. Then, for each watershed, I calculated an area-weighted average of these Gamma values based on the fractional coverage of each soil group. The resulting watershed-specific Gamma values are used as input to the SWY model to better represent spatial differences in subsurface flow connectivity.

















