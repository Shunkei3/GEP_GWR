# Data (base year: 2019)


## Primary Data for InVEST Seasonal Water Yield (SWY) Model

+ Land Use/Land Cover data:
  + ESA, Land Cover CCI
  + 0.002777778 degree resolution (?, approximately 300 m)
  + European Space Agency Climate Change Initiative land cover at 300m resolution?

+ Digital Elevation Model (DEM):
  + [Voild-filled DEM from HydroSHEDS v1](https://www.hydrosheds.org/hydrosheds-core-downloads)
    + 15 arc-seconds resolution (approximately 500 m)
    + based on The Shuttle Radar Topography Mission (SRTM) elevation data

+ **Watershed boundary data**:
  + [HydroBASINS (level 6) from HydroSHEDS v1](https://www.hydrosheds.org/products/hydrobasins)
    + Lehner, B., Grill G. (2013). Global river hydrography and network routing: baseline data and new approaches to study the world’s large river systems. Hydrological Processes, 27(15): 2171–2186. https://doi.org/10.1002/hyp.9740


+ **Precipitation and Evapotranspiration data**:
  + [AgERA5 historic and near real time forcing data (version 2.0) from Copernicus Climate Change Service (C3S) Climate Data Store (CDS)](https://cds.climate.copernicus.eu/datasets/sis-agrometeorological-indicators?tab=overview)
  + 0.1 degree resolution
    + Copernicus Climate Change Service (2020): Agrometeorological indicators from 1979 to present derived from reanalysis. Copernicus Climate Change Service (C3S) Climate Data Store (CDS). DOI: 10.24381/cds.6c68c9bb (Accessed on DD-MMM-YYYY)


+ **Climate zone data**
  + [Köppen-Geiger classification map for the period 1991-2020 downloaded from GloH2O](https://www.gloh2o.org/)
    + 0.1 degree resolution 
    + Beck et al. (2023) : High-resolution (1 km) Köppen-Geiger maps for 1901–2099 based on constrained CMIP6 projections, Scientific Data 10, 724 (2023).

+ **Soil Hydrologic Groups**
  + [HYSOGs250m, v1](https://daac.ornl.gov/SOILS/guides/Global_Hydrologic_Soil_Group.html)
    + 250 m resolution
    + @Ross.etal2018, Ross, C.W., L. Prihodko, J.Y. Anchang, S.S. Kumar, W. Ji, and N.P. Hanan. 2018. Global Hydrologic Soil Groups (HYSOGs250m) for Curve Number-Based Runoff Modeling. ORNL DAAC, Oak Ridge, Tennessee, USA. https://doi.org/10.3334/ORNLDAAC/1566

+ **Biophysical Tables**:
  + **Curve Numbers**
    + [GCN250, global curve number datasets for hydrologic modeling and design](https://figshare.com/articles/dataset/GCN250_global_curve_number_datasets_for_hydrologic_modeling_and_design/7756202)
    + @Jaafar.Ahmad2019
    + Jaafar, Hadi; Ahmad, Farah (2019). GCN250, global curve number datasets for hydrologic modeling and design. figshare. Dataset. https://doi.org/10.6084/m9.figshare.7756202.v1
  + **Monhtly Crop Evapotranspiration Coefficient**
    + FAO Irrigation and Drainage Paper 56 (Chapter 6, 1998)
      + @Allen.etal1998
    + @Nistor.etal2018


## Other Data

+ **Major GW bains boundary data**:
  + Global Geo-processed Data of Aquifer Properties by 0.5° Grid, Country and Water Basins (https://data.msdlive.org/records/9sany-rht26)
  + This dataset is used to identify the watersheds that overlap with major groundwater basins.

+ **A global layer of streams**:
  + [HydroRIVERS v1 from HydroSHEDS](https://www.hydrosheds.org/products/hydrorivers)
    + Lehner, B., Grill G. (2013). Global river hydrography and network routing: baseline data and new approaches to study the world’s large river systems. Hydrological Processes, 27(15): 2171–2186. https://doi.org/10.1002/hyp.9740
  + This is used to calculate the flow accumulation threshold for each watershed in the InVEST SWY modeling process.
  + "Global gidded CN dataset (GCN250) is generated from the ESA CCI-LC maps and the HYSOGs250m soil data based on the USDA curve number tables and plant functional types"

+ **Reflectance Based Kc (NDVI-Kc methods)**
+ Terra Vegetation Indices Monthly L3 Global 0.05Deg CMG V061
+ For example:
  + @Bausch1995 for corn
  + @Choudhury.etal1994
  + See @Glenn.etal2010
    + "Bausch and Neale (1987, 1989) showed that VIs could increase the accuracy of crop coefficients by providing a measurement of the actual state of the crop canopy during development."
    + "Choudhury et al. (1994) explored the theoretical basis of using VIs to replace Kc in Eq. 5."
    + However, on monthly time steps, the NDVI-Kc methods developed by Hunsaker et al. (2005a, b) and others (Neale et al. 2005)
  + @Kamble.etal2013
    + K_{c,ndvi} = 1.4751 * NDVI - 0.1725 for rainfed and irrigated agriculture 
  + @Goffin.etal2022
    + K_{c,NDVI} = NDVI ×2.70
  + @Reyes-Gonzalez.etal2015
    + Kc FAO-56 = 1.62 NDVI – 0.1471 (corn) (3)
    + Kc FAO-56 = 2.11 NDVI – 0.4989 (alfalfa) (4)

| Study                      | Crop/Region         | Equation                            |
| -------------------------- | ------------------- | ----------------------------------- |
| **Bausch (1995)**          | Corn (Colorado, US) | $K_c = 1.37 \cdot NDVI - 0.06$      |
| **Campos et al. (2010s)**  | Grapes, peaches     | $K_c = 1.41 \cdot NDVI - 0.17$      |
| **Duchemin et al. (2006)** | Wheat (Morocco)     | $K_c = 1.33 \cdot NDVI - 0.08$      |
| **Glenn et al. (2011)**    | Global summary      | $K_c \approx 1.2 \cdot NDVI + 0.05$ |



+ **Crop coefficient (Kc) estimation from leaf area index (LAI)** 
  + See https://community.naturalcapitalproject.org/t/crop-coefficient-kc-estimation-from-leaf-area-index-lai/2544/3



# Replicate
+ Global GEP Total in base year
+ A csv of GEP by country in base year
  + Format?
+ A PNG map of that map

+ Create a `run.py` runs the calculation


+ Github organization: NatCapTEEMs
+ "GEP_Dev"
+ 
+ "Project" folder
  + Global

+ EE-spec country names (Ee_r264_id, r264_correspondence) 
  + "M49" id: This is the actual ISO3 code

+ GEP 1000USD


**Methods: Currency conversion**
+ See method used by ESDV: ESDV_Global_Update-FINAL-Report-June-2020.pdf
+ Adjust county-specific currencies to common base year
+ Convert local currencies to 2019 international dollars using OWrld Bank PPP-adjusted exchange rates 