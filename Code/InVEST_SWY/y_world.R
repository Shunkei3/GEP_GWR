library(data.table)
library(sf)
library(ggplot2)

library(countrycode)


countrycode(sourcevar = df[, "country"],
                            origin = "country.name",
                            destination = "continent")

                            
world_bd <- 
  st_read(file.path(base_dt_dir, "Raw/World_Data_Grids/country.shp")) %>%
  setnames(names(.), tolower(names(.)))
st_crs(world_bd) <- 4326
