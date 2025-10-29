+ You should redo the selection of watersheds:
  + Instead of 20% threshold, use the water withdrawal threshold.
+ Also, when you select watersheds, make sure to assign continent names correctly.
+ Move some code from "3_Local_Recharge.qmd" to "Prep_Watersheds.qmd" to make the workflow clearer.
+ Also, be careful about the geometry type in a data. It contains GEOMETRYCOLLECTION, which may cause errors in some operations.
  + table(st_geometry_type(ctry_bd_r))
  + This might be useful

```{r}
ctry_mask <- ctry_bd_r %>%
  st_make_valid() %>%                    
  st_zm(drop = TRUE, what = "ZM") %>%    
  st_collection_extract("POLYGON") %>% 
  st_cast("MULTIPOLYGON") %>%
  filter(!st_is_empty(.))
```


+ terra::disagg
  + https://gis.stackexchange.com/questions/456089/casting-all-geometries-from-spatvector-from-type-multipolygon-to-type-polygon