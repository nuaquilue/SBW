##########################################################################################
## 
##  Creation of location files for BioSim
##
##########################################################################################

rm(list=ls())
library(sp)
library(raster)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select

## Latitude - Longitude - Elevation table for BioSIM, saved as .csv file 
ELEV <- raster(paste0(dirname(getwd()), "/DataIn/Elevation/elev_QC3.grd"))
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV)), Latitude=coordinates(ELEV)[,2],
                  Longitude=coordinates(ELEV)[,1], Elevation=ELEV[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation)) # 3.275.127 to 2.074.993
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, paste0(dirname(getwd()), "/BioSIM/Quebec/Loc/Specific Locations 011d.csv"), 
          quote=F, row.names=F)

## Lower resolution by 2
ELEV2 = raster::aggregate(ELEV, fact=2, fun=mean)
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV2)), Latitude=coordinates(ELEV2)[,2],
                  Longitude=coordinates(ELEV2)[,1], Elevation=ELEV2[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation)) # 819.764 to 523.167
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, paste0(dirname(getwd()), "/BioSIM/Quebec/Loc/Specific Locations 022d.csv"), 
          quote=F, row.names=F)

## Crop Elevation map by the study area 
zones <- readOGR(paste0(dirname(getwd()), "/DataIn/ZonageFeux/2020.11.27/zones_nuria2.shp"))
zones_latlon = spTransform(zones,raster::crs(ELEV))
plot(zones_latlon)
ZONES = rasterize(zones_latlon, ELEV2)
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV2)), Latitude=coordinates(ELEV2)[,2],
                  Longitude=coordinates(ELEV2)[,1], Elevation=ELEV2[], Mask=ZONES[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation) & !is.na(Mask)) # 819.764 to 259.746
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, paste0(dirname(getwd()), "/BioSIM/Quebec/Loc/Specific Locations Mask022d.csv"), 
          quote=F, row.names=F)

