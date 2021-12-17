##########################################################################################
## 
##  Function that reads BioSim outputs containing monthly climatic projections
##  for the study area of 7 variables: min temp, max temp, mean temp, accum precip
##  mean wind speed, max wind speed, and degree days 5ºC.
##  For each year and variable a raster brick is built and saved as a .rds file
##
##########################################################################################

rm(list=ls())
library(sf)
library(sp)
library(raster)
library(rasterVis)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select


biosim.maps = function(year=2021, export.file = "Export (DegreeDay5)", 
                       path="D:/OneDrive - ctfc.cat/QBCMOD/"){
  
  ## Locations
  loc = read.csv(paste0(path, "BioSIM/Quebec/Loc/Specific Locations Mask022d.csv"))
  
  ## Read all BioSIM outputs if no specified
  if(is.null(export.file)){
    dd5 = read.csv(paste0(path, "BioSIM/Quebec/Output/Export (DegreeDay5) ", year, ".csv"))
    names(dd5)[ncol(dd5)] = "Degree.Day.5"
    prec = read.csv(paste0(path, "BioSIM/Quebec/Output/Export (Analysis) sum ", year, ".csv"))
    names(prec)[ncol(prec)] = "Total.Precipitation"
    wind = read.csv(paste0(path, "BioSIM/Quebec/Output/Export (Analysis) max ", year, ".csv"))
    names(wind)[ncol(wind)] = "Wind.Speed.at.10.meters.Highest"
    tw = read.csv(paste0(path, "BioSIM/Quebec/Output/Export (Analysis) mean ", year, ".csv"))
    dta = tw %>%  left_join(wind, by=c("KeyID", "Year", "Month")) %>%
      left_join(prec, by=c("KeyID", "Year", "Month")) %>% left_join(dd5, by=c("KeyID", "Year", "Month"))
    rm(dd5); rm(prec); rm(wind); rm(tw); gc()
  }
  else{
    dta = read.csv(paste0(path, "BioSIM/Quebec/", file, " ", year, ".csv"))
  }
  
  ## Raster map of the study area, Quebec province 
  load(paste0(path, "SBW/data/mask.rda"))  
  # crs(mask)  # --> unfortunatelly, this projection is not supported by EPSG, 
  # neither the projections of the fire regime zones map
  
  ## For each variable and month... build maps?
  clim.vars = names(dta)[-(1:3)]
  for(var in clim.vars){
    brick_clim = NULL
    file.name = ifelse(var=="Minimum.Air.Temperature", "tmn", 
                       ifelse(var=="Air.Temperature", "temp", 
                              ifelse(var=="Maximum.Air.Temperature", "tmx", 
                                     ifelse(var=="Total.Precipitation", "prec",
                                            ifelse(var=="Wind.Speed.at.10.meters", "wind", 
                                                   ifelse(var=="Wind.Speed.at.10.meters.Highest", "windmx", "dday"))))))
    cat(paste0("Build raster brick for ", var), "\n")
    for(month in 1:12){
      cat(paste0("  ", month, "/", year), "\n")
      aux = dta %>% filter(Year==year, Month==month) %>% select(KeyID, all_of(var))
      df_clim = loc %>% select(KeyID, Latitude, Longitude) %>% left_join(aux, by="KeyID")
      # Instead of :
      # x = data.frame(x=clim.dta[,4])
      # spdf_clim = SpatialPointsDataFrame(coords = df_clim[,c("Longitude", "Latitude")], data = x,
      #                                    proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      # Update to sf object:
      spdf_clim = st_as_sf(df_clim, coords = c("Longitude", "Latitude"))
      st_crs(spdf_clim) = 4326
      
      # Transform to Lambert Conformal Conic and rasterize using mask as the raster of reference
      spdf_clim_lcc = st_transform(spdf_clim, crs = crs(mask))  
      raster_clim = rasterize(spdf_clim_lcc, mask, field=var)
      
      ## Build a RasterBrick for each variable, with one raster per month
      if(is.null(brick_clim)){
        brick_clim = brick(raster_clim)
      }
      else{
        brick_clim[[month]] = raster_clim
      }
      
      # write a monthly raster? no, it's better to save a single brick
      # writeRaster(raster_clim, filename=paste0(path, "DataOutSBW/biosim.maps/", 
      #             file.name, "_", year, ".", month, ".tiff"), format="GTiff", overwrite=TRUE)  
    }
    names(brick_clim) = paste0(file.name, "_", year, ".", 1:12)
    saveRDS(brick_clim, paste0(path, "DataOutSBW/biosim.maps/", year, "_", var, ".rds"))
    png(filename = paste0(path, "DataOutSBW/biosim.maps/", year, "_", file.name, ".png"), width=700, height=700)
    print(levelplot(brick_clim, main=paste(year, var), par.settings = BuRdTheme))
    dev.off()
  }
  
}


############ Run the biosim.map function for a single year ############
year=2016; export.file = NULL; path="D:/OneDrive - ctfc.cat/QBCMOD/"
for(year in 1997:2015)
  biosim.maps(year, export.file, path)


############ To restructurate some BioSIM outputs files ############
year=2020; export.file = NULL; path="D:/OneDrive - ctfc.cat/QBCMOD/"
prec = read.csv(paste0(path, "BioSIM/Quebec/Output/DegreeDay5.csv"))
aux = group_by(prec, KeyID, Year, Month) %>% summarise(x=sum(Wind.Speed.at.10.meters.Highest.))
names(aux)[4] = "Wind.Speed.at.10.meters.Highest"
write.csv(aux, paste0(path, "BioSIM/Quebec/Output/Export (Analysis) max 2021.csv"), quote=F, row.names=F)
#KeyID,Year,Month,Total Precipitation(Sum)


