rm(list=ls())
library(sp)
library(rgdal)
library(raster)
library(rgeos)

weather = "DailyWeather"
site = "Canada-USA"
period = "1980-2020" 
path = "C:/WORK/OneDrive - ctfc.cat/QBCMOD/"


buffer_metostat = function(weather="NormalWeather", site="Canada-USA", period="1981-2010", 
                       path="C:/WORK/OneDrive - ctfc.cat/QBCMOD/"){
  
  ## Meteorological stations
  dta = read.csv(paste0(path, "BioSIM/", weather, "/", site, "_", period, "/", 
                        site, " ", period, ifelse(weather=="NormalWeather", ".NormalsHdr.csv", ".DailyHdr.csv")))  
  
  ## Build points layers with elevation
  elevation = as.data.frame(dta[,"Elevation"])
  names(elevation) = "elev"
  sp_meteo = SpatialPointsDataFrame(coords = dta[,c("Longitude", "Latitude")], data = elevation,
                                    proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  plot(sp_meteo, col="blue", pch=19, cex=0.5)
  writeOGR(sp_meteo, dsn = paste0(path, "BioSim/MeteoStations"), layer=paste0(site,"_",period), driver="ESRI Shapefile", overwrite_layer=TRUE)
  


  ## Cartographic projection of layers to Lambert Conical Conforme
  zones = readOGR("C:/WORK/onedrive - ctfc.cat/QBCMOD/DataIn/ZonageFeux/2020.11.27/zones_nuria2.shp")
  raster::crs(zones)
  sp_meteo_lcc = spTransform(sp_meteo, raster::crs(zones))  
  plot(zones, col="grey90")
  plot(sp_meteo_lcc, add=T, col="blue", pch=19, cex=0.5)
  
  # Tryin to plot... it does not work
  # zones_latlon = spTransform(zones, raster::crs(sp_meteo))
  # plot(zones_latlon, col="grey90")
  # plot(sp_meteo, add=T, col="blue", pch=19, cex=0.5)

  dta = read.csv(paste0(path,"BioSIM/NormalWeather/Quebec++1981-2010(adjusted_from_2003-2017)/Quebec++ 1981-2010 (adjusted from 2003-2017).NormalsHdr.csv"))
  elevation = as.data.frame(dta[,"Elevation"])
  names(elevation) = "elev"
  sp_meteo_qbc = SpatialPointsDataFrame(coords = dta[,c("Longitude", "Latitude")], data = elevation,
                                    proj4string=CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  sp_meteo_qbc_lcc = spTransform(sp_meteo_qbc, raster::crs(zones))  
  plot(sp_meteo_qbc_lcc, add=T, col="violet", pch=19, cex=0.5) 
  
  
  ## Build circular buffer around meteo stations of 2.7º ~ 300 km
  buff = gBuffer(sp_meteo_qbc, width=2.7)
  plot(buff, col="forestgreen")
  plot(sp_meteo, add=T, col="blue", pch=19, cex=0.5) 
  plot(sp_meteo_qbc, add=T, col="violet", pch=19, cex=0.5) 
  buff_lcc = spTransform(buff, raster::crs(zones))  
  df_buff_lcc = SpatialPolygonsDataFrame(buff_lcc, data=as.data.frame(1), proj4string = raster::crs(zones))
  raster::shapefile(x=buff_lcc, file=paste0(path, "BioSim/MeteoStations/buff_meteo_qbc.shp"), overwrite=T)
  # writeOGR(df_buff_lcc, dsn = paste0(path, "BioSim/MeteoStations"), 
  #          layer="buff_meteo_qbc", driver="ESRI Shapefile", overwrite_layer=TRUE)  
  
}


