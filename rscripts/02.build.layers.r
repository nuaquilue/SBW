library(lattice)
library(sp)
library(raster)
library(tidyverse)

## HOST SUITABILITY MAP
load(file="inputlyrs/rdata/land.rdata")
breaks <- c(0,20,40,60,80,100,999)
tags <- c("C10","C30", "C50", "C70", "C90", "OLD")
land$age.class <- cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
land$sbw.suscep <- ifelse(land$spp %in% c("SAB", "OTH.RES.N") & land$age.class %in% c("C50", "C70", "C90", "OLD"), 1,
                   ifelse(land$spp %in% c("SAB", "OTH.RES.N") & land$age.class %in% c("C10", "C30"), 0.75,
                   ifelse(land$spp %in% c("EPN", "OTH.RES.S") & land$age.class %in% c("C50", "C70", "C90", "OLD"), 0.5,
                   ifelse(land$spp %in% c("EPN", "OTH.RES.S") & land$age.class %in% c("C10", "C30"), 0.25, 0))))
load(file="inputlyrs/rdata/mask.rdata")
MAP <- MASK
MAP[!is.na(MASK[])] <- land$sbw.suscep
plot(MAP, col=heat.colors(5)[5:1])


## ELEVATION
ELEV <- raster("C:/WORK/OneDrive - ctfc.cat/QBCMOD/DataIn/Elevation/elev_QC3.grd")
crs(ELEV)
ELEVlcc <- projectRaster(ELEV, res=2000, crs=crs(MASK))
ELEVlcc.clip <- crop(ELEVlcc, MASK)
load(file="inputlyrs/rdata/land.rdata")
land$elev <- ELEVlcc.clip[!is.na(MASK[])]
## Assign mean elevation of the neigbour cells to those cells with non informed values
zcells <- filter(land, !is.na(spp) & is.na(elev))
r <- 3
while(nrow(zcells)>0){
  neighs <- nn2(select(land, x, y), select(zcells, x,y), searchtype="priority", k=r^2)
  values <- matrix(land$elev[neighs$nn.idx], ncol=r^2)
  land$elev[!is.na(land$spp) & is.na(land$elev)] <- apply(values, 1, mean, na.rm=T)
  zcells <- filter(land, !is.na(spp) & is.na(elev))
  r <- r+2
}
save(land, file="inputlyrs/rdata/land.rdata")

## Latitude - Longitude - Elevation table for BioSIM, saved as .csv file 
ELEV <- raster("C:/WORK/OneDrive - ctfc.cat/QBCMOD/DataIn/Elevation/elev_QC3.grd")
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV)), Latitude=coordinates(ELEV)[,2],
                  Longitude=coordinates(ELEV)[,1], Elevation=ELEV[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation)) # 3.275.127 to 2.074.993
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, "C:/WORK/OneDrive - ctfc.cat/QBCMOD/BioSIM/Quebec/Loc/Specific Locations 011d.csv", 
          quote=F, row.names=F)
## Lower resolution by 2
ELEV2 = raster::aggregate(ELEV, fact=2, fun=mean)
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV2)), Latitude=coordinates(ELEV2)[,2],
                  Longitude=coordinates(ELEV2)[,1], Elevation=ELEV2[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation)) # 819.764 to 523.167
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, "C:/WORK/OneDrive - ctfc.cat/QBCMOD/BioSIM/Quebec/Loc/Specific Locations 022d.csv", 
          quote=F, row.names=F)

## Crop Elevation map by the study area 
zones <- readOGR("C:/WORK/onedrive - ctfc.cat/QBCMOD/DataIn/ZonageFeux/2020.11.27/zones_nuria2.shp")
zones_latlon = spTransform(zones,raster::crs(ELEV))
plot(zones_latlon)
ZONES = rasterize(zones_latlon, ELEV2)
data = data.frame(Name=paste0("cell_", 1:ncell(ELEV2)), Latitude=coordinates(ELEV2)[,2],
                  Longitude=coordinates(ELEV2)[,1], Elevation=ELEV2[], Mask=ZONES[], State="Quebec", Country="Canada")
data = data %>% filter(!is.na(Elevation) & !is.na(Mask)) # 819.764 to 259.746
data = data %>% add_column(KeyID=1:nrow(data), .before="Name")
write.csv(data, "C:/WORK/OneDrive - ctfc.cat/QBCMOD/BioSIM/Quebec/Loc/Specific Locations Mask022d.csv", 
          quote=F, row.names=F)


## VULNERABILITY MAP
source("mdl/sbw.vulnerability.r")
vuln <- vulnerability(land, MASK, 0.5)
MAP <- MASK
MAP[!is.na(MASK[])] <- vuln$v
# plot(MAP, col=c("grey80", heat.colors(5)[5:1]))
plot(MAP, col=heat.colors(5)[5:1])
group_by(vuln, spp, factor(host.pref)) %>% summarise(v=mean(v))
