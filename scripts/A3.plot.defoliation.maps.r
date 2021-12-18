##########################################################################################
## 
##  Plot maps of 
##  * A. Host suitability
##  * B. Accumulated level of defoliation intensity
##  * C. Number of consecutive  years of defoliaion
##  * D. Current level of defoliation intensity
##  * E. SBW climatic niche
##  * F. SBW climatic niche  under climatic change
##  * G. Annual level of defoliation
##
##########################################################################################

rm(list=ls())
library(sp)
library(raster)
library(viridisLite)
library(viridis)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select

#### A. Host suitability 
load(file="data/land.sbw.rda")
load(file="data/mask.rda")
breaks = c(0,20,40,60,80,100,999)
tags = c("C10","C30", "C50", "C70", "C90", "OLD")
landscape$age.class = cut(landscape$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
landscape$sbw.suscep = ifelse(landscape$spp %in% c("SAB", "OTH.RES.N") & landscape$age.class %in% c("C50", "C70", "C90", "OLD"), 1,
                          ifelse(landscape$spp %in% c("SAB", "OTH.RES.N") & landscape$age.class %in% c("C10", "C30"), 0.75,
                                 ifelse(landscape$spp %in% c("EPN", "OTH.RES.S") & landscape$age.class %in% c("C50", "C70", "C90", "OLD"), 0.5,
                                        ifelse(landscape$spp %in% c("EPN", "OTH.RES.S") & landscape$age.class %in% c("C10", "C30"), 0.25, 0))))
map = mask
map[!is.na(mask[])] = landscape$sbw.suscep
plot(map, col=heat.colors(5)[5:1], main="Host suitability")


#### B. Accumulated level of defoliation intensity
map = mask
map[!is.na(mask[])] = landscape$cum.intens.def
mx.cum = max(landscape$cum.intens.def)
plot(map, col=c(rep("grey90",15), heat.colors(mx.cum-15+1)[(mx.cum-15+1):1]), 
     main="Accumulated level of defoliation intensity", xaxt="n", yaxt="n")
plot(map, col=c("grey90", heat.colors(mx.cum)[mx.cum:1]), xaxt="n", yaxt="n")


#### C. Number of consecutive  years of defoliaion
map[!is.na(mask[])] = landscape$ny.def
mx.ny = max(landscape$ny.def)
plot(map, col=terrain.colors(mx.ny)[mx.ny:1], 
     main="Number of consecutive years of defoliation", xaxt="n", yaxt="n")
plot(map, col=c("grey90", plasma(mx.ny)), xaxt="n", yaxt="n")


#### D. Current level of defoliation intensity
map[!is.na(mask[])] = landscape$curr.intens.def
plot(map, col=c("grey90", "yellowgreen", "orange", "darkred"), 
     main="Current level of defoliation intensity", xaxt="n", yaxt="n")


#### E. SBW climatic niche
landscape$sbw.niche = ifelse(landscape$temp>0.5 & landscape$temp<2.8, 1,
                         ifelse(landscape$temp>-1.5 & landscape$temp<4, 2, 3))
map[!is.na(mask[])] = landscape$sbw.niche
tiff(paste0(dirname(getwd()), "/DataOutSBW/sbw.niche/sbw.niche_2020.tiff"))
plot(map, col=c("yellowgreen", "orange", "darkred"), xaxt="n", yaxt="n")
dev.off()


#### F. SBW climatic niche  under climatic change
clim.scn = "rcp45"
load(file=paste0("data/temp_", clim.scn, "_MIROC_ESM_CHEM.rdata")) 
year = c(rep(NA,3), seq(2020,2095,5))
for(i in 4:19){
  cat(year[i], "\n")
  aux = cc.temp[,c(1,i)] 
  names(aux) = c("cell.id", "temp")
  land = select(land, -temp) %>% left_join(aux, by="cell.id")
  landscape$sbw.niche = ifelse(landscape$temp>0.5 & landscape$temp<2.8, 1, ifelse(landscape$temp>-1.5 & landscape$temp<4, 2, 3))
  map[!is.na(mask[])] = landscape$sbw.niche
  tiff(paste0(dirname(getwd()), "/DataOutSBW/sbw.niche/sbw.niche_",  clim.scn,"_",  year[i],  ".tiff"))
  plot(map, col=c("yellowgreen", "orange", "darkred"), xaxt="n", yaxt="n", main=paste(clim.scn, "-", year[i]))
  dev.off()
}


#### G. Annual level of defoliation
load(paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.intens.y.mask.rdata"))
for(y in 1968:2020){
  cat(y, "\n")
  defol.y = filter(sbw.intens.y, year==y)
  # get the information about defoliation for the target year       
  land.sbw = left_join(landscape, defol.y, by="cell.id")  
  land.sbw$intensity[is.na(land.sbw$intensity)] = 0
  # build the map 
  map = mask
  map[!is.na(mask[])] = land.sbw$intensity
  tiff(paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/sbw.defol.intens_", y, ".tiff"))
  if(all(unique(land.sbw$intensity) %in% 0:1))
    plot(map, col=c("grey90", "yellowgreen"), main=y) 
  if(all(unique(land.sbw$intensity) %in% 0:2))
    plot(map, col=c("grey90", "yellowgreen", "orange"), main=y) 
  if(all(unique(land.sbw$intensity) %in% 0:3))
    plot(map, col=c("grey90", "yellowgreen", "orange", "darkred"), main=y)
  dev.off()
}



