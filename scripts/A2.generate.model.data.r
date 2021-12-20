##########################################################################################
## 
##  Generation of data for SBW model: tables, raster of the study area, and completed
##  'land' data.frame
##
##########################################################################################

rm(list = ls())
library(sp)
library(raster)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select

###### TABLES
input.path = "C:/WORK/OneDrive - ctfc.cat/QBCMOD/QbcLDM/inputfiles/"
post.sbw.reg = read.table(paste0(input.path, "PostSBWRege.txt"), header=T)
forest.succ = read.table(paste0(input.path,"ForestSucc.txt"), header=T)
spp.colonize.persist = read.table(paste0(input.path,"SppColonizePersist.txt"), header=T)
temp.suitability = read.table(paste0(input.path,"ThMeanTemp.txt"), header=T)  
prec.suitability = read.table(paste0(input.path,"ThAnnualPrecip.txt"), header=T)  
soil.suitability = read.table(paste0(input.path,"ThSoil.txt"), header=T)  
default.tables = list(post.sbw.reg=post.sbw.reg, forest.succ=forest.succ, 
           spp.colonize.persist=spp.colonize.persist, temp.suitability=temp.suitability,
           prec.suitability=prec.suitability, soil.suitability=soil.suitability)
usethis::use_data(default.tables, overwrite = T)


###### RASTER MASK: Rename raster of the study area, from MASK to mask
input.path = "c:/work/OneDrive - ctfc.cat/QBCMOD/QbcLDM/inputlyrs/"
load(paste0(input.path, "rdata/mask.rdata"))
mask = MASK
usethis::use_data(mask, overwrite = T)


###### CLIMATE: Rename data frames with climatic projections
input.path = "c:/work/OneDrive - ctfc.cat/QBCMOD/QbcLDM/inputlyrs/"
load(paste0(input.path, "rdata/prec_rcp45_MIROC_ESM_CHEM.rdata"))
prec_rcp45 = cc.prec
usethis::use_data(prec_rcp45, overwrite = T)
load(paste0(input.path, "rdata/prec_rcp85_MIROC_ESM_CHEM.rdata"))
prec_rcp85 = cc.prec
usethis::use_data(prec_rcp85, overwrite = T)
load(paste0(input.path, "rdata/temp_rcp45_MIROC_ESM_CHEM.rdata"))
temp_rcp45 = cc.temp
usethis::use_data(temp_rcp45, overwrite = T)
load(paste0(input.path, "rdata/temp_rcp85_MIROC_ESM_CHEM.rdata"))
temp_rcp85 = cc.temp
usethis::use_data(temp_rcp85, overwrite = T)


###### ELEVATION: Add elevation to 'land' data frame
## Read original 'land' and updated raster mask of the study area
input.path = "c:/work/OneDrive - ctfc.cat/QBCMOD/QbcLDM/inputlyrs/"
load(file = paste0(input.path, "rdata/land.rdata"))
load(file = "data/mask.rda")

## Change projection of Elevation and clip it to the study area
ELEV = raster(paste0(dirname(dirname(getwd())), "/DataIn/Elevation/elev_QC3.grd"))
crs(ELEV)
ELEVlcc = projectRaster(ELEV, res=2000, crs=crs(mask))
ELEVlcc.clip = crop(ELEVlcc, mask)
land$elev = ELEVlcc.clip[!is.na(mask[])]

## Assign mean elevation of the neigbour cells to those cells with non informed values
zcells = filter(land, !is.na(spp) & is.na(elev))
r = 3
while(nrow(zcells)>0){
  neighs = RANN::nn2(select(land, x, y), select(zcells, x,y), searchtype="priority", k=r^2)
  values = matrix(land$elev[neighs$nn.idx], ncol=r^2)
  land$elev[!is.na(land$spp) & is.na(land$elev)] = apply(values, 1, mean, na.rm=T)
  zcells = filter(land, !is.na(spp) & is.na(elev))
  r = r+2
}

## Also, remove unecessary fields in 'land' and 
land = select(land, -tsfire, -tsccut, - tspcut)

## And Temperature in 'land' is x 10, but not in the climatic projections
land$temp = land$temp/10 

## Now compute the accumulated intensity of defoliation, the number of years since the first defoliation,
## and the number of years without defoliation to reset defoliation once it is greater or equal than 5
load(file=paste0(dirname(dirname(getwd())), "/DataOutSBW/sbw.defol.intens/SBW.intens.y.mask.rdata"))
land$ny.def = 0 
land$ny.def0 = 0
land$cum.intens.def = 0
for(y in 1992:2020){
  defol.y = filter(sbw.intens.y, year==y)
  # get the information about defoliation for the target year       
  land = left_join(land, defol.y, by="cell.id")  
  # Increment number of years of defoliation whenever current intensity is >0
  land$ny.def[!is.na(land$intensity)] = land$ny.def[!is.na(land$intensity)]+1
  # Mark that the number of years without defoliation is 0 too
  land$ny.def0[!is.na(land$intensity)] =  0
  # Increment number of years of NO defoliation whenever intensity is NA
  land$ny.def0[is.na(land$intensity)] = land$ny.def0[is.na(land$intensity)]+1
  # Reset number of years since defoliation to 0 when there is no defoliation in the current year
  # and in the 4 preceding years
  land$ny.def[is.na(land$intensity) & land$ny.def0>=5] = 0
  # Increment cumulative intensity whenever current intensity is >0
  land$cum.intens.def[!is.na(land$intensity)] = land$cum.intens.def[!is.na(land$intensity)]+land$intensity[!is.na(land$intensity)]
  # Reset cumulative intenstiy when there is no defoliation in the current year
  # and in the 4 preceding years
  land$cum.intens.def[is.na(land$intensity) & land$ny.def0>=5] = 0
  land = select(land, -intensity, -year)
}
## current level of defoliation
land = left_join(land, defol.y, by="cell.id")  
land$curr.intens.def[!is.na(land$intensity)] = land$intensity[!is.na(land$intensity)]
land$curr.intens.def[is.na(land$intensity)] = 0
land = select(land, -intensity, -year)

## Rename 'land' to 'landscape' to be ready to use in sbw.outbreak() and
## save it in the 'data' folder.
landscape = land
usethis::use_data(landscape, overwrite = T)  # Forest and sbw outbreak data
