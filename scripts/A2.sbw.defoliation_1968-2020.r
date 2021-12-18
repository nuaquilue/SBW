##########################################################################################
## 
##  Read the TBEdefoliation shape files donwloaded from the MFFP 
##  to build data frames with cell-level defoliation intensity per year
##
##########################################################################################

rm(list=ls())
library(sp)
library(raster)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select


## Assign cell.id to values of mask and then convert raster to a polygon 
## and save as a Shapefile
load("data/mask.rda")
map = mask
dta = data.frame(cell.id=1:ncell(mask), x=mask[]) %>% filter(!is.na(x))
map[!is.na(mask[])] = dta$cell.id
poly = rasterToPolygons(map)
shapefile(poly, file="DataIn/MaskId.shp", overwrite=T)
rm(poly); rm(mask); rm(MAP); gc()


## Read layer "levels of defoliation" within the cells of the study area 
## Build a data frame with defolition intensity at the cell level per each year
sbw.defol.y = data.frame(year=NA, cell.id=NA, intensity=NA, area.km=NA)
for(y in 1992:2020){
  SBWmask = shapefile(paste0(dirname(getwd()), "/DataIn/TBEdefoliation/TBEmask/TBE_", y, "mask.shp"))  
  aux = SBWmask@data %>% filter(Ia!="0", z!=0) %>% mutate(year=ANNEE, cell.id=z, intensity=Ia) %>% 
         group_by(year, cell.id, intensity) %>% summarise(area.km=sum(area)/10^6)
  sbw.defol.y = rbind(sbw.defol.y, as.data.frame(aux))
}
for(y in 1968:1991){
  SBWmask = shapefile(paste0(dirname(getwd()), "/DataIn/TBEdefoliation/TBEmask/TBE_", y, "mask.shp"))  
  aux = SBWmask@data %>% filter(NIVEAU!="0", z!=0) %>% mutate(year=ANNEE, cell.id=z,
         intensity=ifelse(NIVEAU<=3, 1, ifelse(NIVEAU<=5, 2, 3))) %>%  
         group_by(year, cell.id, intensity) %>% summarise(area.km=sum(area)/10^6)
  sbw.defol.y = rbind(sbw.defol.y, as.data.frame(aux))
}
sbw.defol.y = sbw.defol.y[-1,]
sbw.defol.y$year = as.numeric(sbw.defol.y$year)
save(sbw.defol.y, file=paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.defol.y.mask.rdata"))


## Count the (close to) true area defoliated by year and plot it
load(file=paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.defol.y.mask.rdata"))
area.defol.y = group_by(sbw.defol.y, year) %>% summarise(a=sum(area.km))
area.defol.y
ggplot(area.defol.y, aes(x=year, a/10^3)) + geom_line() + geom_vline(xintercept=1980, color="red")+
  geom_vline(xintercept=1990, color="red")+ geom_vline(xintercept=2000, color="red")+
  geom_vline(xintercept=2010, color="red")+   geom_hline(yintercept=8, color="blue")+
  geom_point() + scale_color_brewer(palette="Dark2")  + theme_classic()
area.defol.y


## Select those cells that will be actually defoliated
## and assign the most abundant intensity level of defoliation
thresh = data.frame(year=1968:2020, 
                     th=c(1.47, 1.9, 2.09, 1.999, 2.01, 1.928, 2.073, 2.0318, 2.0395, 2.0022, 1.9655, 1.9365,   #1968-1979
                          1.93, 1.938, 1.974, 1.999, 1.965, 1.96, 1.896, 1.88, 1.855, 1.8, 1.785, 1.78,  #1980-1991
                          1.47,0.7,0.94,1.05,.87,.9,1.07,1.235,.93,  #1992-2000
                          1.05,0.85,.7,1.18,1.29, 0.935, 1.04, 1.16, 1.25,1.4215,  #2001-2010
                          1.6075,1.6648, 1.729, 1.744,  1.75151, 1.72765, 1.7633,  1.7873, 1.9178, 1.9299)) #2011-2020
sbw.intens.y = left_join(sbw.defol.y, thresh, by="year") %>% filter(area.km>=th) %>% 
                pivot_wider(id_cols=c(year, cell.id), names_from=intensity, values_from=area.km)
sbw.intens.y$intensity = apply(sbw.intens.y[,3:5], 1, which.max)
cell.defol.y = group_by(sbw.intens.y, year) %>% summarise(count=length(year)*4) %>% 
                left_join(area.defol.y, by="year") %>% mutate(dif=a-count)

# global error
sum(cell.defol.y$dif)
# number of cell per intensity level and year
table(sbw.intens.y$year, sbw.intens.y$intensity)
## Data frame with cell-level defoliation intensity per year
sbw.intens.y = select(sbw.intens.y, year, cell.id, intensity)
save(sbw.intens.y, file=paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.intens.y.mask.rdata"))


