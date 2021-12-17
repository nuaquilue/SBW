##########################################################################################
## 
##  kk
##
##########################################################################################

rm(list=ls())
library(sp)
library(rgdal)
library(raster)
library(viridisLite)
library(viridis)
library(tidyverse)

## Assign cell.id to values of MASK and then convert raster to polygon .shp
load("inputlyrs/rdata/mask.rdata")
MAP <- MASK
dta <- data.frame(cell.id=1:ncell(MASK), x=MASK[]) %>% filter(!is.na(x))
MAP[!is.na(MASK[])] <- dta$cell.id
plot(MAP)
poly <- rasterToPolygons(MAP)
shapefile(poly, file="C:/WORK/QBCMOD/DataIn/MaskId.shp", overwrite=T)
rm(poly); rm(MASK); rm(MAP); gc()


## Read layer "levels of defoliation" within the cells of the study area
sbw.defol.y <- data.frame(year=NA, cell.id=NA, intensity=NA, area.km=NA)
for(y in 1992:2020){
  SBWmask <-readOGR(paste0("C:/WORK/QBCMOD/DataIn/TBEdefoliation/TBEmask/TBE_", y, "mask.shp"))  
  aux <- SBWmask@data %>% filter(Ia!="0", z!=0) %>% mutate(year=ANNEE, cell.id=z, intensity=Ia) %>% 
         group_by(year, cell.id, intensity) %>% summarise(area.km=sum(area)/10^6)
  sbw.defol.y <- rbind(sbw.defol.y, as.data.frame(aux))
}
for(y in 1968:1991){
  SBWmask <-readOGR(paste0("C:/WORK/QBCMOD/DataIn/TBEdefoliation/TBEmask/TBE_", y, "mask.shp"))  
  aux <- SBWmask@data %>% filter(NIVEAU!="0", z!=0) %>% mutate(year=ANNEE, cell.id=z,
         intensity=ifelse(NIVEAU<=3, 1, ifelse(NIVEAU<=5, 2, 3))) %>%  
         group_by(year, cell.id, intensity) %>% summarise(area.km=sum(area)/10^6)
  sbw.defol.y <- rbind(sbw.defol.y, as.data.frame(aux))
}
sbw.defol.y <- sbw.defol.y[-1,]
sbw.defol.y$year <- as.numeric(sbw.defol.y$year)
save(sbw.defol.y, file="C:/work/qbcmod/DataOut/SBW.defol.y.mask.rdata")


## Count the (close to) true area defoliated by year
area.defol.y <- group_by(sbw.defol.y, year) %>% summarise(a=sum(area.km))
area.defol.y
ggplot(area.defol.y, aes(x=year, a/10^3)) + geom_line() + geom_vline(xintercept=1980, color="red")+
  geom_vline(xintercept=1990, color="red")+ geom_vline(xintercept=2000, color="red")+
  geom_vline(xintercept=2010, color="red")+   geom_hline(yintercept=8, color="blue")+
  geom_point() + scale_color_brewer(palette="Dark2")  + theme_classic()

area.defol.y

## select those cells that will be actually defoliated
## and assign the most abundant intensity level of defoliation
thresh <- data.frame(year=1968:2020, 
                     th=c(1.47, 1.9, 2.09, 1.999, 2.01, 1.928, 2.073, 2.0318, 2.0395, 2.0022, 1.9655, 1.9365,   #1968-1979
                          1.93, 1.938, 1.974, 1.999, 1.965, 1.96, 1.896, 1.88, 1.855, 1.8, 1.785, 1.78,  #1980-1991
                          1.47,0.7,0.94,1.05,.87,.9,1.07,1.235,.93,  #1992-2000
                          1.05,0.85,.7,1.18,1.29, 0.935, 1.04, 1.16, 1.25,1.4215,  #2001-2010
                          1.6075,1.6648, 1.729, 1.744,  1.75151, 1.72765, 1.7633,  1.7873, 1.9178, 1.9299)) #2011-2020
sbw.intens.y <- left_join(sbw.defol.y, thresh, by="year") %>% filter(area.km>=th) %>% 
                pivot_wider(id_cols=c(year, cell.id), names_from=intensity, values_from=area.km)
sbw.intens.y$intensity <- apply(sbw.intens.y[,3:5], 1, which.max)
cell.defol.y <- group_by(sbw.intens.y, year) %>% summarise(count=length(year)*4) %>% 
                left_join(area.defol.y, by="year") %>% mutate(dif=a-count)
filter(cell.defol.y, year==1991)


# global error
sum(cell.defol.y$dif)
# number of cell per intensity level and year
table(sbw.intens.y$year, sbw.intens.y$intensity)
## 
sbw.intens.y <- select(sbw.intens.y, year, cell.id, intensity)
save(sbw.intens.y, file="C:/work/qbcmod/DataOut/SBW.intens.y.mask.rdata")




## Maps
load("inputlyrs/rdata/mask.rdata")
load(file="inputlyrs/rdata/land.sbw.rdata")
MAP <- MASK
MAP[!is.na(MASK[])] <- land$cum.intens.def
mx.cum <- max(land$cum.intens.def)
plot(MAP, col=c(rep("grey90",15), heat.colors(mx.cum-15+1)[(mx.cum-15+1):1]), main="Accumulated level of defoliation intensity")
plot(MAP, col=c("grey90", heat.colors(mx.cum)[mx.cum:1]), xaxt="n", yaxt="n")

MAP[!is.na(MASK[])] <- land$ny.def
mx.ny <- max(land$ny.def)
plot(MAP, col=terrain.colors(mx.ny)[mx.ny:1], main="Number of consecutive years of defoliation")
plot(MAP, col=c("grey90", plasma(mx.ny)), xaxt="n", yaxt="n")

MAP[!is.na(MASK[])] <- land$curr.intens.def
plot(MAP, col=c("grey90", "yellowgreen", "orange", "darkred"), main="Current level of defoliation intensity")
plot(MAP, col=c("grey90", "yellowgreen", "orange", "darkred"), xaxt="n", yaxt="n")

## sbw.niche
land$sbw.niche <- ifelse(land$temp/10>0.5 & land$temp/10<2.8, 1,
                         ifelse(land$temp/10>-1.5 & land$temp/10<4, 2, 3))
MAP[!is.na(MASK[])] <- land$sbw.niche
tiff("c:/work/qbcmod/dataout/sbw.niche_2020.tiff")
plot(MAP, col=c("yellowgreen", "orange", "darkred"), xaxt="n", yaxt="n")
dev.off()

## sbw.niche under climate.change
clim.scn <- "rcp45"
load(file=paste0("inputlyrs/rdata/temp_", clim.scn, "_MIROC_ESM_CHEM.rdata")) 
year <- c(rep(NA,3), seq(2020,2095,5))
for(i in 4:19){
  cat(year[i], "\n")
  aux <- cc.temp[,c(1,i)] 
  names(aux) <- c("cell.id", "temp")
  land <- select(land, -temp) %>% left_join(aux, by="cell.id")
  land$sbw.niche <- ifelse(land$temp>0.5 & land$temp<2.8, 1, ifelse(land$temp>-1.5 & land$temp<4, 2, 3))
  MAP[!is.na(MASK[])] <- land$sbw.niche
  tiff(paste0("c:/work/qbcmod/dataout/sbw.niche_",  clim.scn,"_",  year[i],  ".tiff"))
  plot(MAP, col=c("yellowgreen", "orange", "darkred"), xaxt="n", yaxt="n", main=paste(clim.scn, "-", year[i]))
  dev.off()
}



## MAPS OF CURRENT DEFOLIATION FOR EACH YEAR
rm(list = ls())
load("inputlyrs/rdata/mask.rdata")
load(file="inputlyrs/rdata/land.rdata")
load(file="C:/work/qbcmod/DataOut/SBW.intens.y.mask.rdata")
dta.intens <- data.frame(year=NA, intens=NA, ncell=NA, freq=NA)
for(y in 1968:2020){
  cat(y, "\n")
  defol.y <- filter(sbw.intens.y, year==y)
  # get the information about defoliation for the target year       
  land.sbw <- left_join(land, defol.y, by="cell.id")  
  land.sbw$intensity[is.na(land.sbw$intensity)] <- 0
    # MAP <- MASK
    # MAP[!is.na(MASK[])] <- land.sbw$intensity
    # tiff(paste0("c:/work/qbcmod/dataout/sbw.defol.intens_", y, ".tiff"))
    # plot(MAP, col=c("grey90", "yellowgreen", "orange", "darkred"), main=y)
    # plot(MAP, col=c("grey90", "yellowgreen", "orange"), main=y) #2
    # plot(MAP, col=c("grey90", "yellowgreen"), main=y) #1
    # dev.off()
  aux <- data.frame(y, as.data.frame(table(land.sbw$intensity)), 
                    as.data.frame(table(land.sbw$intensity)/sum(land.sbw$intensity>0)))[,-4]
  names(aux) <- c("year", "intens", "ncell", "freq")
  dta.intens <- rbind(dta.intens, aux)
}
dta.intens <- dta.intens[-1,]
dta.intens <- filter(dta.intens, intens>0)
tot <- group_by(dta.intens, year) %>% summarise(ncell=sum(ncell))
## area defoliated
km2.cell <- 4
ggplot(tot, aes(x=year, y=ncell*km2.cell/10^3)) + geom_line() + 
  geom_point() + scale_color_brewer(palette="Dark2")  + theme_classic() +
  geom_vline(xintercept=1970, color="grey60")+geom_vline(xintercept=1980, color="grey80")+
  geom_vline(xintercept=1990, color="grey60")+geom_vline(xintercept=2000, color="grey80")+
  geom_vline(xintercept=2010, color="grey60") + geom_vline(xintercept=2020, color="grey60") + ylab("1000 km2")
  #geom_hline(yintercept=120, color="red")+ geom_hline(yintercept=2000, color="red")
## frequency intensity
ggplot(dta.intens, aes(x=year, y=freq, group=as.factor(intens))) +
  geom_line(aes(colour=as.factor(intens))) + geom_point(aes(color=as.factor(intens))) +
  scale_color_manual(values=c("yellowgreen", "orange", "darkred"))  + theme_classic()  


scn.name <- "test0"
sbw.defol.intens <- read.table(paste0("outputs/", scn.name, "/SBWdefoliation.txt"), header=T)
# years per phase
group_by(sbw.defol.intens, phase) %>% summarise(ny=length(unique(year)))
# amount defoliated
tot.predict <- filter(sbw.defol.intens, curr.intens.def>0) %>% group_by(year) %>% summarise(ncell=sum(ncell))
tot.predict <- rbind(tot, tot.predict)
tot.predict$area <- tot.predict$ncell*4/10^3
none <- filter(sbw.defol.intens, curr.intens.def==0, phase=="collapse", ncell==147464) %>% dplyr::select(year, ncell)
if(nrow(none)>0){
  none$ncell <- none$area <- 0
  tot.predict <- rbind(tot.predict, none)
}
ggplot(tot.predict, aes(x=year, y=log(area))) + geom_line() + 
  geom_point() +  theme_classic() + geom_vline(xintercept=2020, color="grey60")
ggplot(tot.predict, aes(x=year, y=area)) + geom_line() + 
  geom_point() +  theme_classic() + geom_vline(xintercept=2020, color="grey60")
# pct intensity
intens.predict <- filter(sbw.defol.intens, curr.intens.def>0) %>% select(year, curr.intens.def, ncell, pct)
names(intens.predict) <- c("year", "intens", "ncell", "freq")
intens.predict <- rbind(dta.intens, intens.predict)
ggplot(intens.predict, aes(x=year, y=freq, group=as.factor(intens))) +
  geom_line(aes(colour=as.factor(intens))) + geom_point(aes(color=as.factor(intens))) +
  scale_color_manual(values=c("yellowgreen", "orange", "darkred")) + theme_classic() +
  theme(legend.position = "none") + geom_vline(xintercept=2020, color="grey60")



## new incorporations versus current 
load(file="C:/work/qbcmod/DataOut/SBW.intens.y.mask.rdata")
dta <- data.frame(year=NA, previous=NA, ncell=NA, pct=NA)
for(y in 1969:2020){
  cat(y, "\n")
  defol.y <- filter(sbw.intens.y, year==(y-1))
  tot <- filter(sbw.intens.y, year==y) %>% group_by(year) %>% summarise(x=length(year))
  aux <- filter(sbw.intens.y, year==y) %>% mutate(previous=cell.id %in% defol.y$cell.id) %>% 
    group_by(year, previous) %>% summarise(ncell=length(year)) %>% left_join(tot, by="year") %>% 
    mutate(pct=round(ncell/x*100,1)) %>% dplyr::select(-x)
  dta <- rbind(dta, as.data.frame(aux))
}
dta <- dta[-1,]
dta$status <- ifelse(dta$previous, "already", "new")
dta.out <- filter(dta, year>=1993 & year<=2011) %>% filter(year!=2000)
dta.out <- filter(dta, year>=2011)

summary(dta.out$pct[dta.out$previous])
summary(dta.out$pct[!dta.out$previous])

ggplot(dta, aes(x=year, y=pct, group=status)) + geom_line(aes(colour=status)) +
  geom_point(aes(color=status)) + scale_color_brewer(palette="Dark2")  + theme_classic() +
  geom_vline(xintercept=1970, color="grey60")+geom_vline(xintercept=1980, color="grey80")+
  geom_vline(xintercept=1990, color="grey60")+geom_vline(xintercept=2000, color="grey80")+
  geom_vline(xintercept=2010, color="grey60")


tot <- group_by(dta, year) %>% summarise(ncell=sum(ncell))




## count
filter(sbw, cum.intens.def>=15) %>% group_by(spp) %>% summarise(area=length(spp)*4)
filter(sbw, spp=="SAB", cum.intens.def>=15) %>% group_by(spp, ny.def) %>% summarise(area=length(spp)*4)

kk <- filter(sbw, spp=="SAB", cum.intens.def>=15) %>% group_by(ny.def, cum.intens.def) %>% 
      summarise(area=length(spp)*4) 
kk$class.ny <- NA
kk$class.ny[kk$ny.def<5] <- "M0"
kk$class.ny[kk$ny.def==5] <- "M1"
kk$class.ny[kk$ny.def %in% c(6,7) & kk$cum.intens.def>=17] <- "M1"
kk$class.ny[kk$ny.def %in% c(6,7) & kk$cum.intens.def<17] <- "M2"
kk$class.ny[kk$ny.def %in% c(8,9) & kk$cum.intens.def>=20] <- "M1"
kk$class.ny[kk$ny.def %in% c(8,9) & kk$cum.intens.def<20] <- "M2"
kk$class.ny[kk$ny.def<10 & kk$cum.intens.def>20] <- "M1"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def<=20] <- "M3"
kk$class.ny[kk$ny.def>=12 & kk$cum.intens.def<=25] <- "M3"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def>25 & kk$cum.intens.def<=30] <- "M4"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def>=30] <- "M5"

ggplot(kk, aes(as.factor(ny.def), cum.intens.def)) +
  geom_point(aes(size=area, colour=class.ny)) + theme_bw()

group_by(kk, class.ny) %>% summarise(area=sum(area))

a <- filter(sbw, spp=="SAB", cum.intens.def>=15)
 table(a$ny.def, a$cum.intens.def) *4
 b <- table(a$ny.def)*4; b
 plot(names(b), as.numeric(b), xlab="ny.def", ylab="area",
      main="Balsam fir area with accumulated defoliation intensity >=15")
