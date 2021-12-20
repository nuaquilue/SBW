##########################################################################################
## 
##  Explore sbw defoliation data
##
##########################################################################################

rm(list=ls())
library(sp)
library(raster)
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select

## A. Compute the number of cells per level of defoliation
## and the percentage it represents, per each year
load("data/land.sbw.rda")
load(paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.intens.y.mask.rdata"))
dta.intens = data.frame(year=NA, intens=NA, ncell=NA, freq=NA)
for(y in 1968:2020){
  cat(y, "\n")
  defol.y = filter(sbw.intens.y, year==y)
  # get the information about defoliation for the target year       
  land.sbw = left_join(landscape, defol.y, by="cell.id")  
  land.sbw$intensity[is.na(land.sbw$intensity)] = 0
  aux = data.frame(y, as.data.frame(table(land.sbw$intensity)), 
                   as.data.frame(table(land.sbw$intensity)/sum(land.sbw$intensity>0)))[,-4]
  names(aux) = c("year", "intens", "ncell", "freq")
  dta.intens = rbind(dta.intens, aux)
}
dta.intens = dta.intens[-1,]
dta.intens = filter(dta.intens, intens>0)


## B. Plot a line graph with annual area defoliated 
tot = group_by(dta.intens, year) %>% summarise(ncell=sum(ncell))
km2.cell = 4
ggplot(tot, aes(x=year, y=ncell*km2.cell/10^3)) + geom_line() + 
  geom_point() + scale_color_brewer(palette="Dark2")  + theme_classic() +
  geom_vline(xintercept=1970, color="grey60")+geom_vline(xintercept=1980, color="grey80")+
  geom_vline(xintercept=1990, color="grey60")+geom_vline(xintercept=2000, color="grey80")+
  geom_vline(xintercept=2010, color="grey60") + geom_vline(xintercept=2020, color="grey60") + ylab("1000 km2")

## Plot a line graph with frequency of defoliation intensity
ggplot(dta.intens, aes(x=year, y=freq, group=as.factor(intens))) +
  geom_line(aes(colour=as.factor(intens))) + geom_point(aes(color=as.factor(intens))) +
  scale_color_manual(values=c("yellowgreen", "orange", "darkred"))  + theme_classic()  


## C. Compare pair of consecutive years to count new incorporations 
load(paste0(dirname(getwd()), "/DataOutSBW/sbw.defol.intens/SBW.intens.y.mask.rdata"))
dta = data.frame(year=NA, previous=NA, ncell=NA, pct=NA)
for(y in 1969:2020){
  cat(y, "\n")
  defol.y = filter(sbw.intens.y, year==(y-1))
  tot = filter(sbw.intens.y, year==y) %>% group_by(year) %>% summarise(x=length(year))
  aux = filter(sbw.intens.y, year==y) %>% mutate(previous=cell.id %in% defol.y$cell.id) %>% 
    group_by(year, previous) %>% summarise(ncell=length(year)) %>% left_join(tot, by="year") %>% 
    mutate(pct=round(ncell/x*100,1)) %>% dplyr::select(-x)
  dta = rbind(dta, as.data.frame(aux))
}
dta = dta[-1,]
dta$status = ifelse(dta$previous, "already", "new")
dta.out = filter(dta, year>=1993 & year<=2011) %>% filter(year!=2000)
dta.out = filter(dta, year>=2011)

##
summary(dta.out$pct[dta.out$previous])
summary(dta.out$pct[!dta.out$previous])

## Plot the annual percentage of new defoliated cells vs. cells already defoliated
ggplot(dta, aes(x=year, y=pct, group=status)) + geom_line(aes(colour=status)) +
  geom_point(aes(color=status)) + scale_color_brewer(palette="Dark2")  + theme_classic() +
  geom_vline(xintercept=1970, color="grey60")+geom_vline(xintercept=1980, color="grey80")+
  geom_vline(xintercept=1990, color="grey60")+geom_vline(xintercept=2000, color="grey80")+
  geom_vline(xintercept=2010, color="grey60")


## D. Plot of mortality groups
load("data/land.sbw.rda")
filter(landscape, cum.intens.def>=15) %>% group_by(spp) %>% summarise(area=length(spp)*4)
filter(landscape, spp=="SAB", cum.intens.def>=15) %>% group_by(spp, ny.def) %>% summarise(area=length(spp)*4)
kk = filter(landscape, spp=="SAB", cum.intens.def>=15) %>% group_by(ny.def, cum.intens.def) %>% 
  summarise(area=length(spp)*4) 
kk$class.ny = NA
kk$class.ny[kk$ny.def<5] = "M0"
kk$class.ny[kk$ny.def==5] = "M1"
kk$class.ny[kk$ny.def %in% c(6,7) & kk$cum.intens.def>=17] = "M1"
kk$class.ny[kk$ny.def %in% c(6,7) & kk$cum.intens.def<17] = "M2"
kk$class.ny[kk$ny.def %in% c(8,9) & kk$cum.intens.def>=20] = "M1"
kk$class.ny[kk$ny.def %in% c(8,9) & kk$cum.intens.def<20] = "M2"
kk$class.ny[kk$ny.def<10 & kk$cum.intens.def>20] = "M1"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def<=20] = "M3"
kk$class.ny[kk$ny.def>=12 & kk$cum.intens.def<=25] = "M3"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def>25 & kk$cum.intens.def<=30] = "M4"
kk$class.ny[kk$ny.def>=10 & kk$cum.intens.def>=30] = "M5"
ggplot(kk, aes(as.factor(ny.def), cum.intens.def)) +
  geom_point(aes(size=area, colour=class.ny)) + theme_bw()
group_by(kk, class.ny) %>% summarise(area=sum(area))


## E. 
km2.cell = 4
a = filter(landscape, spp=="SAB", cum.intens.def>=15)
b = table(a$ny.def)*km2.cell
plot(names(b), as.numeric(b), xlab="ny.def", ylab="area",
     main="Balsam fir area with accumulated defoliation intensity >=15")
