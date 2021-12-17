##########################################################################################
## 
##  Function that reads BioSim outputs containing monthly climatic projections
##  for the study area of 7 variables: min temp, max temp, mean temp, accum precip
##  mean wind speed, max wind speed, and degree days 5ºC.
##  For each year and variable a raster brick is built and saved as a .rds file
##
##########################################################################################

rm(list=ls())
library(tidyverse)
options(dplyr.summarise.inform=F)
select = dplyr::select


############ Function to plot age class distribution per species over time ############
ageclass.spp = function(scn.name){
  
  ## Read abundance of species per age class and reclassify in 
  ## 'young', 'mature', and 'old'.
  spp.age <- read.table(paste0("outputs/", scn.name, "/SppByAgeClass.txt"), header=T)
  spp.age.year <- mutate(spp.age, ymo=ifelse(age.class %in% c("C10", "C30", "C50"), "young",
    ifelse(age.class %in% c("C70", "C90"), "mature", "old"))) %>% 
    group_by(run, year, spp, ymo) %>% summarise(area=sum(area)) %>% 
    group_by(year, spp, ymo) %>% summarise(area=mean(area)) 
  
  ## Plot and save
  ylab1 = expression(paste(10^3 %.% km^2))
  p = ggplot(spp.age.year, aes(x=year, y=area/10^3, group=ymo)) + geom_line(aes(colour=ymo)) +
      geom_point(aes(color=ymo)) + scale_color_brewer(palette="Dark2") + theme_classic() +
      facet_wrap(~ spp, ncol=5, nrow=2) +  scale_y_continuous(ylab1)
  ggsave(p, filename = paste0(dirname(getwd()), "/DataOutSBW/mdl.outs/ageclass.spp_", scn.name, ".png"),
         width = 10, height = 6)
}


############ Function to plot age class distribution per species over time ############
defol.intens.phase = function(scn.name){
  
  ## Read model data
  track.sbw.defol.intens <- read.table(paste0("outputs/", scn.name, "/SBWdefoliation.txt"), header=T)  
  
  ## Plot % defoliation intenstiy per outbreak phase
  p1 = ggplot(filter(track.sbw.defol.intens, !is.na(pct)), aes(x=year, y=pct, group=as.factor(curr.intens.def))) +
       geom_line(aes(colour=as.factor(curr.intens.def))) + geom_point(aes(color=as.factor(curr.intens.def))) +
       scale_color_brewer(palette="Dark2") + theme_classic() + facet_wrap(~ phase)  
  ggsave(p1, filename = paste0(dirname(getwd()), "/DataOutSBW/mdl.outs/pct.defol.phase_", scn.name, ".png"),
         width = 10, height = 6)
  
  ## Plot of number of cells defoliated per year
  ylab1 = expression(paste(10^3 %.% ncell))
  p2 = ggplot(filter(track.sbw.defol.intens, !is.na(pct)), aes(x=year, y=ncell/10^3, group=as.factor(curr.intens.def))) +
       geom_line(aes(colour=as.factor(curr.intens.def))) + geom_point(aes(color=as.factor(curr.intens.def))) +
       scale_color_brewer(palette="Dark2") + theme_classic()  +   scale_y_continuous(ylab1)
  ggsave(p2, filename = paste0(dirname(getwd()), "/DataOutSBW/mdl.outs/ncell.defol.year_", scn.name, ".png"),
         width = 10, height = 6)
}




############ Execute outputs processing functions ############
scn.name = "Test80"
ageclass.spp(scn.name)
defol.intens.phase(scn.name)
