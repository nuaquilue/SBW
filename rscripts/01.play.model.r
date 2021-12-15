############################################ RUN SBW ##################################################
rm(list=ls())
library(sp); library(raster); library(tidyverse)
source("mdl/define.scenario.r")
scn.name <- "test"; define.scenario(scn.name)
source(paste0("outputs/", scn.name, "/scn.def.r"))
source("mdl/sbw.outbreak.r")
source("mdl/neigh.influence.sbw.spread.r")
load(file="inputlyrs/rdata/mask.rdata")
load(file="inputlyrs/rdata/land.sbw.rdata")
land$temp <- land$temp/10
`%notin%` <- Negate(`%in%`)
scn.name <- "test"
is.sbw <- T
write.maps <- T
## Decide (top-down) duration of the current outbreak
outbreak.duration <- rdunif(1, -1,2)+6+duration.last.outbreak
done <- T
t <- 1; sbw.schedule <- seq(1,80,1)




############################################ RUN A SCN ##################################################
rm(list=ls())
source("mdl/define.scenario.r"); source("mdl/landscape.dyn.r")  
scn.name <- "test0"
define.scenario(scn.name)
nrun <- 1
write.maps <- T
time.horizon <- 80
is.wildfires <-  is.clearcut <- is.partialcut <- F
is.sbw <- T
dump(c("nrun", "write.maps",  "time.horizon", "is.wildfires", "is.clearcut", "is.partialcut", "is.sbw"), 
     paste0("outputs/", scn.name, "/scn.custom.def.r"))
landscape.dyn(scn.name)


############################################ RUN A SET OF SCN ##################################################
library(readxl)
rm(list=ls())
source("mdl/define.scenario.r"); source("mdl/landscape.dyn.r")  
scenarios <- read_xlsx("Scenarios.xlsx", sheet="ScnPaper")
for(i in 17){
  scn.name <- scenarios$scn.name[i]
  define.scenario(scn.name)
  ## general
  nrun <- scenarios$nrun[i]
  write.maps <- F
  plot.fires <- F
  ## processes
  is.wildfires <- as.logical(scenarios$is.wildfires[i])
  is.clearcut <- as.logical(scenarios$is.clearcut[i])
  is.partialcut <- as.logical(scenarios$is.partialcut[i])
  ## modifiers of target area
  is.fuel.modifier <- as.logical(scenarios$is.fuel.modifier[i])
  is.clima.modifier <- as.logical(scenarios$is.clima.modifier[i])
  ## scenario parameters
  clim.scn <- ifelse(scenarios$clim.scn[i]=="NA", NA, scenarios$clim.scn[i])
  replanif <- as.logical(scenarios$replanif[i])
  th.small.fire <- scenarios$th.small.fire[i]
  wflam <- scenarios$wflam[i]
  wwind <- scenarios$wwind[i]
  pigni.opt <- scenarios$pigni[i]
  target.old.pct <- scenarios$target.old.pct[i]
  dump(c("nrun", "write.maps", "plot.fires", "is.wildfires", "is.clearcut", "is.partialcut", 
         "is.fuel.modifier", "is.clima.modifier", "clim.scn", "replanif",
         "th.small.fire", "wflam", "wwind", "pigni.opt", "target.old.pct"), 
       paste0("outputs/", scn.name, "/scn.custom.def.r"))
  landscape.dyn(scn.name)
}



##################################### BUILD INITIAL CONDITIONS ###########################################
rm(list=ls())
work.path <- "C:/WORK/QBCMOD/"
source("mdl/read.state.vars.r")
read.state.vars(work.path)
source("mdl/build.pigni.r")
build.pigni(work.path, lambda=0.04, r=0.95, first.time=F)



