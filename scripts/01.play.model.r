############################################ RUN SBW ##################################################
rm(list=ls())
library(sp); library(raster); library(tidyverse)

params = default.params()
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
source("mdl/sbw.outbreak.r")  
r = sbw.outbreak(custom.params = NA, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 80, nrun = 1, 
                 save.land = FALSE, out.seq = NA, out.path = NA) 
