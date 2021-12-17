############################################ DEBUG sbw.outbreak() ##################################################
rm(list=ls())
library(sp)
library(raster)
library(RANN)
library(tidyverse)
## sbw.outbreak() arguments
custom.params = NA 
rcp = NA  #rcp = "rcp45" or = "rcp85"
prec.proj = NA
temp.proj = NA  
time.horizon = 80
nrun = 1
save.land = FALSE
out.seq = NA
out.path = NA


## Decide (top-down) duration of the current outbreak
outbreak.duration <- rdunif(1, -1,2)+6+duration.last.outbreak
done <- T
t <- 1; sbw.schedule <- seq(1,80,1)




############################################ RUN sbw.outbreak() ##################################################
rm(list=ls())
library(sp)
library(raster)
library(RANN)
library(tidyverse)
source("R/sbw.outbreak.r")  
r = sbw.outbreak(custom.params = NA, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 80, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 
