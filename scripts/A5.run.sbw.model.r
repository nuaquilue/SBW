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
r = sbw.outbreak(custom.params = NULL, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 80, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

## Saving outputs
r = sbw.outbreak(custom.params = NULL, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 2, nrun = 1, save.land = TRUE, out.seq = NA, 
                 out.path = "outputs/algo") 

## Changing default parameters of the list
source("R/default.params.r")  
custom.params = default.params()
custom.params$year.ini = 1990
r = sbw.outbreak(custom.params = custom.params, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 1, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

## Changing values of a parameters table, eg. soil.suitability of SAB
load("data/soil.suitability.rda")
soil.suitability[soil.suitability$spp=="SAB",2:6] = 0
soil.suitability
r = sbw.outbreak(custom.params = NULL, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 10, nrun = 1, save.land = FALSE, out.seq = NA, 
                 out.path = NA, soil.suitability) 

