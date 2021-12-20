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

## Saving outputs
r = sbw.outbreak(custom.params = NULL, custom.tables = NULL, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 80, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

## Save outputs
r = sbw.outbreak(custom.params = NULL, custom.tables = NULL, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 2, nrun = 1, save.land = TRUE, out.seq = NA, out.path = "outputs/algo2") 

## Change default parameters of the list
source("R/default.params.r")  
custom.params = default.params()
custom.params$year.ini = 1990
r = sbw.outbreak(custom.params = custom.params, custom.tables = NULL, 
                 rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 1, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

## Change values of an innput table, eg. soil.suitability of SAB
data(default.tables)
soil.suitability = tbl[["soil.suitability"]]
soil.suitability[soil.suitability$spp=="SAB",2:6] = 0
tbl[["soil.suitability"]] = soil.suitability
r = sbw.outbreak(custom.params = NULL, custom.tables = tbl, rcp = NA, prec.proj = NA, temp.proj = NA,  
                 time.horizon = 10, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA) 

