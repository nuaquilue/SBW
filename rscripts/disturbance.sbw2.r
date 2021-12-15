#+eval=FALSE

######################################################################################
###  disturbance.sbw()
###
###  Description >  Simulates spruce budworm mortality according to a given probability 
###   of mortality after defoliation.
###
###  Arguments >  
###   subland : appropiate selection fo the data frame of the state variables
###   sbw.mortality : data frame of the probability of mortality by bioclimatic domain 
###                   and ecological type
###   disturb.interact : data frame of the time after a distrubance cap happen following
###                      another disturbance 
###   write.tbl.outputs : if TRUE
###   km2.pixel : number of km2 per pixel on the grid 
###   irun : the current replica (used when writing results)
###   t : the current time step  (used when writing results)
###   out.path : directory path to save output text files
###   out.overwrite : if TRUE the output text files are overwritten 
###
###  Details > Trees mortality by spruce budworm defoliation depends on the severity
###            of the epidemics.
###
###  Value >  A vector of the indexes of the killed cells.
######################################################################################

# subland <- subset(land, select = -c(FRZone, MgmtUnit, TSW, Exclus))
# subland <- subset(land, select=c(cell.indx, SDOM_BIO, SppGrp, Temp, TSF, TSD), SppGrp!="NonFor") 
# subland <- land
# subland <- subset(land, select=c(cell.indx, SDOM_BIO, SppGrp, SBWgrowth, TSF))

disturbance.sbw <- function(subland, write.tbl.outputs = TRUE, km2.pixel = km2.pixel, irun, t, out.path = NULL, out.overwrite = TRUE, epn.pheno){
  # Silence
  #options(warn=-1)
  

  # Probability estimated from glm regression model calibrated on SIFORT data
  # glm(as.factor(outbreak) ~ sp_group_1 + extSBW, data = prg123_f, family = binomial(link=logit))
  prob <- ifelse(subland[, "SppGrp"] == "EPN" & subland$TSF !=0, 
                 1/(1 + exp(-(TBEpred_coefs[,1] + TBEpred_coefs[,2] * epn.pheno + TBEpred_coefs[,3] * subland[, "SBWgrowth"]))),
                 
                 ifelse(subland[, "SppGrp"] == "SAB" & subland$TSF !=0, 
                        1/(1 + exp(-(TBEpred_coefs[,1] + TBEpred_coefs[,2] * 0 + TBEpred_coefs[,3] * subland[, "SBWgrowth"]))), 0))
  
  subland$kill <- prob > runif(length(prob))
  
  # Write the results (area killed by domain and species group) in an output text file
  if(write.tbl.outputs){
    aux <- data.frame(sdomain = subland$SDOM_BIO[subland$kill==1],
                      spp.grp = subland$SppGrp[subland$kill==1], 
                      kill.area = km2.pixel)
    aux <- aggregate(aux$kill.area, list(sdomain = aux$sdomain, spp.grp = aux$spp.grp), sum)
    names(aux)[3] <- "kill.area"
    write.table(data.frame(run=irun, time=t, aux),
                file=paste0(out.path, "/SBWmortality.txt"), append =! out.overwrite, quote=F,
                sep="\t", row.names=F , col.names=out.overwrite)
  }
  
  load(file=paste0("inputlyrs/rdata/mask", name.resol, ".rdata"))
  
  MASK[!is.na(MASK)] <- prob
   
  print(range(as.data.frame(MASK), na.rm = TRUE))
  print(paste("Mean SBW mortality probability - ", round(mean(prob), digits = 2)))  
  print(paste("BS phenology coef:", epn.pheno))
  writeRaster(MASK, paste0(out.path, "/asc/TBE_", irun, "_",t, ".asc"), format="ascii", overwrite = TRUE)
  # Return the sbw killed cell.indx 
  return(subland$cell.indx[subland$kill==1])
  
}
