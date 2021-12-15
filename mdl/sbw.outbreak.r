sbw.outbreak = function(scn.name){
  
  ## Load required functions (no need in a R-package)
  source("mdl/buffer.mig.r") 
  source("mdl/forest.transitions.r")  
  source("mdl/sbw.outbreak.r")
  source("mdl/neigh.influence.sbw.spread.r")
  source("mdl/select.others.r")
  source("mdl/suitability.r") 

  
  ## Function to select items not in a vector
  `%notin%` = Negate(`%in%`)
  


  
  ## Set the directory for writing spatial outputs (create if it does not exist yet) 
  if(write.maps){      
    if(!file.exists(paste0(out.path, "/lyr")))
       dir.create(file.path(getwd(), out.path, "/lyr"), showWarnings = F) 
  }

  
  ## Load MASK raster layer of the study area, and compute cell resolution in km2
  load(file="inputlyrs/rdata/mask.rdata")
  km2.pixel = res(MASK)[1] * res(MASK)[2] / 10^6
  
  
  ## Load temperature and precipitation 5-year predictions according to the climatic scenario.
  ## If climate change is not activated, initial temp and precip will be used for the whole simulation.
  if(!is.na(clim.scn)){
    load(file=paste0("inputlyrs/rdata/temp_", clim.scn, "_MIROC_ESM_CHEM.rdata")) 
    load(file=paste0("inputlyrs/rdata/prec_", clim.scn, "_MIROC_ESM_CHEM.rdata"))  
  }
  
  ## GENERAL and VEGETATION DYNAMICS parameters:
  year.ini = 2020
  enable.succ = TRUE
  enfeuil = 0.0
  age.seed = 40     
  p.failure = 0     
  suboptimal = 0.5 

  ## Build the discrete time sequence according to time.step
  time.seq = seq(1, time.horizon, 1) 
  ## Set the scheduling of the processes
  sbw.schedule = seq(sbw.step, time.horizon, sbw.step)
  clim.schedule = seq(clim.step, time.horizon, clim.step)
  
  ## Tracking data.frames 
  track.sbw.defol.intens = NULL #data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA)
  track.sbw.kill = NULL #data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA)
  
  ## Start the simulations
  irun = 1
  # for(irun in 1:nrun){
    
    ## Load dynamic state variables 
    load(file="inputlyrs/rdata/land.sbw.rdata")
    
    ## Temperature in 'land' is x 10, but not in the climatic projections
    land$temp = land$temp/10
    
    
    
    ## Decide (top-down) duration of the current outbreak, epidimic phase and collapse
    outbreak = 12 - current.duration  # testing
    outbreak = rdunif(1,-1,1)+6+duration.last.outbreak - current.duration
    duration.last.outbreak = outbreak + current.duration
    phase = "outbreak"; done = T
    
    
    ## Record initial distributions:
    ncell.def = c(NA, rep(sum(land$curr.intens.def>0),3))
    track.sbw.defol.intens = rbind(track.sbw.defol.intens, 
      data.frame(run=irun, year=year.ini, phase=phase, group_by(land, curr.intens.def) %>% 
                   summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def)))
    
    
    ## Start 
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      print(paste0("scn: ", scn.name, " - run: ", irun, "/", nrun, " - time: ", t, "/", time.horizon, ", from ", year.ini+t-time.step, " to ", t+year.ini))
      
      ## Update climatic variables at each time step if climate change is activated
      ## Column 1 is cell.index, the following columns are temp (precip) in 
      ## 2020-2025, 2025-2030, 2030-2035, ... etc.
      ## The last column (temp95) then corresponds to the period 2095-2100
      ## The first time step (t=5) we start with climate 2020-2025
      if(!is.na(clim.scn) & t %in% clim.schedule){
        # Temp
        aux = cc.temp[,c(1,3+which(clim.schedule==t))] 
        names(aux) = c("cell.id", "temp")
        land = select(land, -temp) %>% left_join(aux, by="cell.id")
        # Precip
        aux = cc.prec[,c(1,3+which(clim.schedule==t))]
        names(aux) = c("cell.id", "prec")
        land = select(land, -prec) %>% left_join(aux, by="cell.id")
       }

      

      
        
      ## 2. SBW (under development)
      kill.cells = integer()
      if(is.sbw & t %in% sbw.schedule){
        cat("SBW outbreak: ")
        
        ## 2.1. Select new cells to be defoliated
        if(preoutbreak>0){
          potential = filter(land, ny.def0>=5, tssbw>=30, temp>0.5, temp<2.8)
          # preoutbreak=5
          sbw.new.sprd = sample(potential$cell.id, size=rdunif(1,1,6-preoutbreak), replace=F,
                                 prob=potential$ny.def0*(1200-potential$elev)/100)
          ## find between 20 to 40 neighs and add to sbw.new.sprd
          neighs = nn2(select(land,x,y), filter(land, cell.id %in% sbw.new.sprd) %>% select(x,y),
                    k=rdunif(1,20,40) , searchtype='priority')
          # neighs = nn2(select(land,x,y), filter(land, cell.id %in% sbw.new.sprd) %>% select(x,y),
          #               k=40 , searchtype='priority')  #test
          nn.indx = neighs[[1]]
          sbw.new.sprd = land$cell.id[nn.indx]
          # just give some intensity to the core of the epicenters
          land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
            sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
          # add some neighs of these epicenters
          sbw.new.sprd = c(sbw.new.sprd, 
              sbw.spread.tonew(land, nc=ncol(MASK), side=res(MASK)[1]/10^3, 
                               radius=12, outbreak, preoutbreak))
          # and finally assign intensity to all of them (rewrite intensity just assigned to epicenter cores)
          land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
            sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
        }
        else if(outbreak>0){
          
          # once in a while vary the radius, to allow spreading furhter away, or reduce it
          # to limit outbreak....
          
          ## 1. Spatial spreading of the current outbreak to cells not yet defoliated, that is,
          ## cells with ny.def0>=5 & tssbw>=30
          ## The function 'sbw.spread.tonew' returns cell.ids.
          radius = rdunif(1,2,15) # 4 to 60 km
          sbw.new.sprd = sbw.spread.tonew(land, nc=ncol(MASK), side=res(MASK)[1]/10^3, 
                                           radius=radius, outbreak, preoutbreak)
          ## Level of defoliation of the cells recently integrated in the outbreak (the sbw.new.spread cells)
          ## It can be 0 (no-defol), 1, 2 or 3!
          land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
            sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.1,0.8,0.1))
        }
        else if(calm>0 | collapse==1){ 
          ## not add new locations to the current outbreak if it is collapsing
          ## when calm, few spontaneously cells will have to be here and there defoliated
          potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"), temp>0.5, temp<2.8)
          sbw.new.sprd = sample(potential$cell.id, size=round(rlnorm(1, 2,1.5)), replace=F,  #mn(lnorm) ~ ifelse(collapse==1, 2, 1.5)
                                 prob=potential$ny.def0*(1200-potential$elev)/100)
          land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
            sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
        }
        
        ## Update SBW tracking variables for newly defoliated cells
        land$ny.def[land$cell.id %in% sbw.new.sprd] = 
          ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]==0, 0, 1)
        land$ny.def0[land$cell.id %in% sbw.new.sprd] = 
          ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]>0, 0, 
                 land$ny.def0[land$cell.id %in% sbw.new.sprd]+1)
        land$cum.intens.def[land$cell.id %in% sbw.new.sprd] = 
          land$curr.intens.def[land$cell.id %in% sbw.new.sprd]

        
        ## 2.2. Update level of defoliation intensity of already defoliated cells, that is, 
        ## cells with ny.def0<5 & tssbw>=30, and ny.def<=18 --> otherwise, the defoliation is forced to stop.
        sbw.curr.sprd = land$cell.id[land$ny.def0<5 & land$tssbw>=30 & land$ny.def<=18 & 
                                        land$cell.id %notin% sbw.new.sprd]
        curr.outbreak = filter(land, cell.id %in% sbw.curr.sprd)
        ## Level of defoliation of these cells. It can be 0 (no-defol), 1, 2 or 3!
        land$curr.intens.def[land$cell.id %in% sbw.curr.sprd] = 
          intens.def.curr(filter(land, cell.id %in% sbw.curr.sprd), niche.opt, niche.good, niche.poor, outbreak, collapse, calm)
        ## Update SBW tracking variables
        land$ny.def[land$cell.id %in% sbw.curr.sprd] = land$ny.def[land$cell.id %in% sbw.curr.sprd]+
          ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]==0, 0, 1)
        land$ny.def0[land$cell.id %in% sbw.curr.sprd] = 
          ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]>0, 0, 
                 land$ny.def0[land$cell.id %in% sbw.curr.sprd]+1)
        land$cum.intens.def[land$cell.id %in% sbw.curr.sprd] = 
          land$cum.intens.def[land$cell.id %in% sbw.curr.sprd]+land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]  

        ## 2.3. Stop defoilating cells that have been defoliated for more than 18 years
        ## Warning, not consider cells that have just been defoliated
        ## And avoid recurrence of defoliation in this oubreak by making tssbw==0
        sbw.stop = land$cell.id[land$ny.def0==0 & land$ny.def>18 & land$cell.id %notin% sbw.curr.sprd]
        land$curr.intens.def[land$cell.id %in% sbw.stop] = 0
        land$ny.def0[land$cell.id %in% sbw.stop] = 1
        land$tssbw[land$cell.id %in% sbw.stop] = 0
        
        ## 2.4. Increase number of year of non-defoliation for the remaing cells
        land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)] = 
          land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)]+1
        
        ## 2.5. For those cells that have just accumulated 5 years of non-defoliation, 
        ## reset the counter of number of years of defoliation and 
        ## cumulative intensity of defoliaion
        land$ny.def[land$ny.def0==5] = 0
        land$cum.intens.def[land$ny.def0==5] = 0
        
        ## 6. Kill cells according to number of years of defoliation and species composition
        ## with probability as cumulative*current intensity of defoliation
        if(outbreak>0 | collapse>0){
          kill.cells = sbw.mortality(land)
          aux = filter(land, cell.id %in% kill.cells) %>% group_by(spp, ny.def, curr.intens.def) %>% 
            summarise(area=length(spp)*km2.pixel) 
          names(aux) = c("spp", "ny.def", "curr.intens.def", "area")
          track.sbw.kill = rbind(track.sbw.kill, data.frame(run=irun, year=t+year.ini, aux))
          ## Mark the killed cells, and reset SBW variables. SBW won't spread to recently killed cells (tssbw<30)
          land$tssbw[land$cell.id %in% kill.cells] = 0
          land$ny.def[land$cell.id %in% kill.cells] = 0
          land$ny.def0[land$cell.id %in% kill.cells] = 0
          land$cum.intens.def[land$cell.id %in% kill.cells] = 0
          land$curr.intens.def[land$cell.id %in% kill.cells] = 0
        }
        
        
        ## Outbreak parameters
        if(preoutbreak>0 & done){ # pre-epidemia
          preoutbreak = preoutbreak-1
          phase = "preoutbreak"
          if(preoutbreak==0)
            duration.last.outbreak = outbreak = rdunif(1,-1,1)+duration.last.outbreak  # +6 Gray2008
          done = FALSE
        }
        else if(outbreak>0 & done){ # epidemia
          outbreak = outbreak-1
          phase = "outbreak"
          if(outbreak==0)
            collapse = rdunif(1,3,5)
          done = FALSE
        }
        else if(collapse>0 & done){  # collapse
          phase = "collapse"
          collapse = collapse-1
          if(collapse==0) #finishing the collapse
            calm = round(rnorm(1, 15, 1))
          done = FALSE
        }
        else if(calm>0 & done){
          phase = "calm"
          calm = calm-1    
          if(calm==0)
            preoutbreak = rdunif(1,3,4)
        }
        cat(phase, "\n")
        sbw.new.sprd = numeric()
        sbw.curr.sprd = numeric()
        sbw.stop = numeric()
        done = TRUE
      }
      
        
      

      


      

      ##################################### VEGETATION DYNAMICS #####################################
      ## First of all, save a vector containing initial forest composition for comparison at the end of the for loop
      initial.forest.comp = land$spp
      if(enable.succ){  
         
        ## Natural regeneration of forest after disturbance depends on the nature of the disturbance, 
        ## the age of the stand at the time the disturbance occurred, and the environmental suitability
        ## according to climate and soils. Compute it:
        suitab = suitability(land, temp.suitability, prec.suitability, soil.suitability, suboptimal)


        ## Regeneration after sbw outbreak
        if(length(kill.cells)>0){
          buffer = buffer.mig4(land, kill.cells, spp.colonize.persist)
          land$spp[land$cell.id %in% kill.cells] = forest.trans(land, kill.cells, post.sbw.reg, buffer, 
                     suitab, spp.colonize.persist, dtype="O", p.failure, age.seed, suboptimal, enfeuil)
        }


      
        ## Natural succession of tree spp at every 40 years starting at Tcomp = 70
        chg.comp.cells = filter(land, (age-age.matu) %in% seq(40,400,40) & tscomp>=70) %>% select(cell.id)
        if(length(unlist(chg.comp.cells))>0){
          buffer = buffer.mig4(land, unlist(chg.comp.cells), spp.colonize.persist)
          land$spp[land$cell.id %in% unlist(chg.comp.cells)] = 
            forest.trans(land, unlist(chg.comp.cells), forest.succ, buffer, suitab, 
                         spp.colonize.persist, dtype="S", p.failure, age.seed, suboptimal, enfeuil)          
        }


      }
      
      
      ## Now, for each cell that has changed composition (because of natural succession or regeneration of
      ## post-distrubance), re-initialize Tcomp
      land$tscomp[land$spp != initial.forest.comp] = 0
      land$age[land$cell.id %in% burnt.cells] = 0
      land$age[land$cell.id %in% kill.cells] = 0
      land$age[land$cell.id %in% cc.cells] = 0
      land$tspcut[land$cell.id %in% c(cc.cells, kill.cells, burnt.cells)] = 0 # -(land$age.matu/2)  # negative values!!!
      

      ## Finally, Aging: Increment Time Since Disturbance and Time Last Forest Composition change by time.step 
      land$age = land$age + time.step
      land$tscomp = land$tscomp + time.step
      land$tsfire = land$tsfire + time.step      
      land$tssbw = land$tssbw + time.step     
      land$tsccut = land$tsccut + time.step
      land$tspcut = land$tspcut + time.step
      
      
      
      ##################################### TRACKING AND SPATIAL OUTS #####################################
      ## Species distribution per fire zone
      track.spp.firezone = rbind(track.spp.firezone, data.frame(run=irun, year=t+year.ini, 
                                group_by(land, frz, spp) %>% summarize(area=length(cell.id)*km2.pixel)))
      ## Fuel type distribution per fire zone
      fuels = fuel.type(land, fuel.types.modif, NA, NA)
      aux = group_by(fuels, frz, type) %>% summarize(n=length(frz)) %>% 
             left_join(zone.size, by="frz") %>% mutate(pct=n/x) %>% select(-n, -x)
      track.fuel = rbind(track.fuel, data.frame(run=irun, year=t+year.ini, aux))
      ## Age classes distribution per species and management unit      
      land$age.class = cut(land$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
      track.spp.age.class = rbind(track.spp.age.class, data.frame(run=irun, year=t+year.ini, 
                              group_by(land, mgmt.unit, spp) %>% count(age.class)) %>%  mutate(area=n*km2.pixel))
      ## Suitability classes distribution per bioclim.domain   
      suitab = suitability(land, temp.suitability, prec.suitability, soil.suitability, suboptimal) 
      aux = left_join(suitab, select(land, cell.id, bioclim.domain), by="cell.id") %>%
              group_by(bioclim.domain, potential.spp) %>% summarize(poor=sum(suit.clim==0)*km2.pixel, 
              med=sum(suit.clim==0.5)*km2.pixel, good=sum(suit.clim==1)*km2.pixel) 
      track.suit.class = rbind(track.suit.class, data.frame(run=irun, year=t+year.ini, aux))
      aux = group_by(fuel.type(land, fuel.types.modif, NA, NA), frz) %>% summarize(x=mean(flam))
      rm(suitab); rm(aux)
      ## Track defoliation intensity
      if(sum(land$curr.intens.def)>0){
        ncell.def = c(NA, rep(sum(land$curr.intens.def>0),length(unique(land$curr.intens.def))-1))
        track.sbw.defol.intens = rbind(track.sbw.defol.intens, 
            data.frame(run=irun, year=t+year.ini, phase=phase, group_by(land, curr.intens.def) %>% 
            summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def)))
      }
      else
        track.sbw.defol.intens = rbind(track.sbw.defol.intens, 
            data.frame(run=irun, year=t+year.ini, phase=phase, curr.intens.def=0, ncell=nrow(land), pct=NA))
      
      

      ## If required, plot maps at each time step 
      if(write.maps){
        MAP = MASK
        cat("... writing output layers", "\n")
          # MAP[!is.na(MASK[])] = land$spp
          # writeRaster(MAP, paste0(out.path, "/lyr/spp_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
          # MAP[!is.na(MASK[])] = land$age
          # writeRaster(MAP, paste0(out.path, "/lyr/Age_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
        if(is.clearcut & t %in% cc.schedule){
          MAP[!is.na(MASK[])] = land$tsccut
          writeRaster(MAP, paste0(out.path, "/lyr/TSCcut_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)  
        }
        if(is.wildfires & t %in% fire.schedule){
          MAP[!is.na(MASK[])] = land$tsfire
          writeRaster(MAP, paste0(out.path, "/lyr/TSF_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
        }
        if(is.sbw & t %in% sbw.schedule){
          MAP[!is.na(MASK[])] = land$curr.intens.def
          writeRaster(MAP, paste0(out.path, "/lyr/curr.def_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
          # MAP[!is.na(MASK[])] = land$cum.intens.def
          # writeRaster(MAP, paste0(out.path, "/lyr/cum.def_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
          # MAP[!is.na(MASK[])] = land$ny.def
          # writeRaster(MAP, paste0(out.path, "/lyr/ny.def_r", irun, "t", t, ".tif"), format="GTiff", overwrite=T)
        }
        # save(land, file=paste0(out.path, "/lyr/land_r", irun, "t", t, ".rdata"))
      }
      
    # } # t
  # } # irun
  
  cat("... writing outputs", "\n")
  write.table(track.sbw.kill[-1,], paste0(out.path, "/SppKilled.txt"), quote=F, row.names=F, sep="\t")
  write.table(track.sbw.defol.intens[-1,], paste0(out.path, "/SBWdefoliation.txt"), quote=F, row.names=F, sep="\t")
  
  
    } # t
} 

