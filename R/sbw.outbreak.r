sbw.outbreak = function(custom.params = NA, rcp = NA, prec.proj = NA, temp.proj = NA,  
                        time.horizon = 80, nrun = 1, save.land = FALSE, out.seq = NA, out.path = NA){
  
  options(dplyr.summarise.inform=F)
  select = dplyr::select
  
  ## Load required functions (no need in a R-package)
  source("R/buffer.mig.r") 
  source("R/default.params.r") 
  source("R/forest.mortality.r") 
  source("R/forest.transitions.r")
  source("R/intens.def.curr.r")  
  source("R/intensity.defoliation.r")  
  source("R/neigh.influence.sbw.spread.r")
  source("R/neighbour.spp.r")
  source("R/select.others.r")
  source("R/spread.tonew.r")
  source("R/suitability.r") 
  
  ## Load model data (no need in a R-package)
  load(file="data/mask.rda")      # Raster mask of the study area
  load(file="data/land.sbw.rda")  # Forest and sbw outbreak data
  load("data/temp.suitability.rda")
  load("data/prec.suitability.rda")
  load("data/soil.suitability.rda")
  load("data/spp.colonize.persist.rda")
  load("data/post.sbw.reg.rda")
  load("data/forest.succ.rda")
  
  ## Initializations and verifications  --------------------------------------------------------------------
  cat("Data preparation ...\n") 
  
  ## Function to select items not in a vector
  `%notin%` = Negate(`%in%`)
  
  ## Compute cell resolution in km2
  km2.pixel = raster::res(mask)[1] * raster::res(mask)[2] / 10^6

  ## Get the list of default parameters and update user-initialized parameters
  params = default.params()
  if(!is.na(custom.params)){
    # Check class of custom.params
    if((!inherits(customParams, "list"))) {
      stop("'custom.params' must be a named list")
    }
    ## Check that the names of the customized parameters are correct
    if(!all(names(custom.params) %in% names(params)))
      stop("Wrong custom parameters names")
    params <- custom.param
  }
  
  ## Set the directory for writing spatial outputs (if indicated) 
  if(save.land){
    if(is.na(out.seq)){
      out.seq = seq(1, time.horizon, 1) 
    } 
    else{
      if(!all(out.seq %in% time.seq)){
        warning("Not all time steps in the output sequence provided are simulated.", call.=F)
      }
    }
    if(is.na(out.path)) stop("Directory path to save outputs not provided")
  }
  
  ## If provided by the user, load temperature and precipitation predictions for the whole study area
  ## Or any other climatic variable included as factor of change in the model !!!!!!!!
  ## Check that all time steps are included and columns names is ok
  prec.chg = temp.chg = F
  if(!is.na(prec.proj)){
    # Check class of prec.proj
    if((!inherits(prec.proj, "data.frame"))) {
      stop("'prec.proj' must be a named data frame")
    }
    ## Check that the names of the customized parameters are correct
    if(colnames(prec.proj)[1] != "cell.id" | ncol(prec.proj) != (length(time.seq)+1) ) 
      stop("Format of the precipitation projections data frame is not correct.")
    prec.chg = T
  } 
  if(!is.na(temp.proj)){
    # Check class of prec.proj
    if((!inherits(temp.proj, "data.frame"))) {
      stop("'temp.proj' must be a named data frame")
    }
    if(colnames(temp.proj)[1] != "cell.id" | ncol(temp.proj) != (length(time.seq)+1) ) 
      stop("Format of the temperature projections data frame is not correct.")
    temp.chg = T
  }
  ## Load precipitation and temperature projections provided with the package according to the climatic scenario.
  if(is.na(prec.proj) & !is.na(rcp)){
    # prec.proj = get(paste0("prec_", rcp))
    load(file=paste0("data/prec_", rcp, "_MIROC_ESM_CHEM.rdata"))  
    prec.proj = select(cc.prec, -x, -y)
    prec.chg = T
  }
  if(is.na(temp.proj) & !is.na(rcp)){   
    # temp.proj = get(paste0("temp_", rcp))
    load(file=paste0("data/temp_", rcp, "_MIROC_ESM_CHEM.rdata")) 
    temp.proj = select(cc.temp, -x, -y)
    temp.chg = T
  }
  
  ## Create the traking data frames
  track.sbw.defol.intens = NULL #data.frame(run=NA, year=NA, phase=NA, curr.intens.def=NA, ncell=NA, pct=NA)
  track.sbw.kill = NULL #data.frame(run=NA, year=NA, spp=NA, ny.def=NA, curr.intens.def=NA, area=NA)
  
  ## Build the time sequence
  time.seq = seq(params$time.step, time.horizon, params$time.step)
  
  
  ## Start the simulations --------------------------------------------------------------------
  cat("\n") 
  cat(paste("Simulations ...\n")) 
  irun = 1
  # for(irun in 1:nrun){
    
    ## Forest and sbw data.frame
    land = landscape
    
    ## Decide (top-down) duration of the current outbreak, epidimic phase and collapse
    outbreak = 12 - params$current.duration  # testing
    collapse = params$collapse
    calm = params$calm
    preoutbreak = params$preoutbreak
    outbreak = rdunif(1,-1,1)+6+params$duration.last.outbreak - params$current.duration
    duration.last.outbreak = outbreak + params$current.duration
    phase = "outbreak"
    done = T
    
    ## Record initial distributions:
    ncell.def = c(NA, rep(sum(land$curr.intens.def>0),3))
    out = data.frame(run=irun, year=params$year.ini, phase=phase, group_by(land, curr.intens.def) %>% 
                      summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def))
    if(is.null(track.sbw.defol.intens)){
      track.sbw.defol.intens = out
    }
    else{
      track.sbw.defol.intens = rbind(track.sbw.defol.intens, out)  
    }
    
    ## Start 
    t = 1
    for(t in time.seq){
      
      ## Track scenario, replicate and time step
      cat(paste0("Replicate ", irun, "/", nrun, ". Time step: ", params$year.ini+t, "/", params$year.ini+time.horizon), "\n")
      
      ## Reiniziale vector of killed cells every time step
      kill.cells = integer()
      
      ## Update climatic variables at each time step if climate change is activated
      ## Column 1 is cell.index, the following columns are temp (precip) in 
      ## 2020-2025, 2025-2030, 2030-2035, ... etc.
      ## The last column (temp95) then corresponds to the period 2095-2100
      ## The first time step (t=5) we start with climate 2020-2025
      if(temp.chg & t < time.horizon){
        cat("  Update temperature projections\n")
        aux = cc.temp[,c(1,3+which(time.seq==t))] 
        names(aux) = c("cell.id", "temp")
        land = select(land, -temp) %>% left_join(aux, by="cell.id")
      }
      if(prec.chg & t < time.horizon){
        cat("  Update precipitation projections\n")
        aux = cc.prec[,c(1,3+which(time.seq==t))]
        names(aux) = c("cell.id", "prec")
        land = select(land, -prec) %>% left_join(aux, by="cell.id")
       }


      cat("SBW outbreak: ", "\n")
      
      ## 1. Defliation in the different phases of the SBW outbreak
      cat("Defoliation in the ")
      if(preoutbreak>0){
        cat("pre-epidemic phase: ", "\n")
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
            spread.tonew(land, nc=ncol(MASK), side=res(MASK)[1]/10^3, 
                             radius=12, outbreak, preoutbreak))
        # and finally assign intensity to all of them (rewrite intensity just assigned to epicenter cores)
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:3, size=length(sbw.new.sprd), replace=T, prob=c(0.2,0.4,0.3,0.1))
      }
      else if(outbreak>0){
        cat("epidemic phase: ", "\n")
        ## 1. Spatial spreading of the current outbreak to cells not yet defoliated, that is,
        ## cells with ny.def0>=5 & tssbw>=30
        ## The function 'spread.tonew' returns cell.ids.
        ## once in a while vary the radius, to allow spreading furhter away, or reduce it to limit outbreak....
        radius = rdunif(1,2,15) # 4 to 60 km
        sbw.new.sprd = spread.tonew(land, nc=ncol(mask), side=res(mask)[1]/10^3, 
                                         radius=radius, outbreak, preoutbreak)
        ## Level of defoliation of the cells recently integrated in the outbreak (the sbw.new.spread cells)
        ## It can be 0 (no-defol), 1, 2 or 3!
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.1,0.8,0.1))
        cat("hi.1", "\n")
      }
      else if(calm>0 | collapse==1){ 
        cat("calm phase: ", "\n")
        ## not add new locations to the current outbreak if it is collapsing
        ## when calm, few spontaneously cells will have to be here and there defoliated
        potential = filter(land, ny.def0>=5, tssbw>=30, spp %in% c("SAB", "EPN"), temp>0.5, temp<2.8)
        sbw.new.sprd = sample(potential$cell.id, size=round(rlnorm(1, 2,1.5)), replace=F,  #mn(lnorm) ~ ifelse(collapse==1, 2, 1.5)
                               prob=potential$ny.def0*(1200-potential$elev)/100)
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd] = 
          sample(0:2, size=length(sbw.new.sprd), replace=T, prob=c(0.5,0.3,0.2))
      }
        
      ## 2. Update SBW tracking variables for newly defoliated cells
      land$ny.def[land$cell.id %in% sbw.new.sprd] = 
        ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]==0, 0, 1)
      land$ny.def0[land$cell.id %in% sbw.new.sprd] = 
        ifelse(land$curr.intens.def[land$cell.id %in% sbw.new.sprd]>0, 0, 
               land$ny.def0[land$cell.id %in% sbw.new.sprd]+1)
      land$cum.intens.def[land$cell.id %in% sbw.new.sprd] = 
        land$curr.intens.def[land$cell.id %in% sbw.new.sprd]

      ## 3. Update level of defoliation intensity of already defoliated cells, that is, 
      ## cells with ny.def0<5 & tssbw>=30, and ny.def<=18 --> otherwise, the defoliation is forced to stop.
      sbw.curr.sprd = land$cell.id[land$ny.def0<5 & land$tssbw>=30 & land$ny.def<=18 & 
                                      land$cell.id %notin% sbw.new.sprd]
      curr.outbreak = filter(land, cell.id %in% sbw.curr.sprd)
      ## Level of defoliation of these cells. It can be 0 (no-defol), 1, 2 or 3!
      land$curr.intens.def[land$cell.id %in% sbw.curr.sprd] = 
        intens.def.curr(filter(land, cell.id %in% sbw.curr.sprd), params, preoutbreak, outbreak, collapse, calm)
      ## Update SBW tracking variables
      land$ny.def[land$cell.id %in% sbw.curr.sprd] = land$ny.def[land$cell.id %in% sbw.curr.sprd]+
        ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]==0, 0, 1)
      land$ny.def0[land$cell.id %in% sbw.curr.sprd] = 
        ifelse(land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]>0, 0, 
               land$ny.def0[land$cell.id %in% sbw.curr.sprd]+1)
      land$cum.intens.def[land$cell.id %in% sbw.curr.sprd] = 
        land$cum.intens.def[land$cell.id %in% sbw.curr.sprd]+land$curr.intens.def[land$cell.id %in% sbw.curr.sprd]  

      ## 4. Stop defoilating cells that have been defoliated for more than 18 years
      ## Warning, not consider cells that have just been defoliated
      ## And avoid recurrence of defoliation in this oubreak by making tssbw==0
      sbw.stop = land$cell.id[land$ny.def0==0 & land$ny.def>18 & land$cell.id %notin% sbw.curr.sprd]
      land$curr.intens.def[land$cell.id %in% sbw.stop] = 0
      land$ny.def0[land$cell.id %in% sbw.stop] = 1
      land$tssbw[land$cell.id %in% sbw.stop] = 0
        
      ## 5. Increase number of year of non-defoliation for the remaing cells
      land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)] = 
        land$ny.def0[land$cell.id %notin% c(sbw.new.sprd, sbw.curr.sprd, sbw.stop)]+1
        
      ## 6. For those cells that have just accumulated 5 years of non-defoliation, 
      ## reset the counter of number of years of defoliation and 
      ## cumulative intensity of defoliaion
      land$ny.def[land$ny.def0==5] = 0
      land$cum.intens.def[land$ny.def0==5] = 0
        
      ## 7. Kill cells according to number of years of defoliation and species composition
      ## with probability as cumulative*current intensity of defoliation
      if(outbreak>0 | collapse>0){
        kill.cells = forest.mortality(land)
        aux = filter(land, cell.id %in% kill.cells) %>% group_by(spp, ny.def, curr.intens.def) %>% 
          summarise(area=length(spp)*km2.pixel) 
        names(aux) = c("spp", "ny.def", "curr.intens.def", "area")
        track.sbw.kill = rbind(track.sbw.kill, data.frame(run=irun, year=t+params$year.ini, aux))
        ## Mark the killed cells, and reset SBW variables. SBW won't spread to recently killed cells (tssbw<30)
        land$tssbw[land$cell.id %in% kill.cells] = 0
        land$ny.def[land$cell.id %in% kill.cells] = 0
        land$ny.def0[land$cell.id %in% kill.cells] = 0
        land$cum.intens.def[land$cell.id %in% kill.cells] = 0
        land$curr.intens.def[land$cell.id %in% kill.cells] = 0
      }
        
      ## 8. Set outbreak parameters
      if(preoutbreak>0 & done){ # pre-epidemic
        preoutbreak = preoutbreak-1
        phase = "preoutbreak"
        if(preoutbreak==0)
          duration.last.outbreak = outbreak = rdunif(1,-1,1)+params$duration.last.outbreak  # +6 Gray2008
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
      # cat(phase, "\n")
      sbw.new.sprd = numeric()
      sbw.curr.sprd = numeric()
      sbw.stop = numeric()
      done = TRUE
      
      ## 9. Track defoliation intensity
      if(sum(land$curr.intens.def)>0){
        ncell.def = c(NA, rep(sum(land$curr.intens.def>0),length(unique(land$curr.intens.def))-1))
        track.sbw.defol.intens = rbind(track.sbw.defol.intens, 
        data.frame(run=irun, year=t+params$year.ini, phase=phase, group_by(land, curr.intens.def) %>% 
                   summarize(ncell=length(cell.id)) %>% mutate(pct=ncell/ncell.def)))
      }
      else
        track.sbw.defol.intens = rbind(track.sbw.defol.intens, 
        data.frame(run=irun, year=t+params$year.ini, phase=phase, curr.intens.def=0, ncell=nrow(land), pct=NA))
      
      

      ## Natural regeneration of forest after disturbance depends on the nature of the disturbance, 
      ## the age of the stand at the time the disturbance occurred, and the environmental suitability
      ## according to climate and soils.   
      if(params$enable.succ){  
        
        cat("Forest regeneration and succesion: ", "\n")
        
        ## First, save a vector containing initial forest composition for comparison
        initial.forest.comp = land$spp
        
        ## Environmental suitability:
        suitab = suitability(land, temp.suitability, prec.suitability, soil.suitability, params$suboptimal)

        ## Regeneration after sbw outbreak
        if(length(kill.cells)>0){
          buffer = buffer.mig(land, kill.cells, spp.colonize.persist)
          land$spp[land$cell.id %in% kill.cells] = forest.trans(land, kill.cells, post.sbw.reg, buffer, 
                     suitab, spp.colonize.persist, dtype="O", params$p.failure, params$age.seed, params$suboptimal, params$enfeuil)
        }

        ## Natural succession of tree spp at every 40 years starting at Tcomp = 70
        chg.comp.cells = filter(land, (age-age.matu) %in% seq(40,400,40) & tscomp>=70) %>% select(cell.id)
        if(length(unlist(chg.comp.cells))>0){
          buffer = buffer.mig(land, unlist(chg.comp.cells), spp.colonize.persist)
          land$spp[land$cell.id %in% unlist(chg.comp.cells)] = 
            forest.trans(land, unlist(chg.comp.cells), forest.succ, buffer, suitab, 
                         spp.colonize.persist, dtype="S", params$p.failure, params$age.seed, params$suboptimal, params$enfeuil)
        }
        
        ## For each cell that has changed composition (because of natural succession or regeneration of
        ## post-distrubance), re-initialize Tcomp and age
        land$tscomp[land$spp != initial.forest.comp] = 0
      }
      
      ## Reset to 0 age of killed cells and do aging 
      ## (i.e. increment of age and time since sbw and forest composition change)
      land$age[land$cell.id %in% kill.cells] = 0
      land$age = land$age + params$time.step
      land$tscomp = land$tscomp + params$time.step
      land$tssbw = land$tssbw + params$time.step     

      ## If required, save landscape data frame at each time step 
      if(save.land & t %in% out.seq){
        if(!file.exists(out.path))
          dir.create(file.path(out.path), showWarnings = T) 
        saveRDS(land, file=paste0(out.path, "landscape_", irun, "t", t, ".rds"))
      }

    } # t
  #} #irun    

  res = list(sbw.defol.intens = track.sbw.defol.intens, sbw.kill = track.sbw.kill)
  return(res)  
      
} 

