sbw.spread.tonew <- function(land, nc, side, radius, outbreak, preoutbreak){
  
  # cat("SBW adjacent spreading", "\n" )

  ## Function to simulate spreading to cells not yet defoliated, that is, cells with ny.def0>=5 and tssbw>=30
  ## The future function sbw.spread.tonew(filter(land, ny.def0>=5)) will return a vector with cell.ids
  ## MB:  Prob.spread in cell c <- Proportion host species in neighborhood x 
  ##      level of defoliation in neighborhood during previous year x some climatic variable 
  ##      (and perhaps wind direction in future versions)
  ## The question is, how many cells do I have to select from the pool of potential cells?
  ## It should be a number of cells proportional to the number of cells defoliated in the previous year¿?
  
  ## Compute neighborhood current defoliation and neighbrohood host preference
  potential <- neigh.influence.sbw.spread(land, nc, side, radius)
  
  ## The probability of sbw spread into a cell is proportional to neigh.curr.def * neigh.host.pref
  ## Be careful that any of these two variables is [0,1], unless I rescale them.
  ## By now, the weight of these two factors is the same (0.5), but it can be a model parameter and 
  ## test the sensibility of the spreading according to it.
  w <- 0.8
  potential$x <- (w*scales::rescale(potential$neigh.curr.def, to=c(0,1)))*
                 ((1-w)*scales::rescale(potential$neigh.host.pref, to=c(0,1)))
                                                        
  ## Select only those cells that at least one neighbor is defoilated
  potential <- filter(potential, neigh.curr.def>0, x>0)
  if(nrow(potential)>1){
    if(outbreak<=10)  ## consolidated part of the outbreak
      nnew <- pmin(nrow(potential), sum(land$curr.intens.def>0)*runif(1, 0.15, 0.45))
    if(preoutbreak>0 | outbreak>10)  
      nnew <- pmin(nrow(potential), sum(land$curr.intens.def>0)*runif(1, 0.75, 1))
    sbw.new.sprd <- sample(potential$cell.id, round(nnew), replace=F, prob=potential$x)
  }
  else
    sbw.new.sprd <- potential$cell.id
  
  return(sbw.new.sprd)
}


intens.def.curr <- function(curr.outbreak, niche.opt, niche.good, niche.poor, 
                            outbreak, collapse, calm){
  
  `%notin%` <- Negate(`%in%`)
  
  ## Intensity of defoliation for those cells that were not defoliated in the previous year
  no.defol <- filter(curr.outbreak, curr.intens.def==0)
  if(nrow(no.defol)>0){
    if(collapse>0 | calm>0)
      no.defol$new.intens <- 0
    else
      no.defol$new.intens <- ifelse(abs(rnorm(nrow(no.defol), 0, 2))<=no.defol$cum.intens.def*no.defol$ny.def0, 0, 
                                    sample(1:3, prob=c(1,0.25,0.05), replace=T))
    # table(no.defol$new.intens)
  }
  
  ## Intensity of defoliation for those cells that were defoliated in the previous year
  active.defol <- filter(curr.outbreak, curr.intens.def>0)
  if(collapse>1){
    to.collapse <- rdunif(1, 0.2*nrow(active.defol), 0.6*nrow(active.defol))
    to.collapse <- sample(active.defol$cell.id, size=to.collapse, replace=F, prob=(1/active.defol$ny.def)^2)
    active.defol$new.intens[active.defol$cell.id %in% to.collapse] <-  0
    active.defol$new.intens[active.defol$cell.id %notin% to.collapse] <- 
      pmax(active.defol$curr.intens.def[active.defol$cell.id %notin% to.collapse] -
           rdunif(nrow(active.defol)-length(to.collapse), 0, 1) , 0)
  }
  else if(collapse==1){
    active.defol$new.intens <- 0
  }
  else if(outbreak>0 | preoutbreak>0){  ## or calm==0 ¿?
    active.defol$new.intens <- ifelse(runif(nrow(active.defol),0,1)<0.75, active.defol$curr.intens.def,
        intensity.defoliation(active.defol, niche.opt, niche.good, niche.poor))
    active.defol$new.intens[active.defol$new.intens==1 & active.defol$ny.def>1] <- 
      sample(1:3, size=sum(active.defol$new.intens==1 & active.defol$ny.def>1), replace=T, prob=c(0.7,0.2,0.1))
  }
  else if(calm>1){
    to.collapse <- rdunif(1, 0.2*nrow(active.defol), 0.6*nrow(active.defol))
    to.collapse <- sample(active.defol$cell.id, size=to.collapse, replace=F)
    active.defol$new.intens[active.defol$cell.id %in% to.collapse] <-  0
    active.defol$new.intens[active.defol$cell.id %notin% to.collapse] <- 
      rdunif(nrow(active.defol)-length(to.collapse), 1, 3) 
  } 
  else if(calm==1){
    active.defol$new.intens <- 0
  }
  # table(active.defol$ny.def, active.defol$new.intens )
  # table(active.defol$curr.intens, active.defol$new.intens)
  if(nrow(no.defol)==0)
    kk <- active.defol
  else if(nrow(active.defol)==0)
    kk <- no.defol
  else
    kk <- rbind(no.defol, active.defol)
  kk <- left_join(select(curr.outbreak, cell.id), select(kk, cell.id, new.intens), by="cell.id")
  # sum(kk$new.intens>0)/length(kk$new.intens)
  # round(table(kk$new.intens)[-1]/sum(kk$new.intens>0),2)
  
  return(kk$new.intens)
  
}




intensity.defoliation <- function(outbreak, niche.opt, niche.good, niche.poor){
  
  options(warn=-1)  
  
  ## Host predisposition according to species and age class
  breaks <- c(0,30,60,999)
  tags <- c("C30", "C60", "OLD")
  outbreak$age.class <- cut(outbreak$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
  outbreak$host.pref <- ifelse(outbreak$spp=="SAB" & outbreak$age.class %in% c("C60", "OLD"), 1,
                              ifelse(outbreak$spp=="SAB" & outbreak$age.class=="C30", 0.75,
                                ifelse(outbreak$spp=="EPN" & outbreak$age.class %in% c("C60", "OLD"), 0.5,
                                  ifelse(outbreak$spp=="EPN" & outbreak$age.class=="C30", 0.25, 0))))
  
  
  ## Soil and climatic suitability per tree spp. 
  ## Overall suitability corresponds to the minimum value between soil and climate suitability
  suboptimal <- 0.5
  list.spp <- levels(outbreak$spp)
  list.spp <- list.spp[list.spp!="NonFor"]
  temp.suitability <- read.table("inputfiles/ThMeanTemp.txt", header=T)  
  prec.suitability <- read.table("inputfiles/ThAnnualPrecip.txt", header=T)  
  soil.suitability <- read.table("inputfiles/ThSoil.txt", header=T)  
  dta <- data.frame(cell.id=NA, spp=NA, site.qual=NA)
  for(ispp in list.spp){
    th.temp <- filter(temp.suitability, spp==ispp)[-1]
    th.prec <- filter(prec.suitability, spp==ispp)[-1] 
    th.soil <- filter(soil.suitability, spp==ispp)[-1]
    aux <- filter(outbreak, spp==ispp) %>%  #data.frrame(cell.id=vuln$cell.id, potential.spp=ispp, temp=vuln$temp, prec=vuln$prec, soil=vuln$soil.type) %>%
           mutate(class.temp=as.numeric(cut(temp, th.temp)),
             class.prec=as.numeric(cut(prec, th.prec)),
             suit.temp=ifelse(is.na(class.temp), 0, ifelse(class.temp==2, 1, suboptimal)),
             suit.prec=ifelse(is.na(class.prec), 0, ifelse(class.prec==2, 1, suboptimal)),
             # suit.soil=as.numeric(th.soil[match(soil.type, c("T","O","R","S","A"))]),
             suit.soil=1,
             suit.clim=pmin(suit.temp, suit.prec),
             site.qual=pmin(suit.soil, suit.clim))  %>%
          select(cell.id, spp, site.qual)
    dta <- rbind(dta, aux)
  }
  dta <- dta[-1,]
  aux <- filter(outbreak, spp=="NonFor") %>% mutate(site.qual=0) %>% select(cell.id, spp, site.qual)
  dta <- rbind(dta, aux) 
  outbreak <- left_join(outbreak, dta, by=c("cell.id", "spp"))
  
  
  ## Time to include the effect of SBW climatic niche in determining defoliation intensity
  ## The reduction of the Likelihood of severe defoliation  according to sbw.niche should be a model parameter
  ## By now I use the value MB used to calculate probability of mortality, but's totally arbitrary.
  outbreak$sbw.niche <- ifelse(outbreak$temp>0.5 & outbreak$temp<2.8, niche.opt,
                                   ifelse(outbreak$temp>-1.5 & outbreak$temp<4, niche.good, niche.poor))
  
  
  ## Level of defoliation will be higher as host.pref is higher, and higher as site.qual is lower
  ## Weight the role of site.quality in determining defoliation intensity with 'wsite'
  ## Likelihood of severe defoliation decreases when sbw.niche is not suitable 
  wsite <- 0.2
  outbreak$x <- (outbreak$host.pref + wsite*(1-dta$site.qual))/(1+wsite)
  outbreak$x <- outbreak$x * outbreak$sbw.niche
  
  
  ## At the end, defoliation intensity can be 0, 1, 2, or 3
  ## 0 means no defoliation, and will be mostly where sbw.niche is NOT suitable
  levels <- c("no","low","moderate","severe")
  a <- cut(outbreak$x, breaks=c(0,0.005,0.25,0.5,1), include.lowest=TRUE, right=TRUE, labels=levels)
  
  return(as.numeric(a)-1)
}


sbw.mortality <- function(land){
  
  ## Prob.mortality in cell c <- number of years of defoliation in cell c, and species composition in cell c   
  ## (exemple de paramètres: par ex for SAB: 5-7 years = prob 10% per year, 7-10 = 15%, >10 20%; 
  ## for EPN: 5-7 years = 0%, 7-10= 5 %, >10 = 10%)
  
  ## Vector to store cell.id of killed cells
  kill.cells <- numeric()
  
  ## Low mortality of Balsam fir
  target <- filter(land, spp=="SAB", ny.def>=5, ny.def<=7)
  # size = pmin(nrow(target)*0.1, sum(target$curr.intens.def>0))
  if(nrow(target)>1)
    aux <- sample(target$cell.id, size=nrow(target)*0.1, replace=F, 
                  prob=pmax(target$cum.intens.def*target$curr.intens.def, target$cum.intens.def))
  else
    aux <- target$cell.id
  kill.cells <- c(kill.cells,aux)
  
  ## Medium mortality of Balsam fir
  target <- filter(land, spp=="SAB", ny.def>=8, ny.def<=10)
  if(nrow(target)>1)
    aux <- sample(target$cell.id, size=nrow(target)*0.15, replace=F, 
                  prob=pmax(target$cum.intens.def*target$curr.intens.def, target$cum.intens.def))
  else
    aux <- target$cell.id
  kill.cells <- c(kill.cells,aux)
  
  ## High mortality of Balsam fir
  target <- filter(land, spp=="SAB", ny.def>=11)
  if(nrow(target)>1)
    aux <- sample(target$cell.id, size=nrow(target)*0.2, replace=F, 
                  prob=pmax(target$cum.intens.def*target$curr.intens.def, target$cum.intens.def))
  else
    aux <- target$cell.id
  kill.cells <- c(kill.cells,aux)
  
  ## Medium mortality of Black spruce
  target <- filter(land, spp=="EPN", ny.def>=8, ny.def<=10)
  if(nrow(target)>1)
    aux <- sample(target$cell.id, size=nrow(target)*0.05, replace=F, 
                  prob=pmax(target$cum.intens.def*target$curr.intens.def, target$cum.intens.def))
  else
    aux <- target$cell.id
  kill.cells <- c(kill.cells,aux)
  
  ## High mortality of Balsam fir
  target <- filter(land, spp=="EPN", ny.def>=11)
  if(nrow(target)>1)
    aux <- sample(target$cell.id, size=nrow(target)*0.1, replace=F, 
                  prob=pmax(target$cum.intens.def*target$curr.intens.def, target$cum.intens.def))
  else
    aux <- target$cell.id
  kill.cells <- c(kill.cells,aux)
  
  return(kill.cells)
  
}

