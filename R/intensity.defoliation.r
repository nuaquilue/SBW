intensity.defoliation = function(outbreak, params, tbl){

  ## Host predisposition according to species and age class
  breaks = c(0,30,60,999)
  tags = c("C30", "C60", "OLD")
  outbreak$age.class = cut(outbreak$age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)
  outbreak$host.pref = ifelse(outbreak$spp=="SAB" & outbreak$age.class %in% c("C60", "OLD"), 1,
                               ifelse(outbreak$spp=="SAB" & outbreak$age.class=="C30", 0.75,
                                      ifelse(outbreak$spp=="EPN" & outbreak$age.class %in% c("C60", "OLD"), 0.5,
                                             ifelse(outbreak$spp=="EPN" & outbreak$age.class=="C30", 0.25, 0))))
  
  
  ## Soil and climatic suitability per tree spp. 
  ## Overall suitability corresponds to the minimum value between soil and climate suitability
  suboptimal = 0.5
  list.spp = levels(outbreak$spp)
  list.spp = list.spp[list.spp!="NonFor"]
  temp.suitability = tbl$temp.suitability
  prec.suitability = tbl$prec.suitability
  soil.suitability = tbl$soil.suitability
  dta = data.frame(cell.id=NA, spp=NA, site.qual=NA)
  for(ispp in list.spp){
    th.temp = filter(temp.suitability, spp==ispp)[-1]
    th.prec = filter(prec.suitability, spp==ispp)[-1] 
    th.soil = filter(soil.suitability, spp==ispp)[-1]
    aux = filter(outbreak, spp==ispp) %>%  #data.frrame(cell.id=vuln$cell.id, potential.spp=ispp, temp=vuln$temp, prec=vuln$prec, soil=vuln$soil.type) %>%
      mutate(class.temp=as.numeric(cut(temp, th.temp)),
             class.prec=as.numeric(cut(prec, th.prec)),
             suit.temp=ifelse(is.na(class.temp), 0, ifelse(class.temp==2, 1, suboptimal)),
             suit.prec=ifelse(is.na(class.prec), 0, ifelse(class.prec==2, 1, suboptimal)),
             # suit.soil=as.numeric(th.soil[match(soil.type, c("T","O","R","S","A"))]),
             suit.soil=1,
             suit.clim=pmin(suit.temp, suit.prec),
             site.qual=pmin(suit.soil, suit.clim))  %>%
      select(cell.id, spp, site.qual)
    dta = rbind(dta, aux)
  }
  dta = dta[-1,]
  aux = filter(outbreak, spp=="NonFor") %>% mutate(site.qual=0) %>% select(cell.id, spp, site.qual)
  dta = rbind(dta, aux) 
  outbreak = left_join(outbreak, dta, by=c("cell.id", "spp"))
  
  
  ## Time to include the effect of SBW climatic niche in determining defoliation intensity
  ## The reduction of the Likelihood of severe defoliation  according to sbw.niche should be a model parameter
  ## By now I use the value MB used to calculate probability of mortality, but's totally arbitrary.
  outbreak$sbw.niche = ifelse(outbreak$temp>0.5 & outbreak$temp<2.8, params$niche.opt,
                               ifelse(outbreak$temp>-1.5 & outbreak$temp<4, params$niche.good, params$niche.poor))
  
  
  ## Level of defoliation will be higher as host.pref is higher, and higher as site.qual is lower
  ## Weight the role of site.quality in determining defoliation intensity with 'wsite'
  ## Likelihood of severe defoliation decreases when sbw.niche is not suitable 
  wsite = 0.2
  outbreak$x = (outbreak$host.pref + wsite*(1-dta$site.qual))/(1+wsite)
  outbreak$x = outbreak$x * outbreak$sbw.niche
  
  
  ## At the end, defoliation intensity can be 0, 1, 2, or 3
  ## 0 means no defoliation, and will be mostly where sbw.niche is NOT suitable
  levels = c("no","low","moderate","severe")
  a = cut(outbreak$x, breaks=c(0,0.005,0.25,0.5,1), include.lowest=TRUE, right=TRUE, labels=levels)
  
  return(as.numeric(a)-1)
}