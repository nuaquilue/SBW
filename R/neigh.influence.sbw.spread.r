neigh.influence.sbw.spread = function(land, nc=804, side=2, radius=12){
  
  # radius = 1
  # y = c(-radius:radius); y
  # z = c(y, ncol(MASK)+y, -(ncol(MASK)+y)); z[z!=0]
  # 
  # radius = 2
  # y = c(-radius:radius); y
  # z = c(y, ncol(MASK)+y, -(ncol(MASK)+y), 2*ncol(MASK)+y, -(2*ncol(MASK)+y)); z[z!=0]
  
  ## Neighbor position and distance to the focal cell in km
  x = c(-radius:radius)
  z = x; z
  w = abs(z)*side; w
  for(i in 1:radius){
    z = c(z, i*nc+x, -(i*nc+x)) 
    w = c(w, sqrt((abs(x)*side)^2+(i*side)^2), sqrt((abs(x)*side)^2+(i*side)^2) )
  }
  z = z[z!=0]
  w = w[w!=0]
  
  ## Look at current level of the defoliation in the neighborhood of the cells not 
  ## currently defoliated and that last mortality by oubreak is at least 30 years.
  potential = filter(land, ny.def0>=5, tssbw>=30)
  nslice = 10
  upper = round(nrow(potential)/nslice)
  ## First slice
  potential.slice = potential[1:upper,]
  neigh.curr.def = .compute.neigh.curr.def(land, potential.slice, z, w)
  neigh.host.pref = .compute.neigh.host.pref(land, potential.slice, z, w)
  ## Second to n-1 slice
  for(i in 1:(nslice-1)){
    # cat(i, "\n")
    potential.slice = potential[(i*upper+1):((i+1)*upper),]
    neigh.curr.def = rbind(neigh.curr.def, .compute.neigh.curr.def(land, potential.slice, z, w))
    neigh.host.pref = rbind(neigh.host.pref, .compute.neigh.host.pref(land, potential.slice, z, w))
  }
  ## Last slice
  # cat("last")
  potential.slice = potential[((nslice-1)*upper+1):nrow(potential),]
  neigh.curr.def = rbind(neigh.curr.def, .compute.neigh.curr.def(land, potential.slice, z, w))
  neigh.host.pref = rbind(neigh.host.pref, .compute.neigh.host.pref(land, potential.slice, z, w))
  
  ## Aggregate all the info
  dta = data.frame(neigh.curr.def, neigh.host.pref$x)
  names(dta) = c("cell.id", "neigh.curr.def", "neigh.host.pref")
  
  return(dta)
}


.compute.neigh.curr.def = function(land, potential.slice, z, w){
  nneigh = length(z)
  neighs = data.frame(cell.id=rep(potential.slice$cell.id, each=nneigh),
                       neigh.id=rep(potential.slice$cell.id, each=nneigh) + rep(z, nrow(potential.slice)),
                       w=rep(w, nrow(potential.slice))) %>% 
    filter(neigh.id %in% land$cell.id)
  neigh.def = left_join(neighs, select(land, cell.id, curr.intens.def), by=c("neigh.id"="cell.id")) %>% 
    group_by(cell.id) %>% summarise(x=sum(curr.intens.def/(w*0.5)))
  return(neigh.def)
}

.compute.neigh.host.pref = function(land, potential.slice, z, w){
  nneigh = length(z)
  breaks = c(0,20,40,60,80,100,999)
  tags = c("C10", "C30", "C50", "C70", "C90", "OLD")
  neighs = data.frame(cell.id=rep(potential.slice$cell.id, each=nneigh),
                       neigh.id=rep(potential.slice$cell.id, each=nneigh) + rep(z, nrow(potential.slice)),
                       w=rep(w, nrow(potential.slice))) %>% 
    filter(neigh.id %in% land$cell.id)
  neigh.host = left_join(neighs, select(land, cell.id, spp, age), by=c("neigh.id"="cell.id")) %>% 
    mutate(age.class=cut(age, breaks=breaks, include.lowest=TRUE, right=TRUE, labels=tags)) %>% 
    mutate(host.pref=ifelse(spp=="SAB" & age.class %in% c("C50", "C70", "C90", "OLD"), 1,
                            ifelse(spp=="SAB" & age.class %in% c("C10", "C30"), 0.75,
                                   ifelse(spp=="EPN" & age.class %in% c("C50", "C70", "C90", "OLD"), 0.5,
                                          ifelse(spp=="EPN" & age.class %in% c("C10", "C30"), 0.25, 0)))) ) %>% 
    group_by(cell.id) %>% summarise(x=sum(host.pref/(w*0.5)))
  return(neigh.host)
}
