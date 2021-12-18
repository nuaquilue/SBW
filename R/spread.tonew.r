spread.tonew = function(land, nc, side, radius, outbreak, preoutbreak){
  
  # cat("SBW adjacent spreading", "\n" )
  
  ## Function to simulate spreading to cells not yet defoliated, that is, cells with ny.def0>=5 and tssbw>=30
  ## The future function sbw.spread.tonew(filter(land, ny.def0>=5)) will return a vector with cell.ids
  ## MB:  Prob.spread in cell c = Proportion host species in neighborhood x 
  ##      level of defoliation in neighborhood during previous year x some climatic variable 
  ##      (and perhaps wind direction in future versions)
  ## The question is, how many cells do I have to select from the pool of potential cells?
  ## It should be a number of cells proportional to the number of cells defoliated in the previous year¿?
  
  ## Compute neighborhood current defoliation and neighbrohood host preference
  potential = neigh.influence.sbw.spread(land, nc, side, radius)
  
  ## The probability of sbw spread into a cell is proportional to neigh.curr.def * neigh.host.pref
  ## Be careful that any of these two variables is [0,1], unless I rescale them.
  ## By now, the weight of these two factors is the same (0.5), but it can be a model parameter and 
  ## test the sensibility of the spreading according to it.
  w = 0.8
  potential$x = (w*scales::rescale(potential$neigh.curr.def, to=c(0,1)))*
    ((1-w)*scales::rescale(potential$neigh.host.pref, to=c(0,1)))
  
  ## Select only those cells that at least one neighbor is defoilated
  potential = filter(potential, neigh.curr.def>0, x>0)
  if(nrow(potential)>1){
    if(outbreak<=10)  ## consolidated part of the outbreak
      nnew = pmin(nrow(potential), sum(land$curr.intens.def>0)*runif(1, 0.15, 0.45))
    if(preoutbreak>0 | outbreak>10)  
      nnew = pmin(nrow(potential), sum(land$curr.intens.def>0)*runif(1, 0.75, 1))
    sbw.new.sprd = sample(potential$cell.id, round(nnew), replace=F, prob=potential$x)
  }
  else
    sbw.new.sprd = potential$cell.id
  
  return(sbw.new.sprd)
}