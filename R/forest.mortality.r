forest.mortality <- function(land){
  
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