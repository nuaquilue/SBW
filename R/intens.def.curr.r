
intens.def.curr <- function(curr.outbreak, params, outbreak, collapse, calm){
  
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
                                      intensity.defoliation(active.defol, params))
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







