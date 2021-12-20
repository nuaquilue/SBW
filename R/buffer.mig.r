#' Buffer migration
#'
#' Determines the presence of species within a buffer around target cells.
#' 
#' @param land A \code{landscape} data frame with forest stand records in rows
#' @param target.cell.ids A vector of \code{cell.id} codes for those cells open for colonization, 
#' and for which the presence of source populations nearby must be assessed
#' @param tbls A list of default input tables as in \code{data(default.tables()} or 
#' a customized list of input tables
#'
#' @return A data frame with the \code{cell.id} of the target cells and the presence (TRUE) or 
#' absence (FALSE) of a sufficient number of source populations for each 
#' species: PET, BOJ, ERS, SAB, EPN, other, NonFor
#' 
#' @export
#' 
#' @examples
#' data(landscape)
#' buffer.mig(landscape, landscape$cell.id[runif(10,1,nrow(landscape))])
#' 

buffer.mig = function(land, target.cells.ids, tbls){  
  
  ## Input table
  spp.colonize.persist = tbls$spp.colonize.persist
  
  ## Coordinates of the target cells
  target.cells.ids.xy = land[land$cell.id %in% target.cells.ids, c("cell.id", "x", "y")]
  
  ## Radius ~ max. colonization distance per species, and min number of soruces to enable colonization
  radius.buff = spp.colonize.persist$rad
  nb.buff = spp.colonize.persist$nneigh
    
  # Source cells (i.e. potential colonizers) per species. Minimal age of 50 years.
  micro.boj = land[land$spp=="BOJ" & land$age>=50 & land$tscomp>=50,]
  micro.pet = land[land$spp=="PET" & land$age>=50 & land$tscomp>=50,]
  micro.ers = land[land$spp=="ERS" & land$age>=50 & land$tscomp>=50,]
  micro.epn = land[land$spp=="EPN" & land$age>=50 & land$tscomp>=50,]
  micro.sab = land[land$spp=="SAB" & land$age>=50 & land$tscomp>=50,]
  
  ### Calculate number of source populations in the neighbohood of each target cell. Colonization distances
  ### are species-specific.
  # PET
  list.cell.buff = nn2(micro.pet[,c("x","y")], target.cells.ids.xy[,c("x","y")], 
                        k=nb.buff[4] , searchtype='priority')
  nn.dists.pet = list.cell.buff$nn.dists[, nb.buff[4]] < radius.buff[4]
  
  # BOJ   
  list.cell.buff = nn2(micro.boj[,c("x","y")], target.cells.ids.xy[,c("x","y")], 
                        k=nb.buff[1], searchtype='priority')
  nn.dists.boj = list.cell.buff$nn.dists[,nb.buff[1]]< radius.buff[1]

  # ERS
  list.cell.buff = nn2(micro.ers[,c("x","y")], target.cells.ids.xy[,c("x","y")], 
                        k=nb.buff[3], searchtype='priority')
  nn.dists.ers = list.cell.buff$nn.dists[,nb.buff[3]]< radius.buff[3]
  
  # SAB   
  list.cell.buff = nn2(micro.sab[,c("x","y")], target.cells.ids.xy[,c("x","y")], 
                        k=nb.buff[5],  searchtype='priority')
  nn.dists.sab = list.cell.buff$nn.dists[,nb.buff[5]]< radius.buff[5]
  
  # EPN   
  list.cell.buff = nn2(micro.epn[,c("x","y")], target.cells.ids.xy[,c("x","y")], 
                        k=nb.buff[2],  searchtype='priority')
  nn.dists.epn = list.cell.buff$nn.dists[,nb.buff[2]]< radius.buff[2]
  
  # Build a data frame with the presence or absence of a sufficient number of
  # source populations of each species around each target cell. Currently set
  # at one in all scenarios, but could be modified.
  target.df = data.frame(target.cells.ids.xy, PET=nn.dists.pet, BOJ=nn.dists.boj, ERS=nn.dists.ers, 
                          SAB=nn.dists.sab, EPN=nn.dists.epn, OTH=TRUE, NonFor=TRUE)
  target.df = target.df[,-c(2,3)]  
  target.df = reshape::melt(target.df,id=c("cell.id"))
  names(target.df)[-1] = c("potential.spp", "press.buffer")
  target.df$potential.spp = as.character(target.df$potential.spp)
  return(target.df)

}

