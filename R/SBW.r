#' SBW: A package for simulating spruce budworm outbreaks in Quebec province
#' 
#' The SBW package provides a set of functions to simulate four different phases of 
#' spruce budworm outbreaks: pre-epidemic, epidemic, collapse, and calm, as well as the 
#' landscape-scale processes associated to these (e.g. neighbour contagion, long-term dispersal).
#' It also includes the effects on forest communities (i.e. mortality induced by SBW recurrent defoliation).
#' 
#' @section Author(s):
#' \bold{Maintainer}: Núria Aquilué \email{nuria.aquilue@ctfc.cat}  
#' 
#' \bold{Authors}: Núria Aquilué, Mathieu Bouchard and Élise Filotas
#'
#' @docType package
#' @name SBW
#' 
#' @importFrom tidyr %>%
#' @importFrom dplyr group_by summarise filter select mutate count left_join
#' @importFrom RANN nn2
NULL
#> NULL