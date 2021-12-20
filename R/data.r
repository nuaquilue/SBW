#' Forest landscape features of Quebec province
#'
#' Forest stands characteristics of Quebec province in 2020 at 2 x 2 km of spatial resolution 
#'
#' @format A data frame with 147464 rows and 19 variables:
#' \describe{
#'   \item{cell.id}{Unique grid cell indentificator}
#'   \item{x}{X coordinate in Lambert Confromal Conic NAD83}
#'   \item{y}{Y coordinate in Lambert Confromal Conic NAD83}
#'   \item{temp}{Historical mean temperature, in ºC}
#'   \item{prec}{Mean annual precipitation, in mm}
#'   \item{frz}{Fire regime zone}
#'   \item{eco.type}{Ecological type}
#'   \item{bioclim.domain}{Bioclimatic domain}
#'   \item{mgmt.unit}{Management unit}
#'   \item{spp}{Tree species or group of species}
#'   \item{age.matu}{Forest stand maturity age}
#'   \item{soil.type}{Soil type}
#'   \item{exclus}{Excluded for timber harvesting}
#'   \item{tsfire}{Time since last fire, in years}
#'   \item{tssbw}{Time since last spruce budworm outbreak, in years}
#'   \item{tsccut}{Time since last clear cut, in years}
#'   \item{tspcut}{Time since last partical cut, in years}
#'   \item{tscomp}{Timse since last species change composition, in years}
#' }
#' @source \url{https://mffp.gouv.qc.ca/le-ministere/cartes-plans/}
"landscape"


#' Mask of the study area
#'
#' Binary raster to identify the study area (1 or NA) 
#'
#' @format Raster of 435 (nrow) x 804 (ncol)
#' \describe{
#' Raster of the study area (1 or NA) in the Lambert Conformal Conic NAD83 projection, at 2 x 2 km spatial resoltuion. Forest stands characteristics of Quebec province in 2020 at 2 x 2 km of spatial resolution.
#' The unique grid cell identificator \code{cell.id} in the \code{landscape} data frame coincides with the position of the location in the \code{mask} raster.
#' }
"mask"


#' Model input tables
#'
#' List of data frames being input tables of the SBW model
#' @format A list with 6 elements:
#' \describe{
#'    \item{\code{post.sbw.reg}: A data frame of species regeneration probability after mortality by SBW defoliation
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{potential.spp}: Code of the species to change.}
#'         \item{\code{ptrans}: Probability of transition.}
#'       }
#'    }
#'    \item{\code{forest.succ}: A data frame of transition probability of species seral succession
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{potential.spp}: Code of the species to change.}
#'         \item{\code{ptrans}: Probability of transition.}
#'       }
#'    }  
#'    \item{\code{spp.colonize.persist}: A data frame of species seral succession
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{rad}: Estimated maximum colonization distance (in m).}
#'         \item{\code{nneigh}: Minimum number of source cells within the colonization distance to enable colonization.}
#'         \item{\code{persist}:It indicates whether we allow the transition probability to remain high locally.}
#'       }
#'    }  
#'    \item{\code{temp.suitability}: A data frame with thresholds of temperature to define optimal and good species niche
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{th1}: Temperature lower threshold to define 'good' species niche (in ºC).}
#'         \item{\code{th2}: Temperature lower threshold to define 'optimal' species niche (in ºC).}
#'         \item{\code{th3}: Temperature upper threshold to define 'good' species niche (in ºC).}
#'         \item{\code{th4}: Temperature upper threshold to define 'optimal' species niche (in ºC).}
#'       }
#'    }  
#'    \item{\code{prec.suitability}: A data frame with thresholds of precipitation to define optimal and good species niche
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{th1}: Precipitation lower threshold to define 'good' species niche (in mm).}
#'         \item{\code{th2}: Precipitation lower threshold to define 'optimal' species niche(in mm).}
#'         \item{\code{th3}: Precipitation upper threshold to define 'good' species niche (in mm).}
#'         \item{\code{th4}: Precipitation upper threshold to define 'optimal' species niche (in mm).}
#'       }
#'    }  
#'    \item{\code{soil.suitability}: A data frame with thresholds of precipitation to define optimal and good species niche
#'      \itemize{
#'         \item{\code{spp}: Code of the species.}
#'         \item{\code{thT}: Species suitability in glacial tills [0,1].}
#'         \item{\code{thO}: Species suitability in organic deposists [0,1].}
#'         \item{\code{thR}: Species suitability in rock outcrops or very thin deposits [0,1].}
#'         \item{\code{thS}: Species suitability in sandy deposists [0,1].}
#'         \item{\code{thA}: Species suitability in clay deposits [0,1].}
#'       }
#'    }  
#'  }
"default.tables"

#' Annual precipitation projections under RCP4.5
#'
#' Projections of annual precipitation for Quebec province at 2 x 2 km of spatial resolution and 5-year temporal 
#' resolution under the RCP4.5 scenario derived from the MIROC_ESM_CHEM global model.
#'
#' @format A data frame with 185479 rows and 17 variables:
#' \describe{
#'   \item{cell.id}{Unique grid cell indentificator}
#'   \item{prec20}{Annual precipitation projection for 2020-2024, in mm}
#'   \item{prec25}{Annual precipitation projection for 2025-2029, in mm}
#'   \item{prec30}{Annual precipitation projection for 2030-2034, in mm}
#'   \item{prec35}{Annual precipitation projection for 2035-2039, in mm}
#'   \item{prec40}{Annual precipitation projection for 2040-2044, in mm}
#'   \item{prec45}{Annual precipitation projection for 2045-2049, in mm}
#'   \item{prec50}{Annual precipitation projection for 2050-2054, in mm}
#'   \item{prec55}{Annual precipitation projection for 2055-2059, in mm}
#'   \item{prec60}{Annual precipitation projection for 2060-2064, in mm}
#'   \item{prec65}{Annual precipitation projection for 2065-2069, in mm}
#'   \item{prec70}{Annual precipitation projection for 2070-2074, in mm}
#'   \item{prec75}{Annual precipitation projection for 2075-2079, in mm}
#'   \item{prec80}{Annual precipitation projection for 2080-2084, in mm}
#'   \item{prec85}{Annual precipitation projection for 2085-2089, in mm}
#'   \item{prec90}{Annual precipitation projection for 2090-2094, in mm}
#'   \item{prec95}{Annual precipitation projection for 2095-2099, in mm}
#' }
#' @source \url{https://www.ouranos.ca/en/}
"prec_rcp45"


#' Annual precipitation projections under RCP8.5
#'
#' Projections of annual precipitation for Quebec province at 2 x 2 km of spatial resolution and 5-year temporal 
#' resolution under the RCP8.5 scenario derived from the MIROC_ESM_CHEM global model.
#'
#' @format A data frame with 185479 rows and 17 variables:
#' \describe{
#'   \item{cell.id}{Unique grid cell indentificator}
#'   \item{prec20}{Annual precipitation projection for 2020-2024, in mm}
#'   \item{prec25}{Annual precipitation projection for 2025-2029, in mm}
#'   \item{prec30}{Annual precipitation projection for 2030-2034, in mm}
#'   \item{prec35}{Annual precipitation projection for 2035-2039, in mm}
#'   \item{prec40}{Annual precipitation projection for 2040-2044, in mm}
#'   \item{prec45}{Annual precipitation projection for 2045-2049, in mm}
#'   \item{prec50}{Annual precipitation projection for 2050-2054, in mm}
#'   \item{prec55}{Annual precipitation projection for 2055-2059, in mm}
#'   \item{prec60}{Annual precipitation projection for 2060-2064, in mm}
#'   \item{prec65}{Annual precipitation projection for 2065-2069, in mm}
#'   \item{prec70}{Annual precipitation projection for 2070-2074, in mm}
#'   \item{prec75}{Annual precipitation projection for 2075-2079, in mm}
#'   \item{prec80}{Annual precipitation projection for 2080-2084, in mm}
#'   \item{prec85}{Annual precipitation projection for 2085-2089, in mm}
#'   \item{prec90}{Annual precipitation projection for 2090-2094, in mm}
#'   \item{prec95}{Annual precipitation projection for 2095-2099, in mm}
#' }
#' @source \url{https://www.ouranos.ca/en/}
"prec_rcp85"


#' Mean temperature projections under RCP4.5
#'
#' Projections of mean temperature for Quebec province at 2 x 2 km of spatial resolution and 5-year temporal 
#' resolution under the RCP4.5 scenario derived from the MIROC_ESM_CHEM global model.
#'
#' @format A data frame with 185479 rows and 17 variables:
#' \describe{
#'   \item{cell.id}{Unique grid cell indentificator}
#'   \item{prec20}{Mean temperature projection for 2020-2024, in ºC}
#'   \item{prec25}{Mean temperature projection for 2025-2029, in ºC}
#'   \item{prec30}{Mean temperature projection for 2030-2034, in ºC}
#'   \item{prec35}{Mean temperature projection for 2035-2039, in ºC}
#'   \item{prec40}{Mean temperature projection for 2040-2044, in ºC}
#'   \item{prec45}{Mean temperature projection for 2045-2049, in ºC}
#'   \item{prec50}{Mean temperature projection for 2050-2054, in ºC}
#'   \item{prec55}{Mean temperature projection for 2055-2059, in ºC}
#'   \item{prec60}{Mean temperature projection for 2060-2064, in ºC}
#'   \item{prec65}{Mean temperature projection for 2065-2069, in ºC}
#'   \item{prec70}{Mean temperature projection for 2070-2074, in ºC}
#'   \item{prec75}{Mean temperature projection for 2075-2079, in ºC}
#'   \item{prec80}{Mean temperature projection for 2080-2084, in ºC}
#'   \item{prec85}{Mean temperature projection for 2085-2089, in ºC}
#'   \item{prec90}{Mean temperature projection for 2090-2094, in ºC}
#'   \item{prec95}{Mean temperature projection for 2095-2099, in ºC}
#' }
#' @source \url{https://www.ouranos.ca/en/}
"temp_rcp45"


#' Mean temperature projections under RCP8.5
#'
#' Projections of mean temperature for Quebec province at 2 x 2 km of spatial resolution and 5-year temporal 
#' resolution under the RCP8.5 scenario derived from the MIROC_ESM_CHEM global model.
#'
#' @format A data frame with 185479 rows and 17 variables:
#' \describe{
#'   \item{cell.id}{Unique grid cell indentificator}
#'   \item{prec20}{Mean temperature projection for 2020-2024, in ºC}
#'   \item{prec25}{Mean temperature projection for 2025-2029, in ºC}
#'   \item{prec30}{Mean temperature projection for 2030-2034, in ºC}
#'   \item{prec35}{Mean temperature projection for 2035-2039, in ºC}
#'   \item{prec40}{Mean temperature projection for 2040-2044, in ºC}
#'   \item{prec45}{Mean temperature projection for 2045-2049, in ºC}
#'   \item{prec50}{Mean temperature projection for 2050-2054, in ºC}
#'   \item{prec55}{Mean temperature projection for 2055-2059, in ºC}
#'   \item{prec60}{Mean temperature projection for 2060-2064, in ºC}
#'   \item{prec65}{Mean temperature projection for 2065-2069, in ºC}
#'   \item{prec70}{Mean temperature projection for 2070-2074, in ºC}
#'   \item{prec75}{Mean temperature projection for 2075-2079, in ºC}
#'   \item{prec80}{Mean temperature projection for 2080-2084, in ºC}
#'   \item{prec85}{Mean temperature projection for 2085-2089, in ºC}
#'   \item{prec90}{Mean temperature projection for 2090-2094, in ºC}
#'   \item{prec95}{Mean temperature projection for 2095-2099, in ºC}
#' }
#' @source \url{https://www.ouranos.ca/en/}
"temp_rcp85"

