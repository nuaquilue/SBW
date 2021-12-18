#' Default model parameters
#'
#' Returns a list with the Quebec Landscape Dynamic Model parameters and the default values
#'
#' @return A list with the following items:
#'  \itemize{
#'    \item{\code{year.ini}: The inital year YYYY of the simulations}
#'    \item{\code{time.step}: Time step in years (default value is 1)}
#'    \item{\code{duration.last.outbreak}: String to indicate the unit of timber supply calculation}
#'    \item{\code{current.duration}: A flag to indicate the systematic forestation of conifer communities that become broadleaves communities following a perturbation}
#'    \item{\code{collapse}: Number of years of the current collapse phase. If 0 the phase is not activated yet}
#'    \item{\code{calm}: Number of years of the current calm phase. If 0 the phase is not activated yet}
#'    \item{\code{preoutbreak}: Number of years of the current preoutbreak phase. If 0 the phase is not activated yet}
#'    \item{\code{niche.opt}: Weight of the optimal niche for sbw}
#'    \item{\code{niche.good}: Weight of the good niche for sbw}
#'    \item{\code{niche.poor}: Weight of the poor niche for sbw}
#'    \item{\code{enable.succ}: A flag to enable natural succession every 40 years (if FLASE, composition remains the same)}
#'    \item{\code{enfeuil}: After clear-cut, for a percentage of conifer species (EPN and SAB), transform to pionner species (PET)}
#'    \item{\code{age.seed}: Below this stand age, seed production is very low, and regeneration failures are more likely}
#'    \item{\code{p.failure}: Probability of regeneration failure in young (< 50 years) burned stands ([0,1])}
#'    \item{\code{suboptimal}: Tolerance for sub optimal conditions ([0,1])}
#'  }
#'    
#' @export
#' 
#' @examples
#' # Change the replanif parameter to not re-calculate sustainable timber supply levels every time step
#' params = default.params()
#' params$replanif = FALSE
#' 

default.params = function(){
  return(list(
  
  ## Time lenght in years of a model simulation
  ## 80 time steps of 1 years, it covers the period 2020-2100.
  year.ini = 2020,
  time.step = 1,
  
  ## SPRUCE BUDWORM parameters:  
  duration.last.outbreak = 9,
  current.duration = 10, # from 2011 to 2020, but 14 if from 2007 to 2020
  collapse = 0,
  calm = 0,
  preoutbreak = 0,
  niche.opt = 1,
  niche.good = 0.6,
  niche.poor = 0.3,
  
  ## VEGETATION DYNAMICS parameters:
  enable.succ = TRUE, # enable natural succession every 40 years (if FLASE, composition remains the same)
  enfeuil = 0.0,
  age.seed = 40,     # below this stand age, seed production is very low, and regeneration failures are more likely
  p.failure = 0,     # probability of regeneration failure in young (< 50 years) burned stands
  suboptimal = 0.5  # tolerance for sub optimal conditions
  
  ))
  
}