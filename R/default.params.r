default.params = function(){
  return(list(
  
  ## Time lenght in years of a model simulation
  ## 80 time steps of 1 years, it covers the period 2020-2100.
  year.ini = 2020,
  
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