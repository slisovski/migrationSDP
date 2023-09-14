# library(fGarch)
# library(scales)
# library(sf)
# library(lwgeom)


## S2 Class
setClass(
  "migrationSDP",
  slots = c(Init    = "list",
            Species = "list",
            Sites   = "list",
            Results = "list")
)

make_migrationSDP <- function(init, species, sites, parms, ...) {
  
  crds_sf <- sites$crds %>% st_as_sf(coords = names(crds), crs = 4326)
  
  distM   <- st_distance(crds_sf, by_element = F)/1000
  
  if(parms$bearing) {

    bearM <- suppressWarnings(distm(crds_sf %>% st_shift_longitude() %>% st_coordinates(), fun = geosphere::bearingRhumb))
    bearM[upper.tri(bearM)] <- (bearM[upper.tri(bearM)] + 180) %% 360
    
    # bearM <- ifelse(bearM<(180-parms$angle) | bearM>(180+parms$angle), 0, 1)
    bearM <- ifelse(bearM>70 & bearM<290, 1, 0)
    
  } else bearM <- matrix(1, ncol = nrow(crds_sf), nrow = nrow(crds_sf))
  
  penalty <- 1-rep(sites$penalty, c(1, nrow(crds)-2, 1))
  
  new(
    "migrationSDP",
    Init  = list(
      MaxT   = length(init$minT:init$maxT)-1,
      NSites = nrow(sites$crds)-1,
      MaxX   = init$MaxX
    ),
    Species = list(
      B0    = species$B0,
      w     = species$w,
      xc    = species$xc,
      c     = species$c,
      speed = species$speed,
      WindAssist = 0,
      WindProb   = 1,
      ZStdNorm   = c(-2.5, -2.0, -1.5, -1.0, -0.5,  0.0,  0.5,  1.0,  1.5,  2.0,  2.5),
      PStdNorm   = c(0.0092, 0.0279, 0.0655, 0.1210, 0.1747, 0.2034, 0.1747, 0.1210, 0.0655, 0.0279, 0.0092),
      xFTReward  = c(init$minT:init$maxT)-init$minT,
      yFTReward  = scales::rescale(do.call(init$tReward[[1]], c(x = list(init$minT:init$maxT), unlist(init$tReward[-c(1, length(init$tReward))]))), 
                                   c(0, init$tReward$factor)),
      decError   = init$dError
    ),
    Sites = list(
      crds  = sites$crds,
      dist  = distM,
      bear  = bearM,
      angle = parms$angle,
      b0    = rep(sites$pred[1], nrow(sites$crds)),
      b1    = rep(sites$pred[2], nrow(sites$crds)),
      b2    = rep(sites$pred[3], nrow(sites$crds)),
      pred_a1 =  sites$pred[4],
      pred_a2 =  sites$pred[4],
      gain    = sites$gain,
      expend  = sites$expend,
      penalty = penalty
    ),
    Results = list(
      FitnessMatrix     = NA,
      DecisionMatrix    = NA,
      ProbMatrix        = NA
    )
  )
  
} 

bwdIteration <- function(obj) {

  Init(1,
       obj@Init$MaxT,          ## integer
       obj@Init$NSites,        ## integer
       obj@Init$MaxX,          ## double
       obj@Species$w,          ## double
       obj@Species$xc,         ## double 
       obj@Species$B0,         ## int
       obj@Sites$b0,           ## Rcpp::NumericVector
       obj@Sites$b1,           ## Rcpp::NumericVector
       obj@Sites$b2,           ## Rcpp::NumericVector
       obj@Sites$pred_a1,      ## double
       obj@Sites$pred_a2,      ## double
       obj@Species$c,          ## double
       obj@Species$speed,      ## double
       obj@Species$WindAssist, ## Rcpp::NumericVector
       obj@Species$WindProb,   ## Rcpp::NumericVector
       obj@Species$ZStdNorm,   ## Rcpp::NumericVector  
       obj@Species$PStdNorm,   ## Rcpp::NumericVector
       obj@Species$xFTReward,  ## Rcpp::NumericVector
       obj@Species$yFTReward,  ## Rcpp::NumericVector
       obj@Species$decError,   ## double
       obj@Sites$dist,         ## arma::mat
       obj@Sites$bear,         ## arma::mat
       obj@Sites$angle,        ## double
       obj@Sites$gain,         ## Rcpp::NumericVector 
       obj@Sites$expend,       ## arma::mat
       obj@Sites$penalty       ## Rcpp::NumericVector 
       )
  
  out <- BackwardIteration()

  obj@Results$FitnessMatrix <- out[[1]]
  DM <- array(dim = c(dim(out[[2]]),2))
  DM[,,,1] <- out[[2]]
  DM[,,,2] <- out[[3]]
  obj@Results$DecisionMatrix <- DM
  PM <- array(dim = c(dim(out[[4]]),2))
  PM[,,,1] <- out[[4]]
  PM[,,,2] <- out[[5]]
  obj@Results$ProbMatrix <- PM

  obj
}


#### Foreward
fwdSimulation <- function(model, NrInd, start_t, start_site, start_x) {
  
  InitSim(1, 
          model@Init$MaxT, 
          model@Init$NSites, 
          model@Init$MaxX,
          model@Species$w,
          model@Species$xc,
          model@Species$B0,
          model@Sites$b0,
          model@Sites$b1,
          model@Sites$b2,
          model@Sites$pred_a1,
          model@Sites$pred_a2,
          model@Species$c,
          model@Species$speed,
          model@Species$WindAssist,
          model@Species$WindProb,
          model@Species$ZStdNorm,
          model@Species$PStdNorm,
          model@Species$xFTReward,
          model@Species$yFTReward,
          model@Species$decError,
          model@Sites$dist,
          model@Sites$bear,
          model@Sites$gain,
          model@Sites$expend)
  
  
  x <- round(runif(NrInd, start_x[1], start_x[2]),0)

  if(length(start_site)>1 & length(start_site)<start_x[1]) {
    stop("start_site must have same length as numbers of individuals or a single site.")
  }
  if(length(start_site)==1) start_site <- rep(start_site, NrInd)

  SimOut = array(dim = c(length(x), 6, dim(model@Results$FitnessMatrix)[1]))

  ### First entry
  for(i in 1:dim(SimOut)[1]) {
    SimOut[i, ,start_t] <- c(start_t, start_site[i], x[i], 0, 0, 0)
  }


  ## SimOut: 1 = time, 2 = site, 3 = x, 4 = decision, 5 = flying, 6 = dead {
  for(time in 1:(dim(SimOut)[3]-1)) {

    for(ind in 1:dim(SimOut)[1]) {

      ## Not dead, not arrived, not flying
      if(!SimOut[ind, 6, time] &
         sum(SimOut[ind, 2, ] >= nrow(model@Sites$crds), na.rm = T)<1 & !SimOut[ind, 5, time]) {

        ## Decision
        if(runif(1) <  model@Results$ProbMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]) {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 1]
        } else {
          decision  <- model@Results$DecisionMatrix[SimOut[ind, 2, time], time, SimOut[ind, 3, time], 2]
        }

        ## Action
        if(decision>=0) { ## Flying

          fl_help = simFlying(decision, time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])

          nextt = fl_help[1] + 1
          if(nextt<=time) nextt <- time+1
          if(nextt>dim(SimOut)[3]) time <- dim(SimOut)[3]

          nextx = fl_help[2] + 1
          if(nextx <  0) {
            nextx = 0
            dead  = 1 } else  dead = 0
          if(nextx > model@Init$MaxX) nextx = model@Init$MaxX

          SimOut[ind,,nextt] = c(nextt, decision+1, nextx, NA, 0, dead)
          if(nextt>(time+1)) SimOut[ind,5:6,(time+1):(nextt-1)] = cbind(1,0)
          if(SimOut[ind, 6, nextt])  SimOut[ind, 6, nextt:dim(SimOut)[3]] = 1 ## if dead make dead till the end

          if(SimOut[ind, 2, nextt]==nrow(model@Sites$crds)) {
            SimOut[ind,2:6, nextt:dim(SimOut)[3]] <- SimOut[ind, 2:6, nextt]
            SimOut[ind,1,   nextt:dim(SimOut)[3]] <- seq(nextt, dim(SimOut)[3])
          }

        } else { ## Feeding

          fo_help = simForaging(abs(decision+1.0), time-1, SimOut[ind, 2, time]-1, SimOut[ind, 3, time])

          newx = fo_help[1]+1
          dead = fo_help[2]
          if(newx<=0) dead = 1

          if(newx > model@Init$MaxX) newx = model@Init$MaxX

          SimOut[ind,,time+1] = c(time+1, SimOut[ind, 2, time], newx, abs(decision+1.0), 0, dead)

          if(SimOut[ind, 6, time+1])  SimOut[ind, 6, (time+1):dim(SimOut)[3]] = 1 ## if dead make dead till the end

        }

      }

    } ## Ind loop
  } ## time loop

  SimOut       <- SimOut[,-5,]
  SimOut[,2,] <- SimOut[,2,]-1
  
  SimOut
}

