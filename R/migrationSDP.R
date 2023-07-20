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
  
  # init =    list(minT    = ,
  #                maxT    = ,
  #                MaxX    = ,
  #                tReward = ),
  #                dError  = ),
  # species = list(B0      = ,
  #                w       = ,
  #                xc      = ,
  #                stdNorm = ,
  #                speed   = ),
  # sites   = list(crds    = ,
  #                expend  = ,
  #                gain    = ,
  #                penalty = ),
  # parms   = list(bearing = FALSE)
  
  crds_sf <- sites$crds %>% st_as_sf(coords = names(crds), crs = 4326)
  distM   <- st_distance(crds_sf, by_element = F)/1000
  if(parms$bearing) {
    bearM <- matrix(apply(expand.grid(1:nrow(crds_sf), 1:nrow(crds_sf)), 1, function(x) round(lwgeom::st_geod_azimuth(crds_sf[as.numeric(x),])*180/pi,2)),
                    ncol = nrow(crds_sf), nrow = nrow(crds_sf), byrow = T)
  } else bearM <- matrix(0, ncol = nrow(crds_sf), nrow = nrow(crds_sf))

  
  new(
    "migrationSDP",
    Init  = list(
      MinT   = init$minT,
      MaxT   = init$maxT,
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
      ZStdNorm = species$stdNorm$x,
      PStdNorm = do.call(species$stdNorm$dstr, species$stdNorm[-1]),
      xFTReward  = 0:init$MaxX,
      yFTReward  = do.call(init$tReward[[1]], c(x = list(0:init$MaxX), unlist(init$tReward[-1]))),
      decError   = init$dError
    ),
    Sites = list(
      crds  = sites$crds,
      dist  = distM,
      bear  = bearM,
      b0    = c(0,0),
      b1    = c(0,0),
      b2    = c(0,0),
      pred_a1 =  2,
      pred_a2 =  2,
      expend  =  sites$expend,
      gain    =  sites$gain
    ),
    Results = list(
      FitnessMatrix     = NA,
      DecisionMatrix    = NA,
      ProbMatrix        = NA
    )
  )
  
} 

bwdIteration <- function(obj) {

  Init(obj@Init$MaxT,          ## integer
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
       obj@Sites$gain,         ## arma::mat  
       obj@Sites$gain,         ## arma::mat  
       obj@Sites$gain,         ## arma::mat  
       obj@Sites$expend)       ## arma::mat
  
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
