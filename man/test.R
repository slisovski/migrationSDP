library(migrationSDP) 
library(dplyr)
# library(sf)

### setup locations
start_end <- list(st_point(c(140, -40)), st_point(c(140, 75))) %>% st_as_sfc() %>% st_set_crs(4326)

point_seq <- st_sample(start_end %>% st_union() %>% st_cast("LINESTRING"), 
                       as.numeric(round(st_distance(start_end)[2,1]/(500*1000),0)),
                       type = "regular")
crds   <- (point_seq %>% st_coordinates())[,1:2] %>% as_tibble() %>% setNames(c("Lon", "Lat"))


### expend/gain (example)
expendM <- matrix(2, nrow = nrow(crds), ncol = length(1:100))
gainM   <- matrix(4, nrow = nrow(crds), ncol = length(1:100))


obj <- make_migrationSDP(init =       list(minT    = 1,
                                           maxT    = 100,
                                           MaxX    = 100,
                                           tReward = list(dstr = "dsnorm", mean = 80, sd = 2, xi = 2),
                                           dError  = 7000),
                            species = list(B0      = 3,
                                           w       = 0.028,
                                           xc      = 10,
                                           stdNorm = list(dstr = "dnorm", x = seq(-2.5, 2.5, .5), mean = 0, sd = 2),
                                           speed   = 1440,
                                           c       = 14200),
                            sites   = list(crds    = crds,
                                           expend  = expendM,
                                           gain    = gainM,
                                           penalty = c(0, -0.05, 0)),
                            parms   = list(bearing = FALSE))


mod <- bwdIteration(obj)
plot(raster::raster(mod@Results$FitnessMatrix[,,100]))
matplot(mod@Results$FitnessMatrix[,nrow(crds),], lwd = 1, lty = 1, type = "l")