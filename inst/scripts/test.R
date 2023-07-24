library(migrationSDP) 
library(tidyverse)
library(sf)
sf_use_s2(FALSE)
library(stars)

### setup locations
start_end <- list(st_point(c(140, -40)), st_point(c(140, 75))) %>% st_as_sfc() %>% st_set_crs(4326)

point_seq <- st_sample(start_end %>% st_union() %>% st_cast("LINESTRING"), 
                       as.numeric(round(st_distance(start_end)[2,1]/(500*1000),0)),
                       type = "regular")
crds      <- (point_seq %>% st_coordinates())[,1:2] %>% as_tibble() %>% setNames(c("Lon", "Lat"))


### species
spParms <- shorebirdScaling(100)

### expend (temperature)
daily_exp_x <- matrix(2, nrow = nrow(crds), ncol = 365)

### gain
daily_gain_x <- rep(4, nrow(crds))

### simulation time
start   <- "02-15"
end     <- "07-15"
doy_seq <-  as.numeric(format(seq(as.POSIXct(glue::glue("2020-{start}")), as.POSIXct(glue::glue("2020-{end}")), by = "day"), "%j"))

expM  <- daily_exp_x[,doy_seq]
gain  <- daily_gain_x

obj <- make_migrationSDP(init      =  list(minT    = min(doy_seq),
                                           maxT    = max(doy_seq),
                                           MaxX    = 100,
                                           tReward = list(dstr = "dsnorm", mean = 170, sd = 3, xi = 2, factor = 2),
                                           dError  = 70000),
                            species = list(B0      = 3,
                                           w       = 0.028,
                                           xc      = 10,
                                           speed   = spParms$speed,
                                           c       = spParms$c),
                            sites   = list(crds    = crds,
                                           pred    = c(1e-3, 1e-4, 1e-6, 2, 2),
                                           expend  = expM,
                                           gain    = gain,
                                           penalty = c(0, 0.005, 0)),
                            parms   = list(bearing = FALSE))


mod <- bwdIteration(obj)

# plot(raster::raster(mod@Results$FitnessMatrix[,,90]))
# matplot(mod@Results$FitnessMatrix[,nrow(crds),], lwd = 1, lty = 1, type = "l")

simu    <- fwdSimulation(mod, 100, start_t = 1, start_site = 1, start_x = c(30,50))

x_site_plot(simu)
serial_network_plot(mod, simu)

