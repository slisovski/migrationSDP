## Scaling functions
shorebirdScaling     <- function(lbm) {
  
  out <- list(lbm   = lbm,
              mbm   = approx(c(15, 300), c(23.05, 643.4), xout = lbm)$y)
  
  ### Fuel deposition rate
  out$X1gX  <- (out$mbm-out$lbm)/100
  out$X1xkJ <- out$X1gX*32.77
  out$FDR   <- 2.37*(out$lbm/1000)^-0.27
  out$FDRg  <- (out$FDR*(out$lbm))/100
  out$FDRx  <- out$FDRg/out$X1gX
  out$daysToRefuel <- (out$mbm-out$lbm)/out$FDRg
  
  ### Metabolic rate
  out$DER    <- 912*(lbm/1000)^0.704
  out$BMR    <- 5.06*(lbm/1000)^0.729
  out$cond   <- 0.0067*(lbm)^0.4384
  out$Kesm   <- 41-(out$BMR/(1.2*out$cond))
  out$W      <- 57.3*(lbm/1000)^0.81
  
  ### Flight parameters
  x1  <- c(250, 144, 110, 20)
  dst <- c(14000, 7500, 6500, 3000) 
  fit1 <- lm(log(dst)~log(x1))
  pr   <- exp(predict(fit1))
  out$FlightCap <- exp(predict(fit1, newdata = data.frame(x1 = lbm))) + 1000
  
  # Energy exp ~ temp
  t  <- seq(-35, 45, length = 100)
  te <- 41 - (out$BMR/(1.2*out$cond))
  b = out$BMR - (-out$cond)*te
  
  tm <- data.frame(tm = seq(-35, te, length = 100), W = -out$cond*seq(-25, te, length = 100) + b)
  tm <- rbind(tm, data.frame(tm = seq(te, 45, length = 100), W = rep(out$BMR, 100)))
  
  out$EEFnc  <- suppressWarnings(approxfun(x = tm[,1], y = 2*tm[,2]*86.4, rule = 3))
  
  out$speed  <- 1440
  out$c      <- out$FlightCap / (1-(1+99/100)^-0.5)
  
  return(out)
}




#### Plotting
x_site_plot <- function(simu) {
  
  require(ggplot2)
  require(zoo)
  require(patchwork)
  
  alive <- apply(simu[,5,], 1, function(x) all(x==0))
  
  x_traj <- t(apply(simu[alive,3,], 1, function(x) zoo::na.approx(x, rule = 2)))
  x_traj_tab <- tibble(id = 1:sum(alive)) %>% bind_cols(x_traj %>% as_tibble()) %>%
    pivot_longer(cols = -id, values_to = "x") %>% mutate(time = as.numeric(gsub("V", "", name))) %>% dplyr::select(id, time, x)
  
  pl1 <- ggplot(x_traj_tab, aes(x = time, y = x, group = id, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    scale_color_manual(values = rainbow(sum(alive))) +
    theme_light()
  
  site_traj <- simu[alive,2,] + 1
  site_traj_tab <- tibble(id = 1:sum(alive)) %>% bind_cols(site_traj %>% as_tibble()) %>%
    pivot_longer(cols = -id, values_to = "site") %>% mutate(time = as.numeric(gsub("V", "", name))) %>% 
    dplyr::select(id, time, site) %>% filter(!is.na(site))
  
  pl2 <- ggplot(site_traj_tab, aes(x = time, y = site, group = id, color = as.factor(id))) +
    geom_path(show.legend = FALSE) +
    scale_color_manual(values = rainbow(sum(alive))) +
    theme_light()
  
  pl1 / pl2
}


serial_network_plot <- function(mod, simu) {
  
  require(ggplot2)
  
  alive  <- apply(simu[,5,], 1, function(x) all(x==0))
  sts    <- simu[alive,2,] + 1  
  
  stay   <- apply(sts, 1, function(x) {
    as.data.frame(table(x)) %>% as_tibble() %>% setNames(c('site', 'days')) %>% mutate(site = as.numeric(as.character(site))) %>%
      full_join(tibble(site = 1:nrow(mod@Sites$crds)), by = "site") %>% arrange(site) %>% mutate(days = ifelse(site==1 | site==max(site), NA, days)) %>% pull(days)
  }) %>% apply(., 1, sum, na.rm = T)
  
  transM <- expand_grid(a = 1:nrow(mod@Sites$crds), b = 1:nrow(mod@Sites$crds)) %>% 
    left_join(apply(sts, 1, function(x) tibble(a = x[!is.na(x)][-sum(!is.na(x))], b = x[!is.na(x)][-1]) %>% filter(a!=b) %>% mutate(n = 1)) %>%
                do.call("rbind",.) %>% group_by(a,b) %>% summarise_at("n", sum), by = join_by(a, b)) %>% filter(n>0) %>%
    mutate(lon1 = mod@Sites$crd$Lon[a], lat1 = mod@Sites$crd$Lat[a], lon2 = mod@Sites$crd$Lon[b], lat2 = mod@Sites$crd$Lat[b])
  
  ggplot() +
    geom_curve(data = transM, mapping = aes(x = lon1, y = lat1, xend = lon2, yend = lat2, linewidth = n), 
               alpha = 0.7, lineend='round') +
    geom_point(data = crds %>% mutate(days = stay), mapping = aes(x = Lon, y = Lat, size = days), 
               shape = 16, color = "orange") +
    theme_light()
  
  
}
