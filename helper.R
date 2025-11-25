library(tidyverse)
library(tmvtnorm)

###############################################################################
# Preprocesses the data in preparation for calibration. Converts latitude and
# longitude to spherical coordinates. Selects only the ENAs with a specific
# energy level. Scales the data to between 0 and 1.
#
# @param md data frame containing data from a computer model
# @param fd data frame containing data from field (e.g. satellite) experiment
# @param map year(s) (aka map(s)) to use as field data
# @param esa_lev energy level of ENAs to use from the field data
#
# @return list with cleaned data with inputs scaled to 0-1: xm, um, ym, xf, yf,
# exposure times, backgrounds, and settings used in preprocessing
###############################################################################
preprocess_data <- function(md, fd, map, esa_lev) {

  md$pmfp <- md$parallel_mean_free_path
  md <- md[md$ESA == esa_lev,]
  fd <- fd[fd$esa==esa_lev & fd$map %in% map & fd$time > 0,]

  um <- as.matrix(md[,c("pmfp", "ratio")])
  xm <- as.matrix(geo_to_spher_coords(md$lat, md$lon))
  ym <- md$blurred_ena_rate
  xf <- as.matrix(geo_to_spher_coords(fd$ecliptic_lat, fd$ecliptic_lon))
  colnames(xm) <- colnames(um) <- colnames(xf) <- NULL
  ## pull counts, exposure times, and background rates from field data
  yf <- fd$counts
  e <- fd$time
  bg <- fd$background

  xall <- rbind(xm, xf)

  for (i in 1:ncol(xall)) {
    xm[,i] <- (xm[,i] - min(xall[,i]))/diff(range(xall[,i]))
    xf[,i] <- (xf[,i] - min(xall[,i]))/diff(range(xall[,i]))
  }

  for (i in 1:ncol(um)) {
    umin <- min(um[,i])
    umax <- max(um[,i])
    urange <- diff(range(um[,i]))
    um[,i] <- (um[,i] - umin)/(urange)
  }

  return(list(xm=xm, um=um, ym=ym, xf=xf, yf=yf, e=e, bg=bg,
    settings=list(map=map, esa=esa_lev)))
}

###############################################################################
# Converts geographical (lat, lon) coordinates to spherical (x, y, z)
# coordinates

# @param lat vector containing latitudes of observations
# @param lon vector containing longitudes of observations
#
# @return data frame containing observations in spherical coordinates
###############################################################################
geo_to_spher_coords <- function(lat, lon) {
  x <- cos((pi/180)*(lon-180))*cos((pi/180)*lat)
  y <- sin((pi/180)*(lon-180))*cos((pi/180)*lat)
  z <- sin((pi/180)*lat)
  return(data.frame("x"=x, "y"=y, "z"=z))
}

###############################################################################
# Converts spherical (x, y, z) coordinates to geographical (lat, lon)
# coordinates

# @param x vector containing x coordinates of observations
# @param y vector containing y coordinates of observations
# @param z vector containing z coordinates of observations
#
# @return data frame containing observations in geographical coordinates
###############################################################################
spher_to_geo_coords <- function(x, y, z) {
  lon <- (180/pi)*atan2(y, x) + 180
  lat <- (180/pi)*asin(z)
  return(data.frame("lat"=lat, "lon"=lon))
}

###############################################################################
# Modifies longitude values such that the nose of the heliosphere falls in the
# center of sky map visuals
#
# @return a vector of nose centered longitudes
###############################################################################
nose_center_lons <- function(lons) {
  return((lons - 85) %% 360)
}
