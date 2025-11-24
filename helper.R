library(tidyverse)
library(tmvtnorm)

###############################################################################
# Preprocesses the data in preparation for calibration. Converts latitude and
# longitude to spherical coordinates. Selects only the ENAs with a specific
# energy level. Scales the data to between 0 and 1.
#
# @param md data frame containing data from a computer model
# @param fd data frame containing data from field (e.g. satellite) experiment
# @param esa_lev energy level of ENAs to use in this calibration
# @param fparams parameters (pmfp, ratio) of simulated field data to use
# @param scales scaling values for calibration parameters
# @param tol tolerance of how much data should change to remain in dataset
# @param quant quantile of data to keep in data set
# @param real flag indicating if real data should be used in this calibration
#
# @return list with cleaned data prepared for calibration: Xmod, Umod, Zmod,
# Xfield, Zfield, Ofield, and settings used in preprocessing
###############################################################################
preprocess_data <- function(md, fd, esa_lev, fparams, scales, tol=NA,
  quant=NA, real=FALSE, disc=FALSE) {

  md$pmfp <- md$parallel_mean_free_path
  if (!real) {
    fd$pmfp <- fd$parallel_mean_free_path
    if (disc) {
      fd <- fd %>% dplyr::filter(near(scl, fparams[1]), time > 0)
    } else {
      fd <- fd %>% dplyr::filter(esa==esa_lev, pmfp==fparams[1],
        ratio==fparams[2], time > 0)
    }
  } else{
    fd <- fd %>% dplyr::filter(esa==esa_lev, map %in% fparams, time > 0)
  }

  Xmod <- md %>% dplyr::filter(ESA==esa_lev)
  Umod <- Xmod %>% dplyr::select(pmfp, ratio)

  Xfield <- fd

  Xmod[,c("x", "y", "z")] <- geo_to_spher_coords(Xmod$lat, Xmod$lon)
  Xmod <- Xmod %>% dplyr::select(x, y, z)
  Zmod <- md %>%
    dplyr::filter(ESA==esa_lev) %>%
    dplyr::pull(blurred_ena_rate)

  Xfield[,c("x", "y", "z")] <- geo_to_spher_coords(Xfield$lat, Xfield$lon)
  Zfield <- Xfield %>% dplyr::pull(sim_counts)
  Ofield <- Xfield %>%
    dplyr::select(time, background) %>%
    dplyr::rename(bg=background)
  Xfield <- Xfield %>% dplyr::select(x, y, z)

  Xall <- rbind(Xmod, Xfield)

  for (i in 1:ncol(Xall)) {
    Xmod[,i] <- (Xmod[,i] - min(Xall[,i]))/diff(range(Xall[,i]))
    Xfield[,i] <- (Xfield[,i] - min(Xall[,i]))/diff(range(Xall[,i]))
  }

  for (i in 1:ncol(Umod)) {
    Umin <- min(Umod[,i])
    Umax <- max(Umod[,i])
    Urange <- diff(range(Umod[,i]))
    Umod[,i] <- (Umod[,i] - Umin)/(Urange/scales[i])
  }

  return(list(Xmod=Xmod, Umod=Umod, Zmod=Zmod, Xfield=Xfield, Zfield=Zfield,
    Ofield=Ofield, settings=list(esa=esa_lev, fparams=fparams, pmfp_sc=scales[1],
    ratio_sc=scales[2])))
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
