library(tidyverse)
library(tmvtnorm)

###############################################################################
# Proposes value(s) for next iteration of an McMC. Proposals can come from a
# number of different distributions (depending on the value of method).
#
# @param curr the current sample in the McMC
# @param method method by which to propose a new value (e.g. uniform, normal)
# @param pmin minimum values the sampled parameters are allowed to me
# @param pmax maximum values the sampled parameters are allowed to me
# @param pcovar covariance matrix of parameters (for truncated mvtnorm method)
#
# @return list with posterior samples of calibration parameters, along with
# likelihoods, proposals, acceptance rates, and covariances.
###############################################################################
propose_u <- function(curr, method, pmin=NULL, pmax=NULL, pcovar=NULL) {

  if (method=="unifsq") {
    prop <- matrix(runif(2), nrow=1, dimnames=list(NULL, c("pmfp", "ratio")))
    return(list(prop=prop, pr=1))
  } else if (method=="unifadj") {
    ## Propose u_prime from a unit square (centered around current value)
    prop <- matrix(c(runif(n=1, min=curr[1]-min(curr[1], 1-curr[1]), max=curr[1]+min(curr[1], 1-curr[1])),
      runif(n=1, min=curr[2]-min(curr[2], 1-curr[2]), max=curr[2]+min(curr[2], 1-curr[2]))),
      nrow=1, dimnames=list(NULL, c("pmfp", "ratio")))
    return(list(prop=prop, pr=1))
  } else if (method=="tmvnorm") {
    prop <- tmvtnorm::rtmvnorm(n=1, mean=curr, sigma=pcovar, lower=pmin,
      upper=pmax, algorithm="gibbs")
    pr <- tmvtnorm::dtmvnorm(curr, mean=drop(prop), sigma=pcovar, lower=pmin,
      upper=pmax, log=TRUE) -
      tmvtnorm::dtmvnorm(drop(prop), mean=curr, sigma=pcovar, lower=pmin,
        upper=pmax, log=TRUE)
    return(list(prop=prop, pr=pr))
  } else {
    stop("specified method not implemented")
  }
}

propose_logscl <- function(curr, sd) {
  prop <- rnorm(n=1, mean=curr, sd=sd)
  pr <- dnorm(curr, mean=prop, sd=sd, log=TRUE) -
    dnorm(prop, mean=curr, sd=sd, log=TRUE)
  return(list(prop=prop, pr=pr))
}

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
    fd <- fd %>% dplyr::filter(esa==esa_lev, map==fparams[1], time > 0,
      ooc=="no")
  }

  Xmod <- md %>% dplyr::filter(ESA==esa_lev)
  Umod <- Xmod %>% dplyr::select(pmfp, ratio)

  Xfield <- remove_unchanged_points(md=Xmod, fd=fd, tol=tol, quant=quant)

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
    Umod[,i] <- (Umod[,i] - min(Umod[,i]))/(diff(range(Umod[,i]))/scales[i])
  }

  return(list(Xmod=Xmod, Umod=Umod, Zmod=Zmod, Xfield=Xfield, Zfield=Zfield,
    Ofield=Ofield, settings=list(esa=esa_lev, fparams=fparams, pmfp_sc=scales[1],
    ratio_sc=scales[2])))
}

###############################################################################
# Filters out points in the dataset that do not change much over the range of
# all calibration paramters.

# @param md data frame containing data from a computer model
# @param fd data frame containing data from field (e.g. satellite) experiment
# @param tol tolerance of how much data should change to remain in dataset
# @param quant quantile of data to keep in data set
# @param real flag indicating if real data should be used in this calibration
#
# @return data frame containing a reduced set of field (e.g. satellite) data
###############################################################################
remove_unchanged_points <- function(md, fd, tol=NA, quant=NA) {
  if (is.na(tol) && is.na(quant)) {
    stop("either tolerance or quantile must be specified")
  }
  ## Calculate rounded field coordinates nearest to model coordinates
  mdlats <- unique(md$lat)
  mdlons <- unique(md$lon)
  lat_dists <- sqrt(distance(fd$lat, mdlats))
  lon_dists <- sqrt(distance(fd$lon, mdlons))

  fd$rlat <- fd$rlon <- NA
  for (i in 1:nrow(fd)) {
    fd[i,]$rlat <- mdlats[order(lat_dists[i,])[1]]
    fd[i,]$rlon <- mdlons[order(lon_dists[i,])[1]]
  }

  ## Calculate model coordinates that don't change over calibration parameters
  change_data <- md %>%
    dplyr::group_by(lat, lon) %>%
    dplyr::summarise(perc_change = (max(blurred_ena_rate) - min(blurred_ena_rate))/mean(blurred_ena_rate)) %>%
    dplyr::ungroup()

  threshold <- ifelse(is.na(quant), tol, quantile(change_data$perc_change, quant))
  change_data$keep <- ifelse(change_data$perc_change >= threshold, 1, 0)
  res <- merge(fd[,c("lat", "lon", "rlat", "rlon")], change_data,
    by.x=c("rlat", "rlon"), by.y=c("lat", "lon"))
  res <- merge(fd, res[c("lat", "lon", "keep")])
  return(res %>% dplyr::filter(keep==1))
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
  x <- cos(pi*lon/180)*cos(pi*lat/180)
  y <- sin(pi*lon/180)*cos(pi*lat/180)
  z <- sin(pi*lat/180)
  return(data.frame(x, y, z))
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
  lat <- (180/pi)*asin(z)
  lon <- (180/pi)*atan2(y, x)
  return(data.frame(lat, lon))
}
