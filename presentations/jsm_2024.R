
###############################################################################
###############################################################################
#### Code used to produe slides for JSM 2024
###############################################################################
###############################################################################

## global settings
ena_range <- c(0.04174779, 0.18753771)
ibex_palette <- read.csv("../ibex_rgb.csv")
noselongitude <- 256
center <- 180 - (360 - noselongitude)
orig_360 <- seq(0,300,60)
new_360 <- orig_360 - center + 0.01
new_360[new_360 < 0.01] <- new_360[new_360 < 0.01] + 360

###############################################################################
###############################################################################
#### Create a sub-grid of simulation outputs (PMFP = {1250,1500,1750}, RATIO =
#### {0.005, 0.01, 0.02})
###############################################################################
###############################################################################

library(ggplot2)

## Read in simulation data
md <- read.csv("../data/sims.csv")

## Filter to restricted range of PMFP and ratio
esa_lev <- 4
min_pmfp <- 1250
max_pmfp <- 1750
min_ratio <- 0.005
max_ratio <- 0.02
plot_data <- md[md$ESA==esa_lev & md$parallel_mean_free_path >= min_pmfp &
  md$parallel_mean_free_path <= max_pmfp & md$ratio >= min_ratio &
  md$ratio <= max_ratio,]

## Plot grid of simulator outputs
pdf("sim_subgrid.pdf", width=12, height=5.75)
ggplot(data=plot_data) +
  facet_grid(ratio ~ parallel_mean_free_path) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=blurred_ena_rate)) +
  scale_fill_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360,
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL,
     name="Parallel Mean Free Path")) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45),
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL, name="Ratio")) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="bottom",
    legend.title=element_text(margin=margin(r=10), size=16),
    plot.title=element_text(hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=20),
    strip.text.y=element_text(size=20),
    axis.title.y.right=element_text(size=20, vjust=1),
    axis.title.x.top=element_text(size=20, vjust=1),
    axis.title=element_text(size=20)) +
  coord_fixed(ratio=1)
dev.off()

###############################################################################
###############################################################################
#### Computational time to fit a GP surrogate on large-scale data
###############################################################################
###############################################################################

###############################################################################
###############################################################################
#### Vecchia approximation visual
###############################################################################
###############################################################################

library(lhs)
library(plgp)

m <- 10
n <- 60

set.seed(8907198)
X <- randomLHS(n=n, k=2)
X <- cbind(X, 1:nrow(X))
randX <- sample((m+1):nrow(X),1)
plot(x=X[,1], y=X[,2], type='n', xlab=expression("x"[1]),
  ylab=expression("x"[2]), cex.lab=1.5, cex.axis=1.25)
points(x=X[randX,1], y=X[randX,2], pch=18, col="lightblue", cex=5.0)
nns <- order(distance(X[randX,1:2,drop=FALSE], X[1:(randX-1),1:2]))[1:m]
points(x=X[nns,1], y=X[nns,2], pch=15, col="lightgrey", cex=4.0)
text(x=X[,1], y=X[,2], labels=X[,3], cex=1.25)

###############################################################################
###############################################################################
#### Scaled Vecchia approximation visual
###############################################################################
###############################################################################

library(lhs)

n <- 30
n1 <- n*0.1
n2 <- n - n1
nbig <- 20*n
m <- 5

set.seed(8907198)
Xinit1 <- randomLHS(n=n1, k=2)
Xbig1 <- randomLHS(n=nbig, k=2)
Xmm <- matrix(data=NA, nrow=0, ncol=2)

for (i in 1:n2) {
  mindist <- 0
  minind <- NA
  for (j in 1:nrow(Xbig1)) {
    xnew <- Xbig1[j,,drop=FALSE]
    Xpot <- rbind(Xinit1, Xmm, xnew)
    dist <- distance(Xpot)
    itermin <- min(dist[upper.tri(dist)])
    if (itermin > mindist) {
      mindist <- itermin
      minind <- j
    }
  }
  Xmm <- rbind(Xmm, Xbig1[minind,,drop=FALSE])
  Xbig1 <- Xbig1[-minind,]
}
Xfin <- rbind(Xinit1, Xmm)
Xfin <- cbind(Xfin, 1:nrow(Xfin))

randX <- 30
pdf("neighborhood_pre.pdf", height=5, width=5)
plot(x=Xfin[,1], y=Xfin[,2], type='n', xlab=expression("x"[1]), ylab="",
 cex.lab=1.25, cex.main=1.25,
 main="neighborhood creation pre-scaling")
title(ylab=expression("x"[2]), line=2.5, cex.lab=1.25)
points(x=Xbig1[,1], y=Xbig1[,2], pch=16, cex=0.3, col="grey")
points(x=Xfin[randX,1], y=Xfin[randX,2], pch=18, col="lightblue", cex=4.0)
nns <- order(distance(Xfin[randX,1:2,drop=FALSE], Xfin[1:(randX-1),1:2]))[1:m]
points(x=Xfin[nns,1], y=Xfin[nns,2], pch=15, col="lightgrey", cex=2.75)
text(x=Xfin[,1], y=Xfin[,2], labels=Xfin[,3])
dev.off()

scl1 <- 1.3
scl2 <- 0.7
Xfinsc <- Xfin
Xfinsc[,1] <- Xfinsc[,1]*scl1
Xfinsc[,2] <- Xfinsc[,2]*scl2

pdf("neighborhood_scl_after.pdf", height=5*scl2, width=5*scl1)
plot(x=Xfinsc[,1], y=Xfinsc[,2], type='n', xlab=expression("x"[1]* " scaled"),
  ylab="", cex.lab=1.25, cex.main=1.25,
  main="scaled inputs after neighborhood creation")
title(ylab=expression("x"[2]* " scaled"), line=2.5, cex.lab=1.25)
points(x=Xbig1[,1]*scl1, y=Xbig1[,2]*scl2, pch=16, cex=0.3, col="grey")
points(x=Xfinsc[randX,1], y=Xfinsc[randX,2], pch=18, col="lightblue", cex=4.0)
points(x=Xfinsc[nns,1], y=Xfinsc[nns,2], pch=15, col="lightgrey", cex=2.75)
text(x=Xfinsc[,1], y=Xfinsc[,2], labels=Xfinsc[,3])
dev.off()

pdf("neighborhood_post.pdf", height=5*scl2, width=5*scl1)
plot(x=Xfinsc[,1], y=Xfinsc[,2], type='n', xlab=expression("x"[1]* " scaled"),
  ylab="", cex.lab=1.25, cex.main=1.25,
  main="neighborhood creation post-scaling")
title(ylab=expression("x"[2]* " scaled"), line=2.5, cex.lab=1.25)
points(x=Xbig1[,1]*scl1, y=Xbig1[,2]*scl2, pch=16, cex=0.3, col="grey")
points(x=Xfinsc[randX,1], y=Xfinsc[randX,2], pch=18, col="lightblue", cex=4.0)
nns <- order(distance(Xfinsc[randX,1:2,drop=FALSE], Xfinsc[1:(randX-1),1:2]))[1:m]
points(x=Xfinsc[nns,1], y=Xfinsc[nns,2], pch=15, col="lightgrey", cex=2.75)
text(x=Xfinsc[,1], y=Xfinsc[,2], labels=Xfinsc[,3])
dev.off()

###############################################################################
###############################################################################
#### Out-of-sample prediction visual for Scaled Vecchia GP surrogate. 5x5 grid
#### of skymaps.
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

md <- read.csv("../data/sims.csv")
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
Xmd <- md[,c("x", "y", "z")]
Umd <- md[,c("parallel_mean_free_path", "ratio")]
Zmd <- md$blurred_ena_rate

for (i in 1:ncol(Xmd)) {
  Xmd[,i] <- (Xmd[,i] - min(Xmd[,i]))/diff(range(Xmd[,i]))
}

for (i in 1:ncol(Umd)) {
  Umin <- min(Umd[,i])
  Umax <- max(Umd[,i])
  Urange <- diff(range(Umd[,i]))
  Umd[,i] <- (Umd[,i] - Umin)/(Urange)
}

fit <- fit_scaled(y=Zmd, inputs=as.matrix(cbind(Xmd, Umd)), nug=1e-4, ms=25)

pmfps_new <- seq(1250, 1750, by=125)
ratios_new <- c(0.005, 0.0075, 0.01, 0.015, 0.02)
calib_new_grid <- as.matrix(expand.grid(pmfps_new, ratios_new))
colnames(calib_new_grid) <- c("pmfp", "ratio")
calib_new_grid <- calib_new_grid[-c(1,3,5,11,13,15,21,23,25),]
calib_new_grid_scl <- calib_new_grid
calib_new_grid_scl[,1] <- (calib_new_grid[,1]-min(calib_new_grid[,1]))/diff(range(calib_new_grid[,1]))
calib_new_grid_scl[,2] <- (calib_new_grid[,2]-min(calib_new_grid[,2]))/diff(range(calib_new_grid[,2]))

X <- distinct(md[,c("lat", "lon", "ecliptic_lon_center")])
res <- matrix(NA, nrow=0, ncol=6)

for (i in 1:nrow(calib_new_grid)) {
  XX <- as.matrix(cbind(distinct(Xmd), calib_new_grid_scl[i,,drop=FALSE]))
  colnames(XX) <- c("x", "y", "z", "parallel_mean_free_path", "ratio")
  preds <- predictions_scaled(fit, XX, m=25, joint=FALSE, predvar=FALSE)
  res <- rbind(res, cbind(X, calib_new_grid[i,,drop=FALSE], preds))
  print(i)
}

sims <- md[md$ESA==4 & md$parallel_mean_free_path %in% c(1250, 1500, 1750),]
sims <- sims[sims$ratio %in% c(0.005, 0.01, 0.02),]
sims <- sims[,c("lat", "lon", "ecliptic_lon_center", "parallel_mean_free_path", "ratio", "blurred_ena_rate")]
colnames(sims) <- c(colnames(sims)[1:3], "pmfp", "ratio", "preds")

res <- rbind(res, sims)

pdf("surr_sim_grid.pdf", width=12, height=5.75)
ggplot(data=res) +
  facet_grid(ratio ~ pmfp) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=preds)) +
  scale_fill_gradientn(colors=ibex_palette$hex,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360,
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL,
     name="Parallel Mean Free Path")) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45),
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL, name="Ratio")) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="bottom",
    legend.title=element_text(margin=margin(r=10), size=16),
    plot.title=element_text(hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    strip.text.x=element_text(size=15),
    strip.text.y=element_text(size=15),
    axis.title.y.right=element_text(size=20, vjust=1),
    axis.title.x.top=element_text(size=20, vjust=1),
    axis.title=element_text(size=20)) +
  coord_fixed(ratio=1)
dev.off()

###############################################################################
###############################################################################
#### Additionally, compute metrics (RMSE, RMSPE, coverage, and execution time)
#### for hold-one-out on simulator runs
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

md <- read.csv("../data/sims.csv")
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
Xmd <- md[,c("x", "y", "z")]
Umd <- md[,c("parallel_mean_free_path", "ratio")]
Zmd <- md$blurred_ena_rate

unique_calibs <- distinct(md[,c("parallel_mean_free_path", "ratio")])
rmses <- rmspes <- covs <- exec_times <- c()

for (k in 1:30) {
  for (i in 1:nrow(unique_calibs)) {
    iter_pmfp <- unique_calibs[i,1]
    iter_ratio <- unique_calibs[i,2]
    iter_pmfp_scl <- (iter_pmfp - 500)/2500
    iter_ratio_scl <- (iter_ratio - 0.001)/(0.1-0.001)

    iter_holdout <- md[md$parallel_mean_free_path==iter_pmfp & md$ratio==iter_ratio,]
    iter_md <- md[!(md$parallel_mean_free_path==iter_pmfp & md$ratio==iter_ratio),]

    Xmd <- iter_md[,c("x", "y", "z")]
    Umd <- iter_md[,c("parallel_mean_free_path", "ratio")]
    Zmd <- iter_md$blurred_ena_rate

    for (j in 1:ncol(Xmd)) {
      Xmd[,j] <- (Xmd[,j] - min(Xmd[,j]))/diff(range(Xmd[,j]))
    }

    for (j in 1:ncol(Umd)) {
      Umin <- min(Umd[,j])
      Umax <- max(Umd[,j])
      Urange <- diff(range(Umd[,j]))
      Umd[,j] <- (Umd[,j] - Umin)/(Urange)
    }

    iter_fit <- fit_scaled(y=Zmd, inputs=as.matrix(cbind(Xmd, Umd)), nug=1e-4, ms=25)
    XX <- cbind(distinct(Xmd), matrix(c(iter_pmfp_scl, iter_ratio_scl), nrow=1))
    tic <- proc.time()[3]
    preds <- predictions_scaled(iter_fit, as.matrix(XX), m=25, joint=FALSE, predvar=TRUE)
    toc <- proc.time()[3]
    exec_times <- c(exec_times, toc-tic)
    rmses <- c(rmses, sqrt(mean((preds$means-iter_holdout$blurred_ena_rate)^2)))
    rmspes <- c(rmspes, sqrt(mean((100*(preds$means-iter_holdout$blurred_ena_rate)/iter_holdout$blurred_ena_rate)^2)))
    covs <- c(covs, mean(iter_holdout$blurred_ena_rate <= preds$means+qnorm(0.975, 0, sqrt(preds$vars)) &
      iter_holdout$blurred_ena_rate >= qnorm(0.025, 0, sqrt(preds$vars))))
    cat("\nFinished calib combo: ", i)
  }
  cat("\nFinished mc rep: ", k)
}

pdf("rmses.pdf", height=3.5, width=5)
hist(rmses, xlab="RMSE", ylab="", main="", cex.lab=1.5)
dev.off()

pdf("rmspes.pdf", height=3.5, width=5)
hist(rmspes, xlab="RMSPE", ylab="", main="", cex.lab=1.5)
dev.off()

pdf("covs.pdf", height=3.5, width=5)
hist(covs, xlab="coverage", ylab="", main="", cex.lab=1.5)
dev.off()

pdf("times.pdf", height=3.5, width=5)
hist(exec_times, xlab="execution time (seconds)", ylab="", main="", cex.lab=1.5)
dev.off()



###############################################################################
###############################################################################
#### Create "unknown" observed sky map based on the truth of one simulator
#### output.
###############################################################################
###############################################################################

source("../helper.R")

md <- read.csv("../data/sims.csv")

sim_pmfp <- 1750
sim_ratio <- 0.02
year <- "2009A"
esa_lev <- 4
real_data <- read.csv("../data/ibex_real.csv")
real_data <- real_data[real_data$map==year & real_data$esa==esa_lev,]
sim_dat <- md[md$ESA==esa_lev & md$parallel_mean_free_path==sim_pmfp &
  md$ratio==sim_ratio,c("ecliptic_lon_center", "lat", "lon",
   "blurred_ena_rate")]

sim_dat[,c("x", "y", "z")] <- geo_to_spher_coords(sim_dat$lat, sim_dat$lon)

real_data$sim_ena_rate <- NA
sim_real_map <- read.csv("../data/ibex_sim_map_2009A.csv")
sim_real_map <- sim_real_map[order(sim_real_map$row),]

for (i in 1:nrow(real_data)) {
  near_sim <- sim_real_map[i,2]
  real_data$sim_ena_rate[i] <- sim_dat[near_sim,c("blurred_ena_rate")]
}

real_data$counts <- rpois(real_data$time*(real_data$sim_ena_rate +
 real_data$background), real_data$time*(real_data$sim_ena_rate +
 real_data$background))
real_data$est_rate <- real_data$counts/real_data$time -
 real_data$background
real_data$est_rate <-
 ifelse(is.nan(real_data$est_rate) | real_data$est_rate < 0, 0,
  real_data$est_rate)

pdf("unknown_sim_obs.pdf", height=3.5, width=6)
ggplot(data=real_data, aes(x=ecliptic_lon_center, y=ecliptic_lat)) +
  geom_point(aes(col=est_rate), size=1.5)+
  scale_color_gradientn(colours=ibex_palette$hex,
    name="Estimated ENA Rate", limits=ena_range,
    oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude") + ylab("Latitude")+
  theme_bw() +
  theme(legend.position="n",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15))+
  ggtitle("'Unknown' IBEX satellite data")+
  coord_fixed(ratio=1)
dev.off()

###############################################################################
###############################################################################
#### Create simulator output at true calibration settings
#### Create surrogate output at estimated calibration settings
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

md <- read.csv("../data/sims.csv")
md[,c("x", "y", "z")] <- geo_to_spher_coords(md$lat, md$lon)
Xmd <- md[,c("x", "y", "z")]
Umd <- md[,c("parallel_mean_free_path", "ratio")]
Zmd <- md$blurred_ena_rate

for (i in 1:ncol(Xmd)) {
  Xmd[,i] <- (Xmd[,i] - min(Xmd[,i]))/diff(range(Xmd[,i]))
}

for (i in 1:ncol(Umd)) {
  Umin <- min(Umd[,i])
  Umax <- max(Umd[,i])
  Urange <- diff(range(Umd[,i]))
  Umd[,i] <- (Umd[,i] - Umin)/(Urange)
}

fit <- fit_scaled(y=Zmd, inputs=as.matrix(cbind(Xmd, Umd)), nug=1e-4, ms=25)
est_pmfp <- (1558.13 - 500)/2500
est_ratio <- (0.0348 - 0.001)/(0.1-0.001)
XX <- cbind(distinct(Xmd), matrix(c(est_pmfp, est_ratio), nrow=1))
colnames(XX) <- c("x", "y", "z", "parallel_mean_free_path", "ratio")
lambda_preds <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)

plot_data <- cbind(md[md$parallel_mean_free_path==1750 & md$ratio==0.02,
  c("ecliptic_lon_center", "lat", "blurred_ena_rate")], lambda_preds)

## Simulator output
pdf("sim_ind_output.pdf", height=2.96, width=5.08)
ggplot(data=plot_data) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=blurred_ena_rate)) +
  scale_fill_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="n",
        axis.title=element_text(size=13),
        plot.title=element_text(hjust=0.5, size=15))+
  ggtitle("Simulator output (pmfp = 1750, ratio = 0.02)")
dev.off()

## Surrogate output at estimated calibration settings
pdf("surr_ind_output.pdf", height=2.96, width=5.08)
ggplot(data=plot_data) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=lambda_preds)) +
  scale_fill_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="n",
        axis.title=element_text(size=13),
        plot.title=element_text(hjust=0.5, size=15))+
  ggtitle("Surrogate output (pmfp = 1558.13, ratio = 0.0348)")
dev.off()

###############################################################################
###############################################################################
#### Create surrogate output at estimated calibration settings for real data
###############################################################################
###############################################################################

library(ggplot2)
library(gridExtra)
library(dplyr)

source("../helper.R")
source("../vecchia_scaled.R")

md <- read.csv("../data/sims.csv")
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real$sim_counts <- ibex_real$counts
ibex_real$lat <- ibex_real$ecliptic_lat
ibex_real$lon <- ibex_real$ecliptic_lon

pd <- preprocess_data(md=md, fd=ibex_real, esa_lev=4,
  fparams=data.frame(year="2009A"), scales=c(1,1), tol=NA, quant=0.0,
  real=TRUE, disc=FALSE)
iter_fit <- fit_scaled(y=pd$Zm, inputs=as.matrix(cbind(pd$Xm, pd$Um)),
 nug=1e-4, ms=25)

years <- c("2009A", "2010A", "2011A")
for (y in years) {
  fp <- paste0("mcmc_res_esa4_nmcmc10000_end25_scNA_quantNA_svecchia_real", y, "_[0-9]*.rds")
  fn <- list.files(path="../results", pattern=fp)
  calib_res <- readRDS(paste0("../results/", fn))
  iter_real <- ibex_real[ibex_real$map==y,]

  est_pmfp <- mean(calib_res$mcmc_res$u[seq(1001,10000,by=10),1])
  est_ratio <- mean(calib_res$mcmc_res$u[seq(1001,10000,by=10),2])
  est_scale <- mean(exp(calib_res$mcmc_res$logscl[seq(1001,10000,by=10)]))

  XX_unit <- cbind(distinct(pd$Xm), matrix(c(est_pmfp, est_ratio), nrow=1))
  colnames(XX_unit) <- c("x", "y", "z", "pmfp", "ratio")
  lhat_pred <- predictions_scaled(fit, as.matrix(XX_unit), m=25, joint=FALSE,
   predvar=FALSE)
  lhat_pred_scale <- lhat_pred*est_scale

  est_plot_data <- cbind(distinct(model_data[,c("lat", "lon")]),
   lhat_pred, lhat_pred_scale)
  obs_plot_data <- cbind(iter_real[,c("lat", "lon")],
    iter_real$sim_counts/iter_real$time - iter_real$background)
  colnames(est_plot_data) <- c("lat", "lon", "rate", "scl_rate")
  colnames(obs_plot_data) <- c("lat", "lon", "rate")
  est_plot_data$nlon <- nose_center_lons(est_plot_data$lon)
  obs_plot_data$nlon <- nose_center_lons(obs_plot_data$lon)

  ## Calculate difference between observed and predicted
  real_sim_map <- read.csv(paste0("../data/ibex_sim_map_", y, ".csv"))
  real_sim_map <- real_sim_map[order(real_sim_map$row),]
  obs_plot_data$near_sim_rate <- obs_plot_data$pred_diff <-
    obs_plot_data$diff_w_disc <- NA
  for (i in 1:nrow(obs_plot_data)) {
    obs_plot_data$near_sim_rate[i] <- lhat_pred_scale[real_sim_map$near_sim[i]]
    obs_plot_data$near_sim_scl_rate[i] <- lhat_pred_scale[real_sim_map$near_sim[i]]*est_scale
  }

  obs_plot_data$pred_diff <- obs_plot_data$near_sim_rate - obs_plot_data$rate
  obs_plot_data$diff_w_disc <- obs_plot_data$near_sim_scl_rate - obs_plot_data$rate

  predrange <- c(0.04174779, max(c(est_plot_data$rate, est_plot_data$scl_rate)))

  ## Estimated map (w/ scale included)
  pdf(paste0("est_w_scale_", y, ".pdf"), width=5.08, height=2.96)
  ggplot(est_plot_data, aes(x=nlon, y=lat)) +
    geom_point(aes(col=scl_rate), size = 1.75) +
    scale_color_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
      limits=predrange) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60), label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(legend.position="right") +
    labs(x = "Longitude", y = "Latitude", col = "lambda") +
    ggtitle(paste0(y, ": estimated ENA rate (after discrepancy)"))
  dev.off()

  ## Estimated map (w/o scale included)
  pdf(paste0("est_wo_scale_", y, ".pdf"), width=5.08, height=2.96)
  ggplot(est_plot_data, aes(x=nlon, y=lat)) +
    geom_point(aes(col=rate), size = 1.75) +
    scale_color_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
      limits=predrange) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60), label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(legend.position="right") +
    labs(x = "Longitude", y = "Latitude", col = "lambda") +
    ggtitle(paste0(y, ": estimated ENA rate"))
  dev.off()

  ## Observed rates
  pdf(paste0("obs_rates_", y, ".pdf"), width=5.08, height=2.96)
  ggplot(obs_plot_data, aes(x=nlon, y=lat)) +
    geom_point(aes(col=rate), size = 0.75) +
    scale_color_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
      limits=predrange, oob=scales::oob_squish) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60),
      label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(legend.position="left") +
    labs(x = "Longitude", y = "Latitude", col = "lambda") +
    ggtitle(paste0(y, ": counts/time - background"))
  dev.off()


  max_diff <- max(abs(range(obs_plot_data$pred_diff)))
  ## Difference
  pdf(paste0("diff_", y, ".pdf"), width=5.08, height=2.96)
  ggplot(obs_plot_data, aes(x=nlon, y=lat)) +
    geom_point(aes(col=pred_diff), size=0.5) +
    scale_color_gradientn(colours=colorRampPalette(c("blue", "white", "red"))(128),
      limits=c(-max_diff, max_diff)) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60),
      label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(legend.position="left") +
    labs(x = "Longitude", y = "Latitude", col = "diff") +
    ggtitle(paste0(y, ": predicted rate - observed rate"))
  dev.off()

  pdf(paste0("diff_w_disc", y, ".pdf"), width=5.08, height=2.96)
  ggplot(obs_plot_data, aes(x=nlon, y=lat)) +
    geom_point(aes(col=diff_w_disc), size=0.5) +
    scale_color_gradientn(colours=colorRampPalette(c("blue", "white", "red"))(128),
      limits=c(-max_diff, max_diff)) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60),
      label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(legend.position="left") +
    labs(x = "Longitude", y = "Latitude", col = "abs diff") +
    ggtitle(paste0(y, ": predicted rate - observed rate"))
  dev.off()
}

###############################################################################
###############################################################################
#### Create simulated discrepancy sky maps
###############################################################################
###############################################################################

source("../helper.R")
md <- read.csv("../data/sims.csv")
ibex_real <- read.csv("../data/ibex_real.csv")

sim_pmfp <- 2250
sim_ratio <- 0.05
disc <- 0.1
year <- "2009A"
esa_lev <- 4
ibex_real <- ibex_real[ibex_real$map==year & ibex_real$esa==esa_lev,]
sim_dat <- md[md$ESA==esa_lev & md$parallel_mean_free_path==sim_pmfp &
  md$ratio==sim_ratio,c("ecliptic_lon_center", "lat", "lon",
   "blurred_ena_rate")]
sim_dat[,c("x", "y", "z")] <- geo_to_spher_coords(sim_dat$lat, sim_dat$lon)

ibex_real$sim_ena_rate <- NA
sim_real_map <- read.csv("../data/ibex_sim_map_2009A.csv")
sim_real_map <- sim_real_map[order(sim_real_map$row),]

for (i in 1:nrow(ibex_real)) {
  near_sim <- sim_real_map[i,2]
  ibex_real$sim_ena_rate[i] <- sim_dat[near_sim,c("blurred_ena_rate")]*disc
}

ibex_real$counts <- rpois(ibex_real$time*(ibex_real$sim_ena_rate +
 ibex_real$background), ibex_real$time*(ibex_real$sim_ena_rate +
 ibex_real$background))
ibex_real$est_rate <- ibex_real$counts/ibex_real$time -
 ibex_real$background
ibex_real$est_rate <-
 ifelse(is.nan(ibex_real$est_rate) | ibex_real$est_rate < 0, 0,
  ibex_real$est_rate)

pdf("sim_ind_output.pdf", height=2.96, width=5.08)
ggplot(data=plot_data) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=blurred_ena_rate)) +
  scale_fill_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="n",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15))+
  ggtitle("Simulator output (pmfp = 1750, ratio = 0.02)")
dev.off()

 pdf("unknown_sim_obs.pdf", height=3.5, width=6)
ggplot(data=real_data, aes(x=ecliptic_lon_center, y=ecliptic_lat)) +
  geom_point(aes(col=est_rate), size=1.5)+
  scale_color_gradientn(colours=ibex_palette$hex,
    name="Estimated ENA Rate", limits=ena_range,
    oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude") + ylab("Latitude")+
  theme_bw() +
  theme(legend.position="n",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15))+
  ggtitle("'Unknown' IBEX satellite data")+
  coord_fixed(ratio=1)
dev.off()
