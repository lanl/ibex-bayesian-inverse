
###############################################################################
###############################################################################
#### Code used to produce graphics for CoDA 2025 poster
###############################################################################
###############################################################################

## global settings
ibex_palette <- read.csv("ibex_rgb.csv")
ibex_palette <- ibex_palette[2:255,]
noselongitude <- 256
center <- 180 - (360 - noselongitude)
orig_360 <- seq(0,300,60)
new_360 <- orig_360 - center + 0.01
new_360[new_360 < 0.01] <- new_360[new_360 < 0.01] + 360

###############################################################################
###############################################################################
#### MOTIVATION SECTION
###############################################################################
###############################################################################

#### Plot raw data from satellite for one map
year <- "2020A"
esa_lev <- 4
real_data <- read.csv("../data/ibex_real.csv")
real_data <- real_data[real_data$map==year & real_data$esa==esa_lev,]
real_data$est_rate <- real_data$counts / real_data$time - real_data$background
ena_range <- quantile(real_data$est_rate, probs=c(0.025, 0.975))

pdf("satellite_data_2020A.pdf", height=5.0, width=7.5)
ggplot(data=real_data, aes(x=ecliptic_lon_center, y=ecliptic_lat)) +
  geom_point(aes(col=est_rate), size=0.95)+
  scale_color_gradientn(colours=ibex_palette$hex,
    name="Estimated ENA Rate", limits=ena_range, oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude") + ylab("Latitude")+
  theme_bw() +
  theme(legend.position="n",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15),
        plot.background=element_rect(fill='transparent', color=NA),
        panel.background=element_rect(fill='transparent'))+
  ggtitle("IBEX satellite data")+
  coord_fixed(ratio=1)
dev.off()


#### Create a sub-grid of simulation outputs (PMFP = {5000,1750,30000}, RATIO =
#### {0.001, 0.02, 0.1})

library(ggplot2)

## Read in simulation data
md <- read.csv("../data/sims.csv")

## Filter to restricted range of PMFP and ratio
esa_lev <- 4
pmfps <- c(500, 1750, 3000)
ratios <- c(0.001, 0.02, 0.1)
plot_data <- md[md$ESA==esa_lev & md$parallel_mean_free_path %in% pmfps &
  md$ratio %in% ratios,]

## Plot grid of simulator outputs
pdf("sim_subgrid.pdf", width=12, height=5.75)
ggplot(data=plot_data) +
  facet_grid(ratio ~ parallel_mean_free_path) +
  geom_raster(aes(x=ecliptic_lon_center, y=lat, fill=blurred_ena_rate)) +
  scale_fill_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate", oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=NULL,
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL,
     name="")) +
  scale_y_continuous(expand=c(.0,.0), labels=NULL, breaks=seq(-45,45,45),
    sec.axis=sec_axis(~., labels=NULL, breaks=NULL, name="")) +
  xlab("")+ ylab("") +
  theme(legend.position="none",
    plot.title=element_text(hjust=0.5),
    panel.grid.major=element_blank(),
    panel.grid.minor=element_blank(),
    panel.background=element_rect(fill='transparent'),
    plot.background=element_rect(fill='transparent', color=NA),
    strip.background=element_rect(fill='transparent'),
    legend.background=element_rect(fill='transparent'),
    strip.text.x=element_blank(),
    strip.text.y=element_blank(),
    axis.ticks=element_blank(),
    axis.title.y.right=element_text(size=20, vjust=1),
    axis.title.x.top=element_text(size=20, vjust=1),
    axis.title=element_text(size=20)) +
  coord_fixed(ratio=1)
dev.off()




###############################################################################
###############################################################################
#### VECCHIA SECTION
###############################################################################
###############################################################################

#### Vecchia approximation visual

library(lhs)
library(plgp)

m <- 5
n <- 40

pdf("vecchia_visual.pdf", height=5.0, width=5.0)
set.seed(8907198)
X <- randomLHS(n=n, k=2)
X <- cbind(X, 1:nrow(X))
randX <- sample((m+1):nrow(X),1)
plot(x=X[,1], y=X[,2], type='n', xlab=expression("x"[1]),
  ylab=expression("x"[2]), cex.lab=1.5, cex.axis=1.25)
points(x=X[randX,1], y=X[randX,2], pch=18, col="lightblue", cex=5.0)
nns <- order(distance(X[randX,1:2,drop=FALSE], X[1:(randX-1),1:2]))[1:m]
points(x=X[nns,1], y=X[nns,2], pch=15, col="lightgrey", cex=4.0)
text(x=X[,1], y=X[,2], labels=X[,3], cex=1.1)
dev.off()

###############################################################################
###############################################################################
#### POISSON COMPUTER MODEL CALIBRATION SECTION
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

sims <- read.csv("../data/coda_ibex2020Asims.csv")
sims$est_rate <- sims$sim_counts / sims$time - sims$background
sims$nlon <- nose_center_lons(sims$lon)

pdf("synth_ibex_data.pdf", height=2.82, width=3.62)
ggplot(data=sims, aes(x=nlon, y=lat)) +
  geom_point(aes(col=est_rate), size=0.36) +
  scale_colour_gradientn(colors=ibex_palette$hex, limits=ena_range,
   name="ENA Rate") +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360) +
  scale_y_continuous(expand=c(.0,.0), breaks=seq(-45,45,45)) +
  xlab("Longitude")+ ylab("Latitude") +
  theme(legend.position="n",
        axis.title=element_text(size=13),
        plot.title=element_text(hjust=0.5, size=15),
        plot.background=element_rect(fill='transparent', color=NA),
        panel.background=element_rect(fill='transparent'))+
  ggtitle("Synthetic IBEX data")
dev.off()

blue_rgb <- col2rgb("#6699CC")
t_col <- rgb(blue_rgb[1], blue_rgb[2], blue_rgb[3],
  max=255, alpha=(100 - 60) * 255 / 100)
pdf("bivariate_dist.pdf", height=2.82, width=3.97)
par(mar=c(3.15, 3.15, 0.5, 0.5))
plot(x=mcmc_res$u[1001:10000,1]*2500+500, y=mcmc_res$u[1001:10000,2]*(0.1-0.001)+0.001,
 xlim=c(500, 3000), ylim=c(0, 0.1), xlab="", ylab="", col=t_col)
abline(v=seq(500, 3000, by=250), col="lightgrey", lty=2)
abline(h=seq(0, 0.1, by=0.02), col="lightgrey", lty=2)
points(x=1750, y=0.02, col=2, lwd=2, pch=8, cex=1.5)
title(ylab="Ratio", line=2.15)
title(xlab="Parallel Mean Free Path", line=2.15)
legend("topleft", c("posterior draws", "truth"), pch=c(1, 8),
  col=c(4, 2), lwd=2, bg="white", lty=NA, cex=0.9)
dev.off()


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

#### Simulator output at true calibration settings
plot_data <- cbind(md[md$parallel_mean_free_path==1750 & md$ratio==0.02,
  c("ecliptic_lon_center", "lat", "blurred_ena_rate")], lambda_preds)

pdf("sim_ind_output.pdf", height=2.82, width=3.62)
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
  ggtitle("Simulator output at 'truth'")
dev.off()

## Surrogate output at estimated calibration settings
pdf("surr_ind_output.pdf", height=2.82, width=3.62)
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
  ggtitle("Surrogate output at estimates")
dev.off()



###############################################################################
###############################################################################
#### POISSON DEEP GAUSSIAN PROCESS SECTION
###############################################################################
###############################################################################

source("../simulate_data.R")
source("poiss_test.R")

booth_barnett <- function(x) {
  if (is.null(nrow(x))) {
    x <- matrix(x, ncol=1)
  }
  res <- drop(2*sin(pi*x*6) + 2*cos(pi*x*12))
  res [x > 0.4615] <- log(4.84)
  return(res)
}
n <- 25
N <- 200
X <- matrix(seq(0, 1, length=n), ncol=1)
XX <- matrix(seq(0, 1, length=N), ncol=1)
t <- 50000
burnin <- 40000
thin <- 10
YY <- exp(booth_barnett(seq(0,1,length=N)))
# Generate training data
exp_time <- c(2, 14, 20, 18, 7, 14, 13, 8, 16, 20, 3, 1, 10, 14, 12, 11, 14,
 10, 14, 20, 4, 5, 12, 19, 15)
set.seed(269535279)
sim <- generate_latent_gp_pois_known_func(booth_barnett, X, exp_time, raw=TRUE)
pgpfit <- thin.pgp(mcmc_pois(y=sim$y, X=X, T=t, nug_decay=TRUE, vb=TRUE),
 burnin=burnin, thin=thin)
pgppreds <- predict.pgp(fit=pgpfit, Xnew=XX, ret_all=TRUE)
ylims <- range(c(YY, sim$y$ymean, sim$yraw))
xlims <- range(XX)

# Plot true booth barnett function
pdf("booth_barnett_true.pdf", height=5.3, width=4.8)
plot(x=XX, y=YY, type="n", xlab="X", ylab="Y", xlim=xlims, ylim=ylims)
points(x=jitter(sim$Xraw[,1]), y=sim$yraw, pch=8, col=3)
points(x=X, y=sim$y$ymean)
lines(x=XX, y=YY, lwd=2)
legend("topright", c("poisson draws", "counts / exposure time", "truth"),
  pch=c(8, 1, NA), col=c(3, 1, 1), lwd=c(NA, NA, 2), lty=c(NA, NA, 1), bg="white")
dev.off()

pdf("booth_barnett_pgp_fit.pdf", height=5.3, width=4.8)
plot(x=XX, y=YY, type="n", xlab="X", ylab="Y", xlim=xlims, ylim=ylims)
matplot(x=XX, y=pgppreds$lambda_all[,seq(1, ncol(pgppreds$lambda_all), by=10)],
  type="l", col="lightgrey", lty=1, add=TRUE)
lines(x=XX, y=pgppreds$lambdap, lwd=2, col=2)
lines(x=XX, y=pgppreds$up, col=2, lty=2)
lines(x=XX, y=pgppreds$lp, col=2, lty=2)
points(x=jitter(sim$Xraw[,1]), y=sim$yraw, pch=8, col=3)
points(x=X, y=sim$y$ymean)
legend("topright", c("poisson gp fit", "MCMC predictive draws",
 "95% credible intervals"), pch=NA, col=c(2, "lightgrey", 2),
  lty=c(1, 1, 2))
dev.off()

pdgpfit <- thin.pdgp(mcmc_pois_deep(y=sim$y, X=X, T=t, g0=1e-6, nug_decay=TRUE, vb=TRUE),
 burnin=burnin, thin=thin)
pdgppreds <- predict.pdgp(fit=pdgpfit, Xnew=XX, ret_all=TRUE, ret_deep=TRUE)

pdf("booth_barnett_pdgp_fit.pdf", height=5.3, width=4.8)
plot(x=XX, y=YY, type="n", xlab="X", ylab="Y", xlim=xlims, ylim=ylims)
matplot(x=XX, y=pdgppreds$lambda_all[,seq(1, ncol(pdgppreds$lambda_all), by=10)],
  type="l", col="lightgrey", lty=1, add=TRUE)
lines(x=XX, y=pdgppreds$lambdap, lwd=2, col=4)
lines(x=XX, y=pdgppreds$up, col=4, lty=2)
lines(x=XX, y=pdgppreds$lp, col=4, lty=2)
points(x=jitter(sim$Xraw[,1]), y=sim$yraw, pch=8, col=3)
points(x=X, y=sim$y$ymean)
legend("topright", c("poisson deep gp fit", "MCMC predictive draws",
  "95% credible intervals"), pch=NA, col=c(4, "lightgrey", 4),
  lty=c(1, 1, 2))
dev.off()

pmwdgpfit <- thin.pdgp(mcmc_pois_deep(y=sim$y, X=X, T=t, monow=TRUE,
 aalign=TRUE, nug_decay=TRUE, sep=TRUE, vb=TRUE), burnin=burnin, thin=thin)
pmwdgppreds <- predict.pdgp(fit=pmwdgpfit, Xnew=XX, ret_all=TRUE,
  ret_deep=TRUE)

###############################################################################
###############################################################################
#### MONOTONIC WARPINGS SECTION
###############################################################################
###############################################################################

library(laGP)
library(mvtnorm)

set.seed(981734222)
X <- seq(0, 1, length=50)
DX <- distance(X)
ls <- 0.3
Sigma <- exp(-DX/ls) + sqrt(.Machine$double.eps)
Y1 <- rmvnorm(n=1, sigma=Sigma)
Y2 <- exp(Y1)
Y3 <- cumsum(Y2)
Y3_range <- range(Y3)
Y4 <- (Y3 - Y3_range[1])/diff(Y3_range)
y01_lims <- range(Y1, Y2)

pdf("monowarp_1.pdf", height=4.0, width=5.0)
par(mar=c(5.1, 5.1, 2.6, 0.2))
plot(x=X, y=Y1, xlab="", ylab="", main="")
title(ylab=expression(W[1]), line=2)
title(xlab="X", line=2)
dev.off()
pdf("monowarp_2.pdf", height=4.0, width=5.0)
par(mar=c(5.1, 5.1, 2.6, 0.2))
plot(x=X, y=Y2, xlab="", ylab="", main="")
title(ylab=expression(W[2]), line=2)
title(xlab="X", line=2)
dev.off()
pdf("monowarp_3.pdf", height=4.0, width=5.0)
par(mar=c(5.1, 5.1, 2.6, 0.2))
plot(x=X, y=Y3, xlab="", ylab="", main="")
title(ylab=expression(W[3]), line=2)
title(xlab="X", line=2)
dev.off()
pdf("monowarp_4.pdf", height=4.0, width=5.0)
par(mar=c(5.1, 5.1, 2.6, 0.2))
plot(x=X, y=Y4, xlab="", ylab="", main="")
title(ylab=expression(W[4]), line=2)
title(xlab="X", line=2)
dev.off()

pdf("booth_barnett_pmwdgp.pdf", height=5.3, width=4.8)
plot(x=XX, y=YY, type="n", xlab="", ylab="", xlim=xlims, ylim=ylims)
matplot(x=XX, y=pmwdgppreds$lambda_all[,seq(1, ncol(pmwdgppreds$lambda_all), by=10)],
  type="l", col="lightgrey", lty=1, add=TRUE)
lines(x=XX, y=pmwdgppreds$lambdap, lwd=2, col=7)
lines(x=XX, y=pmwdgppreds$up, col=7, lty=2)
lines(x=XX, y=pmwdgppreds$lp, col=7, lty=2)
points(x=jitter(sim$Xraw[,1]), y=sim$yraw, pch=8, col=3)
points(x=X, y=sim$y$ymean)
title(ylab="Y", line=2)
title(xlab="X", line=2)
legend("topright", c("poisson mono-warped deep gp fit", "95% credible intervals"),
  pch=NA, col=7, lty=c(1, 2), bg="white")
dev.off()


bd <- matrix(0.4615, ncol=1)
pdf("original_inputs.pdf", height=3.5, width=4.2)
par(mar=c(2.75, 2.75, 0.2, 0.2))
plot(x=XX, y=YY, type="n", xlab="", ylab="", xaxt="n", yaxt="n", xlim=xlims, ylim=ylims)
points(x=X, y=sim$y$ymean, cex=2.25)
points(x=X[c(3,5),], y=sim$y$ymean[c(3,5)], col=7, pch=19, cex=2.25)
points(x=X[c(17,19),], y=sim$y$ymean[c(17,19)], col=14, pch=19, cex=2.25)
text(x=X, y=sim$y$ymean, labels=1:nrow(X), cex=0.7)
abline(v=bd, col=5, lty=2, lwd=2)
title(ylab="Y", line=0.5)
title(xlab="X", line=0.5)
legend("topright", c("counts/exposure time (order in X)", "regime boundary"),
 pch=c(1,NA), lty=c(NA,2), col=c(1,5), cex=0.75, bg="white")
dev.off()

pdgp_Xpreds <- predict.pdgp(fit=pdgpfit, Xnew=rbind(X, bd), ret_deep=TRUE)
pdgp_nonmon_Xdraws <- c()
for (i in 1:dim(pdgp_Xpreds$wps)[3]) {
  inc <- pdgp_Xpreds$wps[2,1,i] > pdgp_Xpreds$wps[1,1,i]
  for (j in 2:(dim(pdgp_Xpreds$wps)[1] - 1)) {
    if (inc) {
       if (pdgp_Xpreds$wps[j,1,i] < pdgp_Xpreds$wps[j-1,1,i]) {
         pdgp_nonmon_Xdraws <- c(pdgp_nonmon_Xdraws, i)
         break
       }
    }
  }
}
rand_iter <- sample(pdgp_nonmon_Xdraws, 1)
wbd <- pdgp_Xpreds$wps[(nrow(X)+1),,rand_iter]

pdf("deep_gp_warped_inputs.pdf", height=3.5, width=4.2)
par(mar=c(2.75, 2.75, 0.2, 0.2))
plot(x=pdgp_Xpreds$wps[1:nrow(X),,rand_iter], y=sim$y$ymean, type="n",
 xaxt="n", yaxt="n", xlab="", ylab="", ylim=ylims)
points(x=pdgp_Xpreds$wps[1:nrow(X),,rand_iter], y=sim$y$ymean, cex=2.25)
points(x=pdgp_Xpreds$wps[1:nrow(X),,rand_iter][c(3,5)], y=sim$y$ymean[c(3,5)], col=7, pch=19, cex=2.25)
points(x=pdgp_Xpreds$wps[1:nrow(X),,rand_iter][c(17,19)], y=sim$y$ymean[c(17,19)], col=14, pch=19, cex=2.25)
text(x=pdgp_Xpreds$wps[1:nrow(X),,rand_iter], y=sim$y$ymean, labels=1:nrow(X), cex=0.7)
abline(v=wbd, col=5, lty=2, lwd=2)
title(xlab="W", line=0.75)
dev.off()

pmwdgp_Xpreds <- predict.pdgp(fit=pmwdgpfit, Xnew=rbind(X, bd), ret_deep=TRUE)
mwbd <- mean(pmwdgp_Xpreds$wps[(nrow(X)+1),,])
mw_mean <- apply(pmwdgp_Xpreds$wps[1:nrow(X),,], 1, mean)

pdf("mwdeep_gp_warped_inputs.pdf", height=3.5, width=4.2)
par(mar=c(2.75, 2.75, 0.2, 0.2))
plot(x=mw_mean, y=sim$y$ymean, type="n", xaxt="n", yaxt="n", xlab="",
 ylab="", ylim=ylims)
points(x=mw_mean, y=sim$y$ymean, cex=2.25)
points(x=mw_mean[c(3,5)], y=sim$y$ymean[c(3,5)], col=7, pch=19, cex=2.25)
points(x=mw_mean[c(17,19)], y=sim$y$ymean[c(17,19)], col=14, pch=19, cex=2.25)
text(x=mw_mean, y=sim$y$ymean, labels=1:nrow(X), cex=0.7)
title(xlab="W", line=0.75)
abline(v=mwbd, col=5, lty=2, lwd=2)
dev.off()

pdf("deep_gp_warpings.pdf", height=3.35, width=4.7)
par(mar=c(2.95, 3.25, 0.2, 0.2))
pdgp_nonmon_samp <- sample(pdgp_nonmon_Xdraws, 20)
matplot(x=XX, y=pdgppreds$wps[,1,seq(1, dim(pdgppreds$wps)[3], by=2)],
 type="l", col="lightgrey", lty=1, xlab="", ylab="")
matplot(x=XX, y=pdgppreds$wps[,1,pdgp_nonmon_samp], type="l",
  col=1:length(pdgp_nonmon_samp), lty=2, lwd=1, add=TRUE)
title(ylab="W", line=2.25)
title(xlab="X", line=2.0)
legend("topright", c("non-injective warpings"),
  pch=NA, lwd=2, col=1, lty=2, bg="white", cex=0.9)
dev.off()

pdf("deep_gp_monowarpings.pdf", height=3.35, width=4.7)
par(mar=c(2.95, 3.25, 0.2, 0.2))
matplot(x=XX, y=pmwdgppreds$wps[,1,seq(1, dim(pmwdgppreds$wps)[3], by=2)],
 type="l", col="lightgrey", lty=1, xlab="", ylab="")
title(ylab="W", line=2.25)
title(xlab="X", line=2.0)
dev.off()



###############################################################################
###############################################################################
#### DISCREPANCY ANALYSIS SECTION
###############################################################################
###############################################################################

source("../helper.R")
source('../vecchia_scaled.R')

md <- read.csv("../data/sims.csv")
ibex_real <- read.csv("../data/ibex_real.csv")
ibex_real$sim_counts <- ibex_real$counts
ibex_real$lat <- ibex_real$ecliptic_lat
ibex_real$lon <- ibex_real$ecliptic_lon
ibex_real <- ibex_real[ibex_real$map=="2020A" & ibex_real$esa==4,]

pd <- preprocess_data(md=md, fd=ibex_real, esa_lev=4,
  fparams=data.frame(year="2020A"), scales=c(1,1), tol=NA, quant=0.0,
  real=TRUE, disc=FALSE)
iter_fit <- fit_scaled(y=pd$Zm, inputs=as.matrix(cbind(pd$Xm, pd$Um)),
 nug=1e-4, ms=25)

est_pmfp <- (1696.111 - 500) / 2500 #
est_ratio <- (0.02355665 - 0.001) / (1-0.001)
XX_unit <- cbind(distinct(pd$Xm), matrix(c(est_pmfp, est_ratio), nrow=1))
colnames(XX_unit) <- c("x", "y", "z", "pmfp", "ratio")
lhat_pred <- predictions_scaled(iter_fit, as.matrix(XX_unit), m=25, joint=FALSE,
 predvar=FALSE)

obs_plot_data <- cbind(ibex_real[,c("lat", "lon")],
  ibex_real$sim_counts/ibex_real$time - ibex_real$background)
est_plot_data <- cbind(distinct(md[,c("lat", "lon")]), lhat_pred)
colnames(obs_plot_data) <- colnames(est_plot_data) <- c("lat", "lon", "rate")
obs_plot_data$nlon <- nose_center_lons(obs_plot_data$lon)
est_plot_data$nlon <- nose_center_lons(est_plot_data$lon)

real_sim_map <- read.csv("../data/ibex_sim_map_2020A.csv")
real_sim_map <- real_sim_map[order(real_sim_map$row),]
obs_plot_data$near_sim_rate <- obs_plot_data$pred_diff <- NA
for (i in 1:nrow(obs_plot_data)) {
  obs_plot_data$near_sim_rate[i] <- lhat_pred[real_sim_map$near_sim[i]]
}
obs_plot_data$pred_diff <- obs_plot_data$near_sim_rate - obs_plot_data$rate
max_diff <- max(abs(quantile(obs_plot_data$pred_diff, probs=c(0.05, 0.95))))

## Observed rates
pdf("obs_rates_2020.pdf", width=5.5, height=3.5)
ggplot(obs_plot_data, aes(x=nlon, y=lat)) +
  geom_point(aes(col=rate), size=0.36) +
  scale_color_gradientn(colours=ibex_palette$hex,
    limits=ena_range, oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  xlab("Longitude") + ylab("Latitude")+
  theme_bw()+
  theme(legend.position="right",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15),
        plot.background=element_rect(fill='transparent', color=NA),
        panel.background=element_rect(fill='transparent'))+
  labs(x="Longitude", y="Latitude", col="lambda") +
  ggtitle("Observed ENA Rate")
dev.off()

## Predicted simulator rates
pdf("est_rates_2020.pdf", width=5.5, height=3.5)
ggplot(est_plot_data, aes(x=nlon, y=lat)) +
  geom_raster(aes(fill=rate)) +
  scale_fill_gradientn(colours=ibex_palette$hex,
    limits=ena_range, oob=scales::oob_squish) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  xlab("Longitude") + ylab("Latitude")+
  theme_bw()+
  theme(legend.position="right",
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15),
        plot.background=element_rect(fill='transparent', color=NA),
        panel.background=element_rect(fill='transparent'),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="Longitude", y="", fill="ENAs/sec") +
  ggtitle("Simulator ENA Rate at Estimate")
dev.off()

## Discrepancy
pdf("discrepancy_2020.pdf", width=5.5, height=3.5)
ggplot(obs_plot_data, aes(x=nlon, y=lat)) +
  geom_point(aes(col=pred_diff), size=0.36) +
  scale_color_gradientn(colours=colorRampPalette(c("blue", "white", "red"))(128),
    limits=c(-max_diff, max_diff)) +
  scale_x_reverse(expand=c(0,0), breaks=new_360, labels=orig_360)+
  xlab("Longitude") + ylab("")+
  theme_bw()+
  theme(legend.position="right",
        legend.title=element_blank(),
        axis.title=element_text(size=15),
        plot.title=element_text(hjust=0.5, size=15),
        plot.background=element_rect(fill='transparent', color=NA),
        panel.background=element_rect(fill='transparent'),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())+
  labs(x="Longitude", y="") +
  ggtitle("Discrepancy (predicted - observed)")
dev.off()



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

pmfp=1325
ratio=0.001
 
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
