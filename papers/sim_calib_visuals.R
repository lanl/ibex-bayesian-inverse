source("helper.R")
source('../vecchia_scaled.R')

library(ggplot2)

## Visuals for varying the dimension of the response
res <- readRDS("sim_calib_results_20250902.rds")
single_index <- NA
single_pmfp <- 1750
single_ratio <- 0.02

for (i in 1:length(res)) {
  if (res[[i]]$truth[1]==single_pmfp && res[[i]]$truth[2]==single_ratio) {
    single_index <- i
    break
  }
}

pred_params <- apply(res[[single_index]]$u[seq(9001, 10000, by=10),], 2, mean)
pred_params[1] <- (pred_params[1] - 500)/2500
pred_params[2] <- (pred_params[2] - 0.001)/(0.1-0.001)

model_data <- read.csv(file="../../data/sims.csv")
field_data <- read.csv(file="../../data/sims_real.csv")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=c(single_pmfp, single_ratio), scales=c(1, 1), tol=NA, quant=0.0,
  real=FALSE, disc=FALSE)

field_data <- field_data[field_data$parallel_mean_free_path==single_pmfp &
  field_data$ratio==single_ratio,]
field_data$est_rate <- field_data$sim_counts/field_data$time - field_data$background
field_data$nlon <- nose_center_lons(field_data$lon)
field_data <- field_data[which(!is.nan(field_data$est_rate)),]
model_data <- model_data[model_data$parallel_mean_free_path==single_pmfp &
  model_data$ratio==single_ratio,]

fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)
XX_ll <- cbind(unique(model_data[,c("lon", "lat")]), matrix(pred_params, nrow=1))
XX_ll[,c("x", "y", "z")] <- geo_to_spher_coords(lat=XX_ll$lat, lon=XX_ll$lon)
XX_ll$x <- (XX_ll$x - min(XX_ll$x)) / diff(range(XX_ll$x))
XX_ll$y <- (XX_ll$y - min(XX_ll$y)) / diff(range(XX_ll$y))
XX_ll$z <- (XX_ll$z - min(XX_ll$z)) / diff(range(XX_ll$z))
XX <- XX_ll[,c("x", "y", "z", "1", "2")]
colnames(XX) <- c("x", "y", "z", "pmfp", "ratio")
lhat_curr <- predictions_scaled(fit, as.matrix(XX), m=25, joint=FALSE,
  predvar=FALSE)
pred_data <- data.frame(XX_ll, lhat_curr)
pred_data$nlon <- nose_center_lons(pred_data$lon)

# predrange <- range(c(model_data$blurred_ena_rate, field_data$est_rate), na.rm=TRUE)
predrange <- range(model_data$blurred_ena_rate, na.rm=TRUE)

ggplot(model_data, aes(x=ecliptic_lon_center, y=lat)) +
  geom_raster(aes(fill=blurred_ena_rate)) +
    scale_fill_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
      limits=predrange, oob=scales::oob_squish) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60),
      label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) +
    labs(x="", y="Latitude", fill="ena_rate")
ggsave("ibex_sim_mod.pdf", dpi=320, width=4, height=3.25)

ggplot(field_data, aes(x=nlon, y=lat)) +
  geom_point(aes(color=est_rate), size=0.4, na.rm=TRUE) +
  scale_colour_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
    limits=predrange, oob=scales::oob_squish) +
  scale_x_continuous(trans = "reverse",
    breaks = seq(325, 25, by = -60),
    label = c(60, 0, 300, 240, 180, 120)) +
  theme_bw() +
  theme(axis.ticks.y=element_blank(),
    axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0))) +
  labs(x="Longitude", y="Latitude", col="ena_rate")
ggsave("ibex_sim_field.pdf", dpi=320, width=4, height=3.25)

ggplot(pred_data, aes(x=nlon, y=lat)) +
  geom_raster(aes(fill=lhat_curr)) +
    scale_fill_gradientn(colours = rev(rainbow(6)[c(6,1:5)]),
      limits=predrange, oob=scales::oob_squish) +
    scale_x_continuous(trans = "reverse",
      breaks = seq(325, 25, by = -60),
      label = c(60, 0, 300, 240, 180, 120)) +
    theme_bw() +
    theme(axis.ticks.y=element_blank(),
      axis.title.x=element_text(margin=margin(t=10, r=0, b=0, l=0)),
      legend.spacing.y=unit(5.5, 'mm')) +
    labs(x="", y="Latitude", fill="ENA rate\n")
ggsave("ibex_sim_est.pdf", dpi=320, width=4, height=3.25)
