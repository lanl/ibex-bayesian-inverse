###############################################################################
#### Code for testing recovery of a multiplicative scale discrepancy in the
#### the context of a Bayesian inverse problem for counts (i.e. Bayesian
#### Computer Model Calibration). Creates synthetic satellite data with an
#### artifical scaling factor. Satellite observation locations and exposure
#### times come from a year specified by the user that the Interstellar
#### Boundary Explorer has been in operation. Rate of ENAs at these
#### locations is determined by the IBEX simulator.
####
#### DATA NEEDED: sims.csv, ibex_real.csv
###############################################################################

setwd("..")
source("pois_bayes_inv.R")
source("helper.R")
setwd("automations")

index <- NA
scale <- NA
seed <- 7198701
map <- 2020

scales <- c(seq(0.2, 1.0, by=0.2), 1.5, 2, 4, 8, 10)

## read in the command line arguments
## run with: R CMD BATCH '--args index=1 scale=0.2 seed=719 map=2020' scale_disc_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}

if (is.na(index) && is.na(scale)) {
  stop("must provide a scale or an index")
} else if (!is.na(index) && is.na(scale)) {
  scale <- scales[index]
}
map <- paste0(map, "A")

settings <- data.frame(index=index, scale=scale, seed=seed, map=map)
print(settings)

## Load model and field data
model_data <- read.csv(file="../data/sims.csv")
pmfps <- unique(model_data$parallel_mean_free_path)
ratios <- unique(model_data$ratio)
u <- matrix(c(1750, 0.02), nrow=1)
u_unit <- matrix(c((u[1]-min(pmfps))/diff(range(pmfps)),
  (u[2]-min(ratios))/diff(range(ratios))), nrow=1)
field_data <- read.csv(file="../data/ibex_real.csv")
field_data <- field_data[field_data$map==map,]

## Fit Scaled Vecchia GP surrogate
pd <- preprocess_data(md=model_data, fd=field_data, map=map, esa_lev=4)
fit <- fit_scaled(y=pd$ym, inputs=cbind(pd$xm, pd$um), nug=1e-4, ms=25)

## Predict at field locations
XX <- cbind(pd$xf, matrix(u_unit, nrow=nrow(pd$xf), ncol=length(u_unit), byrow=TRUE))
field_data$ena_rate <- predictions_scaled(fit, XX, m=25, joint=FALSE)
field_data$counts <-
  rpois((field_data$ena_rate*scale+field_data$background)*field_data$time,
    (field_data$ena_rate*scale+field_data$background)*field_data$time)
pd <- preprocess_data(md=model_data, fd=field_data, map=map, esa_lev=4)

## Create simulated field data with multiplicative discrepancy
res <- pois_bayes_inv(xm=pd$xm, um=pd$um, ym=pd$ym, xf=pd$xf, yf=pd$yf, e=pd$e,
  lam0=pd$bg, T=10000, vb=TRUE, settings=list(sample_scl=TRUE))
saveRDS(list(res=res, field_data=field_data, settings=settings),
  file=paste0("scale_disc_test_scale_", scale,
    format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
