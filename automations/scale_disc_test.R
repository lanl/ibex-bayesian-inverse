###############################################################################
#### Conducts Bayesian Computer Model Calibration for the data collected by
#### the Interstellar Boundary Explorer satellite. Uses field data from
#### 2009-2011.
###############################################################################

source("../helper.R")
source("../mcmc.R")
source("../vecchia_scaled.R")

index <- NA
scale <- NA
seed <- 7198701
map <- "2020A"

scales <- c(seq(0.2, 1.0, by=0.2), 1.5, 2, 4, 8, 10)

## read in the command line arguments
## run with: R CMD BATCH '--args index=1 scale=0.2 seed=719 map="2020A"' scale_disc_test.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)

if (is.na(index) && is.na(scale)) {
  stop("must provide a scale or an index")
} else if (!is.na(index) && is.na(scale)) {
  scale <- scales[index]
}
settings <- data.frame(index=index, scale=scale, seed=seed)
print(settings)

## Load model and field data
model_data <- read.csv(file="../data/sims.csv")
u <- matrix(c(1750, 0.02), nrow=1)
u_unit <- matrix(c((u[1]-500)/2500, (u[2]-0.001)/(0.1-0.001)), nrow=1)
field_data <- read.csv(file="../data/ibex_real.csv")
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"
field_data <- field_data[field_data$map==map,]

## Fit Scaled Vecchia GP surrogate
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=map, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
fit <- fit_scaled(y=pd$Zmod, inputs=as.matrix(cbind(pd$Xmod, pd$Umod)),
 nug=1e-4, ms=25)

## Predict at field locations
XX <- cbind(pd$Xfield, u_unit)
colnames(XX) <- NULL
field_data$ena_rate <- predictions_scaled(fit, as.matrix(XX), m=25,
  joint=FALSE, predvar=FALSE)
field_data$sim_counts <-
  rpois((field_data$ena_rate*scale+field_data$background)*field_data$time,
    (field_data$ena_rate*scale+field_data$background)*field_data$time)
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=map, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)

## Create simulated field data with multiplicative discrepancy
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=10000, step=0.05,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=NA, betashape=2)
res <- list(mcmc_res=mcmc_res, field_data=field_data, settings=settings)
saveRDS(res, file=paste0("results/scale_disc_test_scale_", scale,
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
