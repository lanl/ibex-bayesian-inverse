###############################################################################
#### Conducts Bayesian Computer Model Calibration for the data collected by
#### the Interstellar Boundary Explorer satellite. Uses field data from
#### 2009-2011.
###############################################################################

source("../helper.R")
source("../mcmc.R")

step_sizes <- c(0.005, 0.01, 0.02, 0.05)
logscl <- c(0, NA)
endyear <- c(2011, 2013)
beta <- c(1.1, 2, 4, 8)
design <- expand.grid(step_sizes, logscl, endyear, beta)
colnames(design) <- c("step", "scl", "year", "shape")
index <- 1
seed <- 7198701

## read in the command line arguments
## run with: R CMD BATCH '--args index=1' sun_cycle_calib.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)
print(design[index,])
settings <- design[index,]

model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"

fparams <- paste0(2009:settings$year, "A")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=fparams, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=10000, step=settings$step,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=settings$scl,
  betashape=settings$shape)
res <- list(mcmc_res=mcmc_res, settings)
saveRDS(res, file=paste0("../results/sun_cycle_index_", index, ".rds"))
