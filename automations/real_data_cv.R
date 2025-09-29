###############################################################################
#### Cross-validation code for Bayesian Computer Model Calibration for data
#### collected by the Interstellar Boundary Explorer satellite.
#### Uses field data from 2009-2011. Results will be used to calculate CRPS on
#### held out data.
###############################################################################

source("../helper.R")
source("../mcmc.R")

seed <- 711930
ncvs <- 10
fid <- 1
nmcmcs <- 10000

## read in the command line arguments
## run with: R CMD BATCH '--args seed=711930 ncvs=10 fid=1 nmcmcs=10000' real_data_cv.R
args <- commandArgs(TRUE)
if (length(args) > 0) {
  for(i in 1:length(args)) {
    eval(parse(text=args[[i]]))
  }
}
print(args)

model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")
field_data <- field_data[field_data$map %in% paste0(2009:2011, "A"),]
nf <- nrow(field_data)
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"

set.seed(seed)
field_cv_inds <- sample(1:nf, nf)
field_cv_inds <- split(field_cv_inds, cut(seq_along(field_cv_inds), ncvs,
  labels=FALSE))
field_data <- field_data[-field_cv_inds[[fid]],]

fparams <- paste0(2009:2011, "A")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=fparams, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=nmcmcs, step=0.05,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=NA,
  betashape=2, adapt=FALSE)
res <- list(mcmc_res=mcmc_res, settings=list(seed=seed, ncvs=ncvs, fid=fid))
saveRDS(res, file=paste0("../results/real_data_cv_fold", fid, "_seed", seed,
  format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
