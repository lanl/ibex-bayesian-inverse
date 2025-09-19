source("../helper.R")
source("../mcmc.R")

###############################################################################
#### Conducts Bayesian Computer Model Calibration for the data collected by
#### the Interstellar Boundary Explorer satellite. Uses field data from
#### 2009-2011.
###############################################################################

model_data <- read.csv(file="../data/sims.csv")
field_data <- read.csv(file="../data/ibex_real.csv")
colnames(field_data)[colnames(field_data)== "counts"] <- "sim_counts"
colnames(field_data)[colnames(field_data)== "ecliptic_lat"] <- "lat"
colnames(field_data)[colnames(field_data)== "ecliptic_lon"] <- "lon"

fparams <- paste0(2009:2011, "A")
pd <- preprocess_data(md=model_data, fd=field_data, esa_lev=4,
  fparams=fparams, scales=c(1, 1), tol=NA, quant=0.0, real=TRUE)
mcmc_res <- mcmc(Xm=pd$Xmod, Um=pd$Umod, Zm=pd$Zmod, Xf=pd$Xfield,
  Zf=pd$Zfield, Of=pd$Ofield, m=25, nmcmcs=30000, step=0.06517432,
  gpmeth="svecchia", vb=TRUE, true_u=NA, true_logscl=NA)
res <- list(mcmc_res=mcmc_res, esa=4, nmcmcs=30000,
 params=fparams, real=TRUE, step=0.06517432)
saveRDS(res, file=paste0("../results/mcmc_res_esa", 4, "_nmcmc", 30000,
  paste0("_real", paste0(fparams, collapse="")), format(Sys.time(), "_%Y%m%d%H%M%S"), ".rds"))
