#### ACCRE simulation script for all adaptive twowave ODS designs
args = commandArgs(TRUE)
njob = as.integer(args[1])
wd = args[2]

# load libraries
library(MMLB)
library(MASS)
library(mitools)

setwd(wd)

# load functions
source("~/research/mphase_mm/aODS_fcts_RT20191010.R")

# twowave() arguments
DATA = NULL
SEEDS = 12345+5000*njob
YSEED = NULL
NLOW = 4
NHIGH = 6
N = 2500
PXE = 0.25
BETAZ = data.frame('Z'=1)
GAMMA = 2
SIGMA = NULL
OMEGAZ.Xe = 1  
BETAM = data.frame('Int'=-1.50, 'time'=-0.25, 'Xe'=0.25, 'Xe.i'=0.25)
LVMODEL = NULL
TMODEL = ~1
MODEL = 'Y~time+Xe+Xe.i+Z'
MARG.EXP.FORMULA = 'Xe~Z'
WAVE2MI = FALSE
MIREP = 10
GRIDREP = 1
GRIDINC = 20
INC = 10
RNGp = c(0.05, 0.95)
Q = 2
ITERLIM = 500
VERBOSE = FALSE

# Adaptive sample size (Xei)
adapt_sample_size = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N,
                                          PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe, 
                                          MARG.EXP.FORMULA=MARG.EXP.FORMULA, METHOD='ACML',
                                          N2FIXED=FALSE, WAVED1=c(0,150,0), WAVED2=350, TARGET='Xei',
                                          TARGETVAR=0.003, MIREP=MIREP, GRIDREP=GRIDREP,
                                          GRIDINC=GRIDINC, INC=INC, Q=Q, ITERLIM=ITERLIM,
                                          VERBOSE=VERBOSE, RNGp=RNGp))

if(!is.null(adapt_sample_size)) {
  saveRDS(adapt_sample_size, file=paste0('adaptive_sample_size_Xei_acml_', njob, '.RDS'))
}

adapt_sample_size = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N,
                                          PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe,
                                          BETAM=BETAM, LVMODEL=LVMODEL, TMODEL=TMODEL, MODEL=MODEL,
                                          MARG.EXP.FORMULA=MARG.EXP.FORMULA, METHOD='MI',
                                          N2FIXED=FALSE, WAVED1=c(0,150,0), WAVED2=350, TARGET='Xei',
                                          TARGETVAR=0.003, MIREP=MIREP, GRIDREP=GRIDREP,
                                          GRIDINC=GRIDINC, INC=INC, Q=Q, ITERLIM=ITERLIM,
                                          VERBOSE=VERBOSE, RNGp=RNGp))

if(!is.null(adapt_sample_size)) {
  saveRDS(adapt_sample_size, file=paste0('adaptive_sample_size_Xei_mi_', njob, '.RDS'))
}

# Adaptive design
### Xe
adapt_design_Xe = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N, 
                                        PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe,
                                        BETAM=BETAM, LVMODEL=LVMODEL, TMODEL=TMODEL, MODEL=MODEL,
                                        MARG.EXP.FORMULA=MARG.EXP.FORMULA, N2FIXED=TRUE,
                                        WAVED1=c(20,110,20), WAVED2=350, TARGET='Xe', TARGETVAR=NULL,
                                        MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, INC=INC,
                                        Q=Q, ITERLIM=ITERLIM, VERBOSE=VERBOSE, METHOD='ACML', RNGp=RNGp))

if(!is.null(adapt_design_Xe)) {
  saveRDS(adapt_design_Xe, file=paste0('adaptive_design_Xe_acml_', njob, '.RDS'))
}

adapt_design_Xe = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N,
                                        PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe,
                                        BETAM=BETAM, LVMODEL=LVMODEL, TMODEL=TMODEL, MODEL=MODEL,
                                        MARG.EXP.FORMULA=MARG.EXP.FORMULA, N2FIXED=TRUE,
                                        WAVED1=c(20,110,20), WAVED2=350, TARGET='Xe', TARGETVAR=NULL,
                                        MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, INC=INC,
                                        Q=Q, ITERLIM=ITERLIM, VERBOSE=VERBOSE, METHOD='MI', RNGp=RNGp))

if(!is.null(adapt_design_Xe)) {
  saveRDS(adapt_design_Xe, file=paste0('adaptive_design_Xe_mi_', njob, '.RDS'))
}


### XeXei
adapt_design_XeXei = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N,
                                           PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe,
                                           BETAM=BETAM, LVMODEL=LVMODEL, TMODEL=TMODEL, MODEL=MODEL,
                                           MARG.EXP.FORMULA=MARG.EXP.FORMULA, N2FIXED=TRUE,
                                           WAVED1=c(20,110,20), WAVED2=350, TARGET='XeXei', TARGETVAR=NULL,
                                           MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, INC=INC, 
                                           Q=Q, ITERLIM=ITERLIM, VERBOSE=VERBOSE, METHOD='ACML', RNGp=RNGp))

if(!is.null(adapt_design_XeXei)) {
  saveRDS(adapt_design_XeXei, file=paste0('adaptive_design_XeXei_acml_', njob, '.RDS'))
}

adapt_design_XeXei = tryCatch(twowave_aODS(DATA=DATA, SEEDS=SEEDS, YSEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N,
                                           PXE=PXE, BETAZ=BETAZ, GAMMA=GAMMA, SIGMA=SIGMA, OMEGAZ.Xe=OMEGAZ.Xe,
                                           BETAM=BETAM, LVMODEL=LVMODEL, TMODEL=TMODEL, MODEL=MODEL,
                                           MARG.EXP.FORMULA=MARG.EXP.FORMULA, N2FIXED=TRUE,
                                           WAVED1=c(20,110,20), WAVED2=350, TARGET='XeXei', TARGETVAR=NULL,
                                           MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, INC=INC,
                                           Q=Q, ITERLIM=ITERLIM, VERBOSE=VERBOSE, METHOD='MI', RNGp=RNGp))

if(!is.null(adapt_design_XeXei)) {
  saveRDS(adapt_design_XeXei, file=paste0('adaptive_design_XeXei_mi_', njob, '.RDS'))
}
