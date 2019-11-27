# Adaptive design

remove(list=ls())

# Load libraries (this is the location of the R packages on ACCRE)

.libPaths("~/R/R3.3.3lib/") ## where libraries stored.
#library(devtools)
#install_github('mercaldo/MMLB',force=TRUE)
library(MMLB)
library(MASS)

# Update data, result, tmp paths (dpath, rpath, tpath)
# The results and tmp directories are automatically created using the .slurm file

path      <- 'rsch/OutDepSamp/Projects/AODSBinary/LHSAnalysis/ACCRE'
dir_a_des <- 'adaptive_ss'
dir_a_ss  <- 'adaptive_des'
dir_f     <- 'fixed2stage'
dir_lhs   <- 'lhs'

dpath       <- file.path('~', path) # location of aODS_functions.r, lhs.Rdata files

rpath_a_des <- file.path('~', path, dir_a_des, 'results', 'seed_bbb')
tpath_a_des <- file.path('~', path, dir_a_des, 'tmp', 'seed_bbb')

rpath_a_ss  <- file.path('~', path, dir_a_ss, 'results', 'seed_bbb')
tpath_a_ss  <- file.path('~', path, dir_a_ss, 'tmp', 'seed_bbb')

rpath_f     <- file.path('~', path, dir_f, 'results', 'seed_bbb')
tpath_f     <- file.path('~', path, dir_f, 'tmp', 'seed_bbb')

rpath_lhs   <- file.path('~', path, dir_lhs, 'results/')
dir.create(rpath_lhs)

source(paste(dpath,'aODS_functions.R',sep='/'))

LHS <- FALSE # Use code to perform LHS analysis or simulation

# Arguments
SEEDS     <- bbb          # seed updated in slurm file
ITERLIM   <- 500          # 
VERBOSE   <- TRUE         # print output
Q         <- 30           # number of quadrature points for mLV, mTLV models

# Fixed two-stage design arguments
SAMP    <- 500  # fixed 2 stage sample size
SS_DES  <- list(c(0,500,0), c(25,450,25), c(50,400,50),c(75,350,75), c(100,300,100), c(125,250,125), c(150,200,150)  )
TS_DES1 <- list(c(0,100,0), c(25,50,25))
TS_DES2 <- list(c(0,400,0), c(25,350,25), c(50,300,50), c(75,250,75), c(100,200,100), c(125,150,125) ) 

# Adaptive two-stage design arguments
PROFILEEV <- 0.01         # kappa used for adaptive sample size
S1D       <- c(25,50,25)  # stage 1 design
BREP      <- 50           # number of imputed data sets when estimating D2[0,N2,0] or n2*D2[a,2a,a]
MOD       <- 10           # increments of balanced design: D2[alpha, 2alpha, alpha]
N2        <- 400          # stage 2 sample size (fixed sample size)

#PVAL_XeXei <- c(0,0)     # Ask J about this, currently we are plugging in the full cohort (true) estimates of these values
#PVAL_Xe    <- 0          # Profiling not performed! I think it makes sense to change the estimates of Xe, Xe.i=0.     
#PVAL_Xei   <- 0
#   
# if(!LHS) {
#   # Simulation arguments
#   YSEED     <- NULL         # seed if we want to fix Y
#   NLOW      <- 3            # Min number of repeated measurements
#   NHIGH     <- 5            # Max number of repeated measurements 
#   N         <- 5000         # Sample size
#   PXE       <- .25          # P(Xe=1)
#   BETAZ     <- 1            # binary confounder
#   BETAM     <- c(-1.50,-0.25,-.2,.2)
#   GAMMA     <- 2            # 
#   SIGMA     <- NULL         # 
#   LVMODEL   <- NULL         # lv model
#   TMODEL    <- ~1           # transition model
#   lv_mod    <- LVMODEL      # lv model 
#   t_mod     <- TMODEL       # transition model
#   MODEL     <- Y~time+Xe+Xe.i+Z
#   Y_mod     <- MODEL        # y model
#   me_mod    <- Xe~Z         # marginal exposure model - see ImputeData function
# 
#   ########################
#   # Adaptive sample size #
#   ########################
#   
#   aODS_ss_sim_g_gt <- tryCatch( adaptive_ss(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                             PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                             GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                             MODEL = MODEL, PMODEL = "Y~time+Z", Q = Q, ITERLIM = ITERLIM, 
#                                             PROFILEVAR = c("Xe","Xe.i"), VERBOSE = VERBOSE, S1D = S1D, 
#                                             BREP = BREP, PROFILEEV=PROFILEEV))
#   
#   if(!is.null(aODS_ss_sim_g_gt)) saveRDS(aODS_ss_sim_g_gt, file=paste(rpath_a_ss,paste('sim_g_gt_',bbb,'.RDS',sep=''),sep='/'))
#   
#   aODS_ss_sim_g <- tryCatch( adaptive_ss(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                          PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, GAMMA = GAMMA, 
#                                          SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                          MODEL = MODEL, PMODEL = "Y~time+Z", Q = Q, ITERLIM = ITERLIM, 
#                                          PROFILEVAR = c("Xe"), VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, 
#                                          PROFILEEV=PROFILEEV))
#   
#   if(!is.null(aODS_ss_sim_g)) saveRDS(aODS_ss_sim_g, file=paste(rpath_a_ss,paste('sim_g_',bbb,'.RDS',sep=''),sep='/'))
#   
#   aODS_ss_sim_gt <- tryCatch( adaptive_ss(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                           PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                           GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                           MODEL = MODEL, PMODEL = "Y~time+Z", Q = Q, ITERLIM = ITERLIM, 
#                                           PROFILEVAR = c("Xe.i"), VERBOSE = VERBOSE, S1D = S1D, BREP = BREP,
#                                           PROFILEEV=PROFILEEV))
#   
#   if(!is.null(aODS_ss_sim_gt)) saveRDS(aODS_ss_sim_gt, file=paste(rpath_a_ss,paste('sim_gt_',bbb,'.RDS',sep=''),sep='/'))
#   
#   ###################
#   # Adaptive design #
#   ###################
#   
#   aODS_des_sim_g_gt <- tryCatch( adaptive_des(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                               PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                               GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                               MODEL = MODEL, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe","Xe.i"), 
#                                               VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD, MARG.EXP.FORMULA=me_mod))
#   
#   if(!is.null(aODS_des_sim_g_gt)) saveRDS(aODS_des_sim_g_gt, file=paste(rpath_a_des,paste('sim_g_gt_',bbb,'.RDS',sep=''),sep='/'))
#   
#   aODS_des_sim_g <- tryCatch( adaptive_des(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                            PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                            GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                            MODEL = MODEL, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe"),
#                                            VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD, MARG.EXP.FORMULA=me_mod))
#   
#   if(!is.null(aODS_des_sim_g)) saveRDS(aODS_des_sim_g, file=paste(rpath_a_des,paste('sim_g_',bbb,'.RDS',sep=''),sep='/'))
#   
#   aODS_des_sim_gt <- tryCatch( adaptive_des(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                             PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                             GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                             MODEL = MODEL, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe.i"), 
#                                             VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD, MARG.EXP.FORMULA=me_mod))
#   
#   if(!is.null(aODS_des_sim_gt)) saveRDS(aODS_des_sim_gt, file=paste(rpath_a_des,paste('sim_gt',bbb,'.RDS',sep=''),sep='/'))
#   
#   #################
#   # Fixed designs #
#   #################
#   
#   f2s <- tryCatch( fixed2stage(SEEDS = SEEDS, YSEED = YSEED, NLOW = NLOW, NHIGH = NHIGH, N = N, 
#                                PXE = PXE, BETAZ = BETAZ, BETAM = BETAM, 
#                                GAMMA = GAMMA, SIGMA = SIGMA, LVMODEL = LVMODEL, TMODEL = TMODEL,
#                                MODEL = MODEL, Q = Q, ITERLIM = ITERLIM, 
#                                VERBOSE = VERBOSE, SAMP=SAMP, SS_DES=SS_DES, TS_DES1=TS_DES1, TS_DES2=TS_DES2))
#   
#   if(!is.null(f2s)) saveRDS(f2s, file=paste(rpath_f,paste('sim_',bbb,'.RDS',sep=''),sep='/'))
#   
#   
#   ## Remove tmp directories
#   system( paste( 'rm -r ', tpath_a_ss, sep='') )
#   system( paste( 'rm -r ', tpath_a_des, sep='') )
#   system( paste( 'rm -r ', tpath_f, sep='') ) 
#   
# } else { 
#   
  # Arguments
  Y_mod  <- Y~fevpct0.2+cigs0.10+pack0.20+smk.b+smk.w+time+female+bmi0.5+age0.10+Xe+Xe.i#+factor(site)
  lv_mod <- ~1
  t_mod  <- ~1
  marg.mean.exp.formula <- Xe ~ pack0.20+cigs0.10+bmi0.5+age0.10+smk.b+female#+factor(site)
 
  ades     <- c(25,50,25)
  S1D      <- ades
  kappa    <- 0.04^2
  ades_ss  <- c(25,50,25) 
  N2       <- 400 # number to sample when performing an adaptive ODS designs
  skip_int <- 10
  MOD      <- 10
  M        <- 25
  
  load(paste(dpath,'lhs.RData',sep='/')) # load lhs object
 # load("/Users/nmercaldo/Dropbox (Partners HealthCare)/Nate/aODS/LHS/lhs.RData")
  
  dat     <- lhs
  dat$id0 <- dat$id
  dat$id  <- as.numeric(factor(dat$id)) # convert to 1:N
  
  dat$visit1           <- dat$visit-1
  dat$fevpct_diff      <- dat$fevpct-dat$fevpct0
  dat$fevpct_abschange <- dat$fevpct_diff/dat$visit # visit is already post baseline
  dat$fevpct0.2        <- (dat$fevpct0-75)/2
  
  dat$Y <- (dat$fevpct_abschange < -2)*1 
  
  dat$age0         <- dat$age
  dat$age0.10      <- (dat$age0-50)/10   
  dat$cigs0.10     <- (dat$f31cigs-30)/10
  dat$pack0.20     <- (dat$packyear-40)/20
  dat$bmi5         <- (dat$bmi-26)/5
  dat$bmi0.5       <- (dat$bmi0-26)/5
  dat$fvc0.1       <- (dat$fvc0-5)
  dat$fevfvcrat0.1 <- (dat$fevfvcrat0-0.65)
  
  dat$smk.b     <- cluster.summary(dat$id, dat$smk, mean)
  dat$smk.w     <- dat$smk - dat$smk.b
  
  # Add recessive model for snp (see analysis.R to see how this snp was identified; models did NOT account for site differences)
  dat$snp <- (dat$fev.rs177852!='CC')*1
  
  dat$Xe   <- dat$snp
  dat$time <- dat$visit1
  dat$Xe.i <- dat$Xe * dat$time
  
  adat <- dat[dat$nobs>3,]
  adat <- adat[order(adat$id, adat$visit), ]
  
  ss <- tapply(adat$Y, adat$id, mean)*2
  ss[ss<2 & ss>0] <- 1
  adat <- merge(adat, data.frame(ss, id=names(ss)), by='id')
  #adat <- adat[order(adat$id, adat$visit), ]
  adat <- adat[order(adat$id, adat$time), ]
  adat$ss <- adat$ss + 1
  
  # Full cohort analysis
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing full cohort analysis')
  fc.tlv <- mm(Y_mod, dat=adat, id=id, lv.formula=lv_mod,  t.formula=t_mod, verbose=VERBOSE, iter.lim=ITERLIM, q=Q)
  fc_ests <- cbind( unlist(coef(fc.tlv)), sqrt(diag(fc.tlv$mod.cov)))
  initial.vals <- c(fc.tlv$beta,fc.tlv$alpha)
  
  saveRDS(fc_ests, file=paste(rpath_lhs,paste('fc_lhs_',bbb,'.RDS',sep=''),sep='/'))
  
  # Fixed two stage designs
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing two stage fixed ODS design(s)')
  lhs_f2s <- tryCatch( fixed2stage(DATA=adat, SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                               MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, 
                               VERBOSE = VERBOSE, SAMP=SAMP, SS_DES=SS_DES, TS_DES1=TS_DES1, TS_DES2=TS_DES2))
  
  saveRDS(lhs_f2s, file=paste(rpath_lhs,paste('f2s_lhs_',bbb,'.RDS',sep=''),sep='/'))
  
  # Adpative design 
  lhs_aODS_des_g_gt <- tryCatch( adaptive_des(DATA=adat,SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                              MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe","Xe.i"), 
                                              VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD,
                                              MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_des_g_gt)) saveRDS(lhs_aODS_des_g_gt, file=paste(rpath_lhs,paste('lhs_aODS_des_g_gt_',bbb,'.RDS',sep=''),sep='/'))
  
  lhs_aODS_des_g <- tryCatch( adaptive_des(DATA=adat, SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                           MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe"),
                                           VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD,
                                           MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_des_g)) saveRDS(lhs_aODS_des_g, file=paste(rpath_lhs,paste('lhs_aODS_des_g_',bbb,'.RDS',sep=''),sep='/'))
  
  lhs_aODS_des_gt <- tryCatch( adaptive_des(DATA=adat, SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                            MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe.i"), 
                                            VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, N2 = N2, MOD=MOD,
                                            MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_des_gt)) saveRDS(lhs_aODS_des_gt, file=paste(rpath_lhs,paste('lhs_aODS_des_gt_',bbb,'.RDS',sep=''),sep='/'))
  
  # Adaptive sample size
  lhs_aODS_ss_g_gt <- tryCatch( adaptive_ss(DATA=adat,SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                              MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe","Xe.i"), 
                                              VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, PROFILEEV = PROFILEEV,
                                              MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_ss_g_gt)) saveRDS(lhs_aODS_ss_g_gt, file=paste(rpath_lhs,paste('lhs_aODS_ss_g_gt_',bbb,'.RDS',sep=''),sep='/'))
  
  lhs_aODS_ss_g <- tryCatch( adaptive_ss(DATA=adat, SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                           MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe"),
                                           VERBOSE = VERBOSE, S1D = S1D, BREP = BREP, PROFILEEV = PROFILEEV,
                                           MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_ss_g)) saveRDS(lhs_aODS_ss_g, file=paste(rpath_lhs,paste('lhs_aODS_ss_g_',bbb,'.RDS',sep=''),sep='/'))
  
  lhs_aODS_ss_gt <- tryCatch( adaptive_ss(DATA=adat, SEEDS = SEEDS, LVMODEL = lv_mod, TMODEL = t_mod,
                                            MODEL = Y_mod, Q = Q, ITERLIM = ITERLIM, PROFILEVAR = c("Xe.i"), 
                                            VERBOSE = VERBOSE, S1D = S1D, BREP = 5, PROFILEEV = PROFILEEV,
                                            MARG.EXP.FORMULA=marg.mean.exp.formula))
  
  if(!is.null(lhs_aODS_ss_gt)) saveRDS(lhs_aODS_ss_gt, file=paste(rpath_lhs,paste('lhs_aODS_ss_gt_',bbb,'.RDS',sep=''),sep='/'))
# }


