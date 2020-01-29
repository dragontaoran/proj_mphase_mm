#### Author: Ran Tao (r.tao@vanderbilt.edu) ####
#### Date: 11/22/2019						####
rm(list = ls())
gc()

#### ACCRE simulation script for all adaptive twowave ODS designs
args = commandArgs(TRUE)
njob = as.integer(args[1])
wd = args[2]
METHOD = args[3]
goal = args[4]


snp = "rs7911302"

# load functions
source("RT_aODS_fcts.R")

# load and process LHS data
load("/home/taor2/research/mphase_mm/lhs/lhsMWODSBinary.Rdata")
dat3 = dat2[,c("id", "cigs0.10", "pack0.20", "female", "age.10", "bmi0.5", "bmi5.change")]
dat3$time = dat2$visit-1
dat3$Xe = dat2[,snp]
dat3$Xe.i = dat3$Xe*dat3$time
dat3$Y = dat2$fevpctdec5
dat3 = dat3[which(!is.na(dat3$Y)),] # The outcome and covariates cannot be missing.
dat3_s = lapply(split(dat3, dat3$id), function(ZZ) {
    nY = sum(ZZ$Y)
    nn = nrow(ZZ)
    ZZ$ss = (nY==0)*1+(nY>0 & nY<nn)*2+(nY==nn)*3
    ZZ
})
dat3 = do.call(rbind, dat3_s)

# # check if all the covariates are time-invariant
# ids= unique(dat3$id)
# for (id in ids) {
#     idx = which(dat3$id == id)
#     if (length(unique(dat3[idx, "cigs0.10"])) > 1) print(paste0("cigs0.10 is time-varying for id ", id))
#     if (length(unique(dat3[idx, "pack0.20"])) > 1) print(paste0("pack0.20 is time-varying for id ", id))
#     if (length(unique(dat3[idx, "female"])) > 1) print(paste0("female is time-varying for id ", id))
#     if (length(unique(dat3[idx, "age.10"])) > 1) print(paste0("age.10 is time-varying for id ", id))
#     if (length(unique(dat3[idx, "bmi0.5"])) > 1) print(paste0("bmi0.5 is time-varying for id ", id))
#     if (length(unique(dat3[idx, "bmi5.change"])) > 1) print(paste0("bmi5.change is time-varying for id ", id))
# }
# # bmi5.change is time-varying

# twowave() arguments
SEEDS = 12345+5000*njob
LVMODEL = ~1
TMODEL = ~1
MODEL = "Y~time+Xe+Xe.i+cigs0.10+pack0.20+female+age.10+bmi0.5+bmi5.change"
MARG.EXP.FORMULA = 'Xe~cigs0.10+pack0.20+female+age.10+bmi0.5'
MIREP = 5
GRIDREP = 1
GRIDINC = 20
INC = 10
RNGp = c(0.05, 0.95)
RNGp.xei = c(0.05, 1)
BREP = 5
MIREP.final = 50
Q = 10
WAVED1 = c(10, 120, 20)
WAVED1.xei = c(0, 150, 0)
WAVED2 = 350

if (goal == 1) {
    res = tryCatch(twowave_aODS(DATA = dat3, SEEDS = SEEDS, LVMODEL = LVMODEL, TMODEL = TMODEL,
                                MODEL = MODEL, MARG.EXP.FORMULA = MARG.EXP.FORMULA, METHOD = METHOD, 
                                N2FIXED = TRUE, WAVED1 = WAVED1, WAVED2 = WAVED2, 
                                TARGET = 'Xe', TARGETVAR = NULL, MIREP = MIREP, GRIDREP = GRIDREP, GRIDINC = GRIDINC, 
                                INC = INC, RNGp = RNGp, Q = Q, BREP = BREP, MIREP.final = MIREP.final))
    if (!is.null(res)) {
        subwd = paste0(wd, "/adaptive_design_Xe_", METHOD)
        dir.create(subwd, showWarnings=FALSE)
        saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
    }   
} else if (goal == 2) {
    res = tryCatch(twowave_aODS(DATA = dat3, SEEDS = SEEDS, LVMODEL = LVMODEL, TMODEL = TMODEL,
                                MODEL = MODEL, MARG.EXP.FORMULA = MARG.EXP.FORMULA, METHOD = METHOD, 
                                N2FIXED = TRUE, WAVED1 = WAVED1, WAVED2 = WAVED2, 
                                TARGET = 'XeXei', TARGETVAR = NULL, MIREP = MIREP, GRIDREP = GRIDREP, GRIDINC = GRIDINC, 
                                INC = INC, RNGp = RNGp, Q = Q, BREP = BREP, MIREP.final = MIREP.final))
    if (!is.null(res)) {
        subwd = paste0(wd, "/adaptive_design_XeXei_", METHOD)
        dir.create(subwd, showWarnings=FALSE)
        saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
    }
} else if (goal == 3) {
    res = tryCatch(twowave_aODS(DATA = dat3, SEEDS = SEEDS, LVMODEL = LVMODEL, TMODEL = TMODEL,
                                MODEL = MODEL, MARG.EXP.FORMULA = MARG.EXP.FORMULA, METHOD = METHOD, 
                                N2FIXED = FALSE, WAVED1 = WAVED1.xei, WAVED2 = WAVED2, 
                                TARGET = 'Xei', TARGETVAR = 0.0025, MIREP = MIREP, GRIDREP = GRIDREP, GRIDINC = GRIDINC, 
                                INC = INC, RNGp = RNGp.xei, Q = Q, BREP = BREP, MIREP.final = MIREP.final))
    if (!is.null(res)) {
        subwd = paste0(wd, "/adaptive_sample_size_Xei_", METHOD)
        dir.create(subwd, showWarnings=FALSE)
        saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
    }    
}
