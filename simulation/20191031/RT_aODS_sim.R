#### ACCRE simulation script for all adaptive twowave ODS designs
args = commandArgs(TRUE)
njob = as.integer(args[1])
wd = args[2]
METHOD = args[3]

# load functions
source("RT_aODS_fcts.R")

# twowave() arguments
SEEDS = 12345+5000*njob
NLOW = 4
NHIGH = 6
N = 2500
PXE = 0.25
BETAZ = data.frame('Z'=1)
GAMMA = 2
OMEGAZ.Xe = 1  
BETAM = data.frame('Int'=-1.50, 'time'=-0.25, 'Xe'=0.25, 'Xe.i'=0.25)
MIREP = 10
GRIDREP = 1
GRIDINC = 20
INC = 10
RNGp = c(0.05, 0.95)
BREP = 10
MIREP.final = 50

if (METHOD == "ACML") {
	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, 
		BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM, 
		METHOD=METHOD, N2FIXED=TRUE, WAVED1=c(20, 110, 20), WAVED2=350, 
		TARGET='Xe', TARGETVAR=NULL, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, 
		INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))
	if (!is.null(res)) {
		subwd = paste0(wd, "/adaptive_design_Xe_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}



	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, 
		BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM, 
		METHOD=METHOD, N2FIXED=TRUE, WAVED1=c(20, 110, 20), WAVED2=350,
		TARGET='XeXei', TARGETVAR=NULL, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC,
		INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))
	if (!is.null(res)) {
		subwd = paste0(wd, "/adaptive_design_XeXei_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}



	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, 
							  BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM,
							  METHOD=METHOD, N2FIXED=FALSE, WAVED1=c(0, 150, 0), WAVED2=350, 
							  TARGET='Xei', TARGETVAR=0.003, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC, 
							  INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))
	if (!is.null(res)) {
		subwd = paste0(wd, "/adaptive_sample_size_Xei_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}
  
} else if (METHOD == "MI") {
  
  
  
	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, 
		BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM,
		METHOD=METHOD, N2FIXED=TRUE, WAVED1=c(20, 110, 20), WAVED2=350, 
		TARGET='Xe', TARGETVAR=NULL, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC,
		INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))
	if (!is.null(res)) {
		subwd = paste0(wd, "/adaptive_design_Xe_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}



	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE,
		BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM, 
		METHOD=METHOD, N2FIXED=TRUE, WAVED1=c(20, 110, 20), WAVED2=350,
		TARGET='XeXei', TARGETVAR=NULL, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC,
		INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))

	if (!is.null(res)) {
		subwd = paste0(wd, "/adaptive_design_XeXei_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}



	res = tryCatch(twowave_aODS(SEEDS=SEEDS, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE,
		BETAZ=BETAZ, GAMMA=GAMMA, OMEGAZ.Xe=OMEGAZ.Xe, BETAM=BETAM,
		METHOD=METHOD, N2FIXED=FALSE, WAVED1=c(0, 150, 0), WAVED2=350, 
		TARGET='Xei', TARGETVAR=0.003, MIREP=MIREP, GRIDREP=GRIDREP, GRIDINC=GRIDINC,
		INC=INC, RNGp=RNGp, BREP=BREP, MIREP.final=MIREP.final))
	if(!is.null(adapt_sample_size)) {
		subwd = paste0(wd, "/adaptive_sample_size_Xei_", METHOD)
		dir.create(subwd, showWarnings=FALSE)
		saveRDS(res, file=paste0(subwd, '/', njob, '.RDS'))
	}  
}
