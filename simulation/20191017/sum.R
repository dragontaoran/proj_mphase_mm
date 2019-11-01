rm(list=ls())
gc()

NSIM = 1000

fn.in = paste("res/adaptive_sample_size_Xei_acml_", 1, ".RDS")
fi.adaptive_sample_size_Xei_acml = loadRDS(fn.in)

