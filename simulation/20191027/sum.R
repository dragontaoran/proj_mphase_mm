method.set = c("ACML")
analysis.set = paste0("adaptive_", c("sample_size_Xei", "design_Xe", "design_XeXei"))
NSIM = 1000
out.dir = "results"
dir.create(out.dir, showWarnings=FALSE)

for (method in method.set) {
	prefix.set = paste0(analysis.set, "_", method)
	for (prefix in prefix.set) {
		res = list()
		for (nsim in 1:NSIM) {
			res[[nsim]] = readRDS(paste0("res/", prefix, "/", nsim, ".RDS"))
		}
		saveRDS(res, paste0(out.dir, "/", prefix, ".RDS"))
	}
}
