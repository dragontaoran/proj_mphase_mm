method.set = c("MI", "ACML")
analysis.set = paste0("adaptive_", c("sample_size_Xei"))
NSIM = 500
out.dir = "results"
dir.create(out.dir, showWarnings=FALSE)

for (method in method.set) {
	prefix.set = paste0(analysis.set, "_", method)
	for (prefix in prefix.set) {
		res = list()
		j = 1
		for (nsim in 1:NSIM) {
			fn = paste0("res/", prefix, "/", nsim, ".RDS")
			if (file.exists(fn)) {
				res[[j]] = readRDS(fn)
				j = j+1
			}
		}
		saveRDS(res, paste0(out.dir, "/", prefix, ".RDS"))
	}
}
