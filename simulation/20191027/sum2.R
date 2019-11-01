variable = c("Intercept", "Time", "Xe", "Xe.i", "Z", "Gamma")
true.value = c(-1.5, -0.25, 0.25, 0.25, 1, 2)
method.set = c("ACML")
analysis.set = paste0("adaptive_", c("sample_size_Xei", "design_Xe", "design_XeXei"))
q975 = qnorm(0.975)

est_extract = function (res) {
    c(res$ests$coef$beta, res$ests$coef$alpha)
}

se_extract = function (res) {
    c(sqrt(diag(res$ests$vcov$beta)), sqrt(res$ests$vcov$alpha))
}

design_extract = function (res) {
    c(res$design_info$best)
}

p = length(variable)
for (method in method.set) {
    for (analysis in analysis.set) {
        
        prefix = paste0(analysis, "_", method)
        fin = readRDS(paste0("results/", prefix, ".RDS"))
        
        #### Table ###########################################################
        est.set = matrix(unlist(lapply(fin, est_extract)), ncol=p, byrow=TRUE)
        se.set = matrix(unlist(lapply(fin, se_extract)), ncol=p, byrow=TRUE)
        
        bias = colMeans(est.set)-true.value
        se = apply(est.set, 2, sd)
        see = colMeans(se.set)
        cp = rep(NA, p)
        for (i in 1:p) {
            cp[i] = mean(est.set[,i]-q975*se.set[,i] <= true.value[i] & est.set[,i]+q975*se.set[,i] >= true.value[i])
        }
        
        out = data.frame(Variable=variable, Bias=bias, SE=se, SEE=see, CP=cp)
        
        write.table(out, file=paste0("sum2_", prefix, ".tab"), row.names=FALSE, quote=FALSE, sep="\t")
        #### Table ###########################################################
        
#         #### Triplot #########################################################
#         design.set = matrix(unlist(lapply(fin, design_extract)), ncol=3, byrow=TRUE)
#         total = rowSums(design.set)
#         p1 = design.set[,1]/total
#         p2 = design.set[,2]/total
#         p3 = design.set[,3]/total
#         triplot(p1, p2, p3, label=c("V=1", "V=2", "V=3"), main)
#         #### Triplot #########################################################
    }
}

