library(klaR)

method.set = c("ACML")
NSIM = 200
n.digits = 3
table.sep = "\t"

est_extract = function (res) {
    c(res$ests$coef$beta, NA, res$ests$coef$alpha)
}

se_extract = function (res) {
    c(sqrt(diag(res$ests$vcov$beta)), NA, sqrt(diag(res$ests$vcov$alpha)))
}

design_extract = function (res) {
    c(res$design_info$best)
}

est_extract_fixed = function (res) {
    c(res$single$coef$beta, NA, res$single$coef$alpha)
}

se_extract_fixed = function (res) {
    c(sqrt(diag(res$single$vcov$beta)), NA, sqrt(diag(res$single$vcov$alpha)))
}

est_extract_fc = function (res) {
    c(res$fc$coef$beta, NA, res$fc$coef$alpha)
}

se_extract_fc = function (res) {
    c(sqrt(diag(res$fc$vcov$beta)), NA, sqrt(diag(res$fc$vcov$alpha)))
}

#### Design Table ###########################################################
covariate = c("Intercept", "Visit", "SNP", "SNP $\\times$ Visit",
              "Cigarettes/day (per 10)", "Packs/years (per 20)",
              "Female", "Age (per 10 years)", "Baseline BMI (per 5 $\\text{kg}/\\text{m}^2$)",
              "BMI change (per 5 $\\text{kg}/\\text{m}^2$)", "Dependence parameter",
              "{  }$\\gamma$", "{  }$\\log(\\sigma)$")

for (method in method.set) {
    fn.out = paste0("sum2_", method, ".tab")
    sink(fn.out)
    cat("", "", "Adaptive design", "", "", sep=table.sep); cat("\n")
    cat("", "Fixed design", "Efficiency $\\beta_e$", "Efficiency $(\\beta_e, \\beta_{et})$", "Full cohort", sep=table.sep); cat("\n")
    sink()
    
    prefix = paste0("adaptive_design_Xe_", method)
    fin = readRDS(paste0("results/", prefix, ".RDS"))
    
    est.fixed = colMeans(matrix(unlist(lapply(fin, est_extract_fixed)), nrow = NSIM, byrow = TRUE))
    se.fixed = colMeans(matrix(unlist(lapply(fin, se_extract_fixed)), nrow = NSIM, byrow = TRUE))
    fixed = paste0("$", round(est.fixed, n.digits), "$ $(", round(se.fixed, n.digits), ")$")
    fixed[11] = ""
    
    est.xe = colMeans(matrix(unlist(lapply(fin, est_extract)), nrow = NSIM, byrow = TRUE))
    se.xe = colMeans(matrix(unlist(lapply(fin, se_extract)), nrow = NSIM, byrow = TRUE))
    xe = paste0("$", round(est.xe, n.digits), "$ $(", round(se.xe, n.digits), ")$")
    xe[11] = ""

    prefix = paste0("adaptive_design_XeXei_", method)
    fin = readRDS(paste0("results/", prefix, ".RDS"))
    
    est.xexei = colMeans(matrix(unlist(lapply(fin, est_extract)), nrow = NSIM, byrow = TRUE))
    se.xexei = colMeans(matrix(unlist(lapply(fin, se_extract)), nrow = NSIM, byrow = TRUE))
    xexei = paste0("$", round(est.xexei, n.digits), "$ $(", round(se.xexei, n.digits), ")$")
    xexei[11] = ""
    
    est.fc = colMeans(matrix(unlist(lapply(fin, est_extract_fc)), nrow = NSIM, byrow = TRUE))
    se.fc = colMeans(matrix(unlist(lapply(fin, se_extract_fc)), nrow = NSIM, byrow = TRUE))
    fc = paste0("$", round(est.fc, n.digits), "$ $(", round(se.fc, n.digits), ")$")
    fc[11] = ""
    
    sink(fn.out, append=TRUE)
    for (i in seq(length(covariate))) {
        cat(covariate[i], fixed[i], xe[i], xexei[i], fc[i], sep=table.sep); cat("\n")
    }
    sink() 
}
#### Design Table ###########################################################

#### Triplot ################################################################
# png("sum2_triplot.png", units = "in", res = 500, height = 4.4, width = 8)
# layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), heights = c(4, 0.4))

png("sum2_triplot.png", units = "in", res = 500, height = 8.8, width = 8)
layout(matrix(c(1, 2), 2, 1, byrow = TRUE), heights = c(8, 0.8))

fin = readRDS("results/adaptive_design_Xe_ACML.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
total = rowSums(design.set)
p1 = design.set[,1]/total
p2 = design.set[,2]/total
p3 = design.set[,3]/total
triplot(p1, p2, p3, label = c("V = 1", "V = 2", "V = 3"), main = "ACML")

fin = readRDS("results/adaptive_design_XeXei_ACML.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
total = rowSums(design.set)
p1 = design.set[,1]/total
p2 = design.set[,2]/total
p3 = design.set[,3]/total
tripoints(p1, p2, p3, col = "grey")
p.fixed = c(40, 70, 40)/150
tripoints(p.fixed, pch = 13)

# fin = readRDS("results/adaptive_design_Xe_MI.RDS")
# design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
# total = rowSums(design.set)
# p1 = design.set[,1]/total
# p2 = design.set[,2]/total
# p3 = design.set[,3]/total
# triplot(p1, p2, p3, label = c("V = 1", "V = 2", "V = 3"), main = "MI")
# 
# fin = readRDS("results/adaptive_design_XeXei_MI.RDS")
# design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
# total = rowSums(design.set)
# p1 = design.set[,1]/total
# p2 = design.set[,2]/total
# p3 = design.set[,3]/total
# tripoints(p1, p2, p3, col = "grey")
# p.fixed = c(40, 70, 40)/150
# tripoints(p.fixed, pch = 13)
#
# par(mar = c(0,0,0,0))
# plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
# legend("top", col = c("black", "black", "grey"), pch = c(13, 1, 1), bty = "n",
# 	   legend=c("Fixed design",
# 				expression(paste("Adaptive design for ", italic(beta[e]))),
# 				expression(paste("Adaptive design for (", italic(beta[e]), ", ", italic(beta[et]), ")"))), horiz = TRUE)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top", col = c("black", "black", "grey"), pch = c(13, 1, 1), bty = "n",
	   legend=c("Fixed design",
				expression(paste("Adaptive design for ", italic(beta[e]))),
				expression(paste("Adaptive design for (", italic(beta[e]), ", ", italic(beta[et]), ")"))), horiz = TRUE)
dev.off()
#### Triplot ################################################################

# #### Sample Size Plot #######################################################
# png("sum2_adaptive_sample_size_Xei.png", units = "in", res = 500, height = 8, width = 8)
# par(mfrow = c(2, 1))
# 
# fin = readRDS("results/adaptive_sample_size_Xei_ACML.RDS")
# design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
# hist(design.set[,2], main = "ACML", xlab = expression(italic(n[2])), ylab = "", breaks = 50)
# m.n2 = mean(design.set[,2])
# abline(v = m.n2, lwd = 2, lty = 2)
# legend("topright", lty = 2, lwd = 2, bty = "n", legend=paste("mean = ", round(m.n2)))
# 
# fin = readRDS("results/adaptive_sample_size_Xei_MI.RDS")
# design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
# m.n2 = mean(design.set[,2])
# hist(design.set[,2], main = "MI", xlab = expression(italic(n[2])), ylab = "", breaks= 50)
# abline(v = m.n2, lwd = 2, lty = 2)
# legend("topright", lty = 2, lwd = 2, bty = "n", legend=paste("mean = ", round(m.n2)))
# 
# dev.off()
# #### Sample Size Plot #######################################################
