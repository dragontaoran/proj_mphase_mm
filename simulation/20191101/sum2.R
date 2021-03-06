library(klaR)

parameter = c("$\\beta_0$", "$\\beta_t$", "$\\beta_e$", "$\\beta_{et}$", "$\\beta_c$", "$\\gamma$")
true.value = c(-1.5, -0.25, 0.25, 0.25, 1, 2)
method.set = c("ACML", "MI")
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

est_extract_fixed = function (res) {
    c(res$single$coef$beta, res$single$coef$alpha)
}

se_extract_fixed = function (res) {
    c(sqrt(diag(res$single$vcov$beta)), sqrt(res$single$vcov$alpha))
}

#### Design Table ###########################################################
p = length(parameter)
fn.out.design = paste0("sum2_adaptive_design.tab")
sink(fn.out.design)
cat("", "", "Fixed design", rep("", 3), "", "Adaptive design for $\\beta_e$", rep("", 4), "", "Adaptive design for $(\\beta_{e}, \\beta_{et}$", rep("", 4), sep="\t"); cat("\n")
cat("Method", "Parameter", "Bias", "SE", "SEE", "CP", "", "Bias", "SE", "SEE", "CP", "RE", "", "Bias", "SE", "SEE", "CP", "RE", sep="\t"); cat("\n")
sink()
for (method in method.set) {
        
        prefix = paste0("adaptive_design_Xe_", method)
        fin = readRDS(paste0("results/", prefix, ".RDS"))
           
        est.set.fixed = matrix(unlist(lapply(fin, est_extract_fixed)), ncol=p, byrow=TRUE)
        se.set.fixed = matrix(unlist(lapply(fin, se_extract_fixed)), ncol=p, byrow=TRUE)
		
        bias.fixed = colMeans(est.set.fixed)-true.value
        se.fixed = apply(est.set.fixed, 2, sd)
        see.fixed = colMeans(se.set.fixed)
        cp.fixed = rep(NA, p)
        for (i in 1:p) {
            cp.fixed[i] = mean(est.set.fixed[,i]-q975*se.set.fixed[,i] <= true.value[i] & est.set.fixed[,i]+q975*se.set.fixed[,i] >= true.value[i])
        }
		
        est.set.xe = matrix(unlist(lapply(fin, est_extract)), ncol=p, byrow=TRUE)
        se.set.xe = matrix(unlist(lapply(fin, se_extract)), ncol=p, byrow=TRUE)
        
        bias.xe = colMeans(est.set.xe)-true.value
        se.xe = apply(est.set.xe, 2, sd)
        see.xe = colMeans(se.set.xe)
        cp.xe = rep(NA, p)
        for (i in 1:p) {
            cp.xe[i] = mean(est.set.xe[,i]-q975*se.set.xe[,i] <= true.value[i] & est.set.xe[,i]+q975*se.set.xe[,i] >= true.value[i])
        }
		
		re.xe = (se.fixed/se.xe)^2
		
		prefix = paste0("adaptive_design_XeXei_", method)
        fin = readRDS(paste0("results/", prefix, ".RDS"))
		
        est.set.xexei = matrix(unlist(lapply(fin, est_extract)), ncol=p, byrow=TRUE)
        se.set.xexei = matrix(unlist(lapply(fin, se_extract)), ncol=p, byrow=TRUE)
        
        bias.xexei = colMeans(est.set.xexei)-true.value
        se.xexei = apply(est.set.xexei, 2, sd)
        see.xexei = colMeans(se.set.xexei)
        cp.xexei = rep(NA, p)
        for (i in 1:p) {
            cp.xexei[i] = mean(est.set.xexei[,i]-q975*se.set.xexei[,i] <= true.value[i] & est.set.xexei[,i]+q975*se.set.xexei[,i] >= true.value[i])
        }
		
		re.xexei = (se.fixed/se.xexei)^2
		
		sink(fn.out.design, append=TRUE)
		for (i in 1:p) {
			cat(ifelse(i == 1, method, ""), parameter[i], bias.fixed[i], se.fixed[i], see.fixed[i], cp.fixed[i], 
				"", bias.xe[i], se.xe[i], see.xe[i], cp.xe[i], re.xe[i],
				"", bias.xexei[i], se.xexei[i], see.xexei[i], cp.xexei[i], re.xexei[i], sep="\t"); cat("\n")
		}
		sink() 
}
#### Design Table ###########################################################

#### Sample Size Table ######################################################
prefix = paste0("adaptive_sample_size_Xei_ACML")
fin = readRDS(paste0("results/", prefix, ".RDS"))

est.set.acml = matrix(unlist(lapply(fin, est_extract)), ncol=p, byrow=TRUE)
se.set.acml = matrix(unlist(lapply(fin, se_extract)), ncol=p, byrow=TRUE)

bias.acml = colMeans(est.set.acml)-true.value
se.acml = apply(est.set.acml, 2, sd)
see.acml = colMeans(se.set.acml)
cp.acml = rep(NA, p)
for (i in 1:p) {
	cp.acml[i] = mean(est.set.acml[,i]-q975*se.set.acml[,i] <= true.value[i] & est.set.acml[,i]+q975*se.set.acml[,i] >= true.value[i])
}

prefix = paste0("adaptive_sample_size_Xei_MI")
fin = readRDS(paste0("results/", prefix, ".RDS"))

est.set.mi = matrix(unlist(lapply(fin, est_extract)), ncol=p, byrow=TRUE)
se.set.mi = matrix(unlist(lapply(fin, se_extract)), ncol=p, byrow=TRUE)

bias.mi = colMeans(est.set.mi)-true.value
se.mi = apply(est.set.mi, 2, sd)
see.mi = colMeans(se.set.mi)
cp.mi = rep(NA, p)
for (i in 1:p) {
	cp.mi[i] = mean(est.set.mi[,i]-q975*se.set.mi[,i] <= true.value[i] & est.set.mi[,i]+q975*se.set.mi[,i] >= true.value[i])
}

fn.out.sample.size = paste0("sum2_adaptive_sample_size_Xei.tab")
sink(fn.out.sample.size)
cat("", "ACML", rep("", 4), "MI", rep("", 3), sep="\t"); cat("\n")
cat("Parameter", "Bias", "SE", "SEE", "CP", "", "Bias", "SE", "SEE", "CP", sep="\t"); cat("\n")
for (i in 1:p) {
	cat(parameter[i], bias.acml[i], se.acml[i], see.acml[i], cp.acml[i], "", bias.mi[i], se.mi[i], see.mi[i], cp.mi[i], sep="\t"); cat("\n")
}
sink()
#### Sample Size Table ######################################################

#### Triplot ################################################################
png("sum2_triplot.png", units = "in", res = 500, height = 4.4, width = 8)
layout(matrix(c(1, 2, 3, 3), 2, 2, byrow = TRUE), heights = c(4, 0.4))
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

fin = readRDS("results/adaptive_design_Xe_MI.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
total = rowSums(design.set)
p1 = design.set[,1]/total
p2 = design.set[,2]/total
p3 = design.set[,3]/total
triplot(p1, p2, p3, label = c("V = 1", "V = 2", "V = 3"), main = "MI")

fin = readRDS("results/adaptive_design_XeXei_MI.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
total = rowSums(design.set)
p1 = design.set[,1]/total
p2 = design.set[,2]/total
p3 = design.set[,3]/total
tripoints(p1, p2, p3, col = "grey")
p.fixed = c(40, 70, 40)/150
tripoints(p.fixed, pch = 13)

par(mar = c(0,0,0,0))
plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
legend("top", col = c("black", "black", "grey"), pch = c(13, 1, 1), bty = "n",
	   legend=c("Fixed design",
				expression(paste("Adaptive design for ", italic(beta[e]))),
				expression(paste("Adaptive design for (", italic(beta[e]), ", ", italic(beta[et]), ")"))), horiz = TRUE)
dev.off()
#### Triplot ################################################################

#### Sample Size Plot #######################################################
png("sum2_adaptive_sample_size_Xei.png", units = "in", res = 500, height = 8, width = 8)
par(mfrow = c(2, 1))

fin = readRDS("results/adaptive_sample_size_Xei_ACML.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
hist(design.set[,2], main = "ACML", xlab = expression(italic(n[2])), ylab = "", breaks = 50)
m.n2 = mean(design.set[,2])
abline(v = m.n2, lwd = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, bty = "n", legend=paste("mean = ", round(m.n2)))

fin = readRDS("results/adaptive_sample_size_Xei_MI.RDS")
design.set = matrix(unlist(lapply(fin, design_extract)), ncol = 3, byrow = TRUE)
m.n2 = mean(design.set[,2])
hist(design.set[,2], main = "MI", xlab = expression(italic(n[2])), ylab = "", breaks= 50)
abline(v = m.n2, lwd = 2, lty = 2)
legend("topright", lty = 2, lwd = 2, bty = "n", legend=paste("mean = ", round(m.n2)))

dev.off()
#### Sample Size Plot #######################################################
