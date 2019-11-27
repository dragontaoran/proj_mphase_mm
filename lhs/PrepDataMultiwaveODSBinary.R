remove(list=ls())
mean.not.na <- function(x) mean(x, na.rm=TRUE)
cluster.summary <- function( id, x, fun ){ xlist <- split( x, id )
                                           nj    <- unlist( lapply( xlist, length ) )
                                           xj <- unlist( lapply( xlist, fun) )
                                           xsummary <- rep( xj, nj )
                                           xsummary}
FirstValue <- function(Var, id){ xlist  <- split( Var, id )
                                 nj     <- unlist( lapply( xlist, length ) )
                                 dup    <-  duplicated(id)
                                 val    <- Var[!dup]
                                 out    <- rep(val, nj)
                                 out}
length.unique<- function(x) length(unique(x))


path.work <- "~/rsch/OutDepSamp/Projects/ODSforEpi/DataManagement/"
path.import <- "~/rsch/OutDepSamp/Projects/ODSforEpi/ODSContinuous/Analysis"
setwd(path.work)

## load phenotype data.  Object is called dat
load("AnalysisData.Rdata")
names(dat.pheno)[which(names(dat.pheno)=="geneva_id")] <- "GenevaID"
dat.pheno <- dat.pheno[order(dat.pheno$GenevaID, dat.pheno$visit),]
length(unique(dat.pheno$id))

## load genotype data created by the find_snps program.
load("blah.Rdata")

## Bring in rs177852 from Hansel 2013
## Only keep the 'wild-type' columns (wild-type or 0 risk alleles, one risk allele, two risk alleles)
COL <- ncol(dat$hansel2013)
NewHansel <- data.frame(t(dat$hansel2013[,c(seq(4, COL, 3))]))
names(NewHansel) <- dat$hansel2013[,1]
NewHansel$GenevaID <- as.numeric(as.character(gsub("X","", row.names(NewHansel))))
NewHansel <- NewHansel[order(NewHansel$GenevaID),]

## We want this to be a dominant model for the risk allele according tho table 5 of the Lee paper
NewHansel$rs177852 <- abs(1-NewHansel$rs177852)

# Merge LHS and Hansel 2013 SNP data
lhs.new <- merge(dat.pheno, NewHansel, by="GenevaID")
lhs.new <- lhs.new[order(lhs.new$id, lhs.new$visit),]
dim(lhs.new)
length(unique(lhs.new$GenevaID))

lhs.new <- lhs.new[order(lhs.new$id, lhs.new$visit),]

# remove those without the snp data
lhs.new <- lhs.new[!is.na(lhs.new$rs177852),]
length(unique(lhs.new$id))

#lhs.new <- lhs.new[!(lhs.new$ivquad %in% c(1,2,3)),]

dim(lhs.new)
length(unique(lhs.new$id))
table(table(lhs.new$id))
lhs.new <- lhs.new[ !is.na(lhs.new$fev) & !is.na(lhs.new$bmi),]
dim(lhs.new)
length(unique(lhs.new$id))
table(table(lhs.new$id))
lhs.new <- lhs.new[    !is.na(lhs.new$packyear),]
dim(lhs.new)
length(unique(lhs.new$id))
table(table(lhs.new$id))
lhs.new <- lhs.new[!(lhs.new$ivquad %in% c(1,2,3)),]
dim(lhs.new)
length(unique(lhs.new$id))
table(table(lhs.new$id))

lhs.new[,c("id","visit","ivquad","ivquad_f")]

# 
# lhs.new <- lhs.new[!(lhs.new$ivquad %in% c(1,2,3)) & 
#                      !is.na(lhs.new$packyear) & 
#                      !is.na(lhs.new$fev) & 
#                      !is.na(lhs.new$fev) & 
#                      !is.na(lhs.new$bmi),]
dim(lhs.new)
dat <- lhs.new
length(unique(dat$id))
#dat$n.obs <- cluster.summary(dat$id, dat$id, length)
#dat <- dat[dat$n.obs>3,]

# Since we are only looking at people who never quit smoking, we will use baseline data for the outcome.
dat <- dat[order(dat$id, dat$visit),]
dup <- duplicated(dat$id)
dat.1 <- dat[!dup,]
dat.1$visit <- 0
dat.1$fev <- dat.1$fev0
dat.1$fvc <- dat.1$fvc0
dat.1$fevpct <- dat.1$fevpct0
dat.1$bmi <- dat.1$bmi0
dat <- rbind(dat, dat.1)
dat <- dat[order(dat$id, dat$visit),]
#dat <- dat[,c("id","visit","fev","fevpct","fvc","f31cigs","packyear","age","bmi0","bmi","female","rs177852","rs4330400","rs12194741","rs7331728","rs8179183","rs7911302","rs10761570","site", "fevpct0")]
dim(dat[is.na(dat$fevpct),])
dat <- dat[!is.na(dat$fevpct),]

dat$fevpctdec <- dat$fevpct0-dat$fevpct
dat$fevpctdec3 <- ifelse(dat$fevpctdec>3,1,0)
dat$fevpctdec4 <- ifelse(dat$fevpctdec>4,1,0)
dat$fevpctdec5 <- ifelse(dat$fevpctdec>5,1,0)
dat$fevpctdec6 <- ifelse(dat$fevpctdec>6,1,0)
dat$fevpctdec7 <- ifelse(dat$fevpctdec>7,1,0)
dat$fevpctdec8 <- ifelse(dat$fevpctdec>8,1,0)
# Create variables used for modeling
dat$age.10      <- (dat$age-50)/10
dat$cigs0.10    <- (dat$f31cigs-30)/10
dat$pack0.20    <- (dat$packyear-40)/20
dat$bmi5        <- (dat$bmi-25)/5
dat$bmi0.5      <- (dat$bmi0-25)/5
dat$bmi5.change <- dat$bmi5-dat$bmi0.5
dat$site2       <- as.integer(dat$site==2); dat$site3 <- as.integer(dat$site==3);dat$site4 <- as.integer(dat$site==4);
dat$site5       <- as.integer(dat$site==5); dat$site6 <- as.integer(dat$site==6);dat$site7 <- as.integer(dat$site==7);
dat$site8       <- as.integer(dat$site==8); dat$site9 <- as.integer(dat$site==9);dat$site10 <- as.integer(dat$site==10);

library(MMLB)


dat <- dat[dat$visit>0,]
dat       <- dat[order(dat$id, dat$visit),]
dat$n.obs <- cluster.summary(dat$id, dat$id, length)
dat2 <- dat[dat$n.obs>1,]
length(unique(dat2$id))

hist(tapply(dat2$fevpctdec3, dat2$id, sum))
hist(tapply(dat2$fevpctdec4, dat2$id, sum))
hist(tapply(dat2$fevpctdec5, dat2$id, sum))
hist(tapply(dat2$fevpctdec6, dat2$id, sum))
hist(tapply(dat2$fevpctdec7, dat2$id, sum))

L <- length(grep("rs", names(dat), value=TRUE))
blah <- NULL
for (i in c(37:40)){
    dat2$snp <- dat2[[grep("rs", names(dat2), value=TRUE)[i]]]
    print(date())
    FC <- mm(fevpctdec5~snp*I(visit-1)+cigs0.10+pack0.20+female+age.10+bmi0.5+bmi5.change,#+site2+site3+site4+site5+site6+site7+site8+site9+site10, 
             lv.formula=~1, t.formula=~1, id=id, data=dat2)
    print(date())
    blah <- rbind(blah, c(grep("rs", names(dat), value=TRUE)[i], summary(FC)$mean.table[2,c(1,2)],summary(FC)$mean.table[10,c(1,2)]))
    print(blah)
    
    
}
# library(lme4)
# FC <- lmer(fevpct~(rs177852+female+age.10)*(t2+t3+t4+t5)+site2+site3+site4+site5+site6+site7+site8+site9+site10+bmi0.5+bmi5.change+cigs0.10+pack0.20+
#                (visit|id), REML=FALSE, data=dat2)
# FC1 <- lmer(fevpct~rs177852*factor(visit)+#rs177852*(I(visit==2)+I(visit>2))+
#                 (female+age.10+site2+site3+site4+site5+site6+site7+site8+site9+site10)*visit+bmi0.5+bmi5.change+cigs0.10+pack0.20+
#                (visit|id), REML=FALSE, data=dat2)


lhsMWODSBinary <- dat2
save(lhsMWODSBinary, file="lhsMWODSBinary.Rdata", compress=TRUE)

Estimate     Model SE   Estimate     Model SE  
[1,] "rs10443259" 0.04840107   0.08515746 -0.02041119  0.0258186 
[2,] "rs1045895"  -0.001277156 0.08810314 0.01940604   0.02677694
[3,] "rs10493377" 0.05970404   0.09276562 -0.01985144  0.02805038
[4,] "rs10889558" -0.1104373   0.1625535  -0.008962492 0.04956517
[5,] "rs10889562" 0.01024388   0.08758477 -0.0137338   0.02651095
[6,] "rs1137100"  0.0584116    0.08559477 -0.01671199  0.02595209
[7,] "rs1171265"  -0.1695996   0.1306607  0.04492673   0.0396231 
[8,] "rs1171271"  -0.1514283   0.1642485  0.00769283   0.04971601
[9,] "rs1171274"  -0.1131146   0.1626329  -0.008648636 0.04957423
[10,] "rs1171279"  -0.01987334  0.0855731  0.007019839  0.0259306 
[11,] "rs12564626" -0.05195038  0.09280629 0.0199583    0.02812242
[12,] "rs1327118"  0.05698994   0.09431817 -0.02341469  0.02848537
[13,] "rs1327121"  -0.1402657   0.1399166  0.008087922  0.04298884
[14,] "rs17097182" 0.2887818    0.1510802  -0.06769698  0.0457815 
[15,] "rs17127652" -0.1286161   0.2730721  0.03155018   0.08350536
[16,] "rs17412682" -0.004284123 0.09478242 0.01397187   0.0288583 
[17,] "rs1782754"  -0.1540975   0.1639196  0.007492286  0.04957795
[18,] "rs1782763"  -0.144896    0.1366931  0.04506898   0.0415847 
[19,] "rs1805096"  -0.01116954  0.08836126 0.0472151    0.02691161
[20,] "rs1892535"  0.01223971   0.08947775 0.01711927   0.02708417
[21,] "rs1938484"  0.0284351    0.09028473 0.004059158  0.02732163
[22,] "rs2025804"  -0.1478082   0.140504   0.01036191   0.04312518
[23,] "rs2275456"  0.1572411    0.1589563  -0.01383564  0.04774741
[24,] "rs3828034"  -0.05607199  0.08917182 0.02863965   0.02682585
[25,] "rs4655680"  0.05252538   0.0852507  -0.02176925  0.0258477 
[26,] "rs4655811"  -0.1162568   0.1401877  0.01169747   0.04311403
[27,] "rs6588153"  -0.009531174 0.08862171 0.04743313   0.02698959
[28,] "rs6657868"  -0.1141726   0.1218431  0.03873759   0.03738041
[29,] "rs6691346"  0.0456104    0.08516669 -0.01940959  0.02582236
[30,] "rs6694528"  -0.137654    0.09817587 0.03883856   0.02977446
[31,] "rs7531867"  -0.0112904   0.08835835 0.04727188   0.0269107 
[32,] "rs8179183"  -0.07258007  0.0900444  0.04380355   0.02704381
[33,] "rs9436299"  -0.07238039  0.1536809  0.01226351   0.04734576
[34,] "rs9436746"  -0.09910354  0.1247899  0.03635028   0.0380371 
[35,] "rs970468"   -0.1593982   0.1242726  0.0460933    0.03813818
[36,] "rs4330400"  0.07225487   0.08721097 -0.0480955   0.02640657

[37,] "rs12194741" -0.09112188  0.0871274  0.07083717   0.02626083
[38,] "rs10761570" -0.1764366   0.1114001  0.08897165   0.03427893
[39,] "rs4948444"  -0.1448572   0.1091613  0.06775831   0.03344413
[40,] "rs7911302"  -0.1699972   0.1112409  0.09078587   0.03426601

[41,] "rs7922793"  -0.07279685  0.1300411  0.07629815   0.0404093 
[42,] "rs7331728"  0.06110416   0.09202705 -0.05797351  0.02801407
[43,] "rs177852"   -0.05138858  0.08562164 0.05319245   0.02589687
