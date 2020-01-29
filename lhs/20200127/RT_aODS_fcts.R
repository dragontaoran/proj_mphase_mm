# aODS functions
# Nate Mercaldo, Jonathan Schildcrout, Ran Tao
# Code modified by Ran Tao is labeled with "#### This part is rewritten by RT ####"
# load libraries
library(MMLB)
library(MASS)
library(mitools)

gen_XI = function(SEED, NLOW, NHIGH, N, PXE, OMEGAZ.Xe) {
  # Arguments
  #
  # SEED  = seed for rng
  # NLOW  = minumum cluster size
  # NHIGH = maximum cluster size
  # N     = number of subjects to sample
  # PXE   = proprotion of Xe in population
  # OMEGAZ.Xe = Mean difference in Z for Xe=1 vs Xe=0
  
  # Output: design matrix
  
  # Set seed
  set.seed(SEED)
  
  # Generate number of clusters per subjects
  nclust = sample(seq(NLOW, NHIGH), N, replace=TRUE)
  if (NLOW == NHIGH) {
    nclust = rep(NLOW, N)
  }
  
  # Generate key exposure (Xe) and a possibly related binary variable Z 
  Xe = rbinom(N, size=1, prob=PXE)
  Z_Xe0 = rnorm(sum(Xe==0), 0, 1)
  Z_Xe1 = rnorm(sum(Xe==1), OMEGAZ.Xe, 1)
  Z = rep(NA, N)
  Z[Xe==0] = Z_Xe0
  Z[Xe==1] = Z_Xe1
  
  XI = vector('list', N)
  for (i in seq(N)) {
	#### This part is changed by RT ####
	time.tmp = seq(nclust[i])-1
    XI[[i]] = cbind(i, rep(1, nclust[i]), time.tmp, Xe[i], time.tmp*Xe[i], Z[i])
	#### This part is changed by RT ####
  }
  XI = data.frame(do.call(rbind, XI))
  colnames(XI) = c('id', 'int', 'time', 'Xe', 'Xe.i', 'Z')
  XI
}

gen_data = function(SEED, NLOW, NHIGH, N, PXE, MEANFORMULA, LVFORMULA, TFORMULA, BETAZ, BETAM, GAMMA, SIGMA, Q, XI=NULL) {
  # Set seed
  set.seed(SEED)
  
  if (is.null(XI)) {
    stop('Specify XI')
  }
  nclust = as.numeric(table(XI$id))
  XI$Y = 0 # Just add junk column for Y, this will be rewritten later!
  
  if (is.character(MEANFORMULA)) {
    MEANFORMULA = as.formula(MEANFORMULA)
  }
  if (is.character(LVFORMULA)) {
    LVFORMULA = as.formula(LVFORMULA)
  }
  if (is.character(TFORMULA)) {
    TFORMULA = as.formula(TFORMULA)
  }
  
  # Generate Y
  ests = BETAM
  if (!is.null(GAMMA)) {
    ests = c(ests, GAMMA)
  }
  if (!is.null(SIGMA)) {
    ests = c(ests, log(SIGMA))
  }
  
  if (is.data.frame(BETAM) | is.data.frame(BETAZ) ) {
    beta.mz = c(as.numeric(BETAM), as.numeric(BETAZ))
  } else {
    beta.mz = c(BETAM, BETAZ)
  }
  
  sdat = GenBinaryY(mean.formula=MEANFORMULA, lv.formula=LVFORMULA, t.formula=TFORMULA, 
                    beta=beta.mz, sigma=SIGMA, gamma=GAMMA, id=id, data=XI, q=Q, Yname = "Y") 
  sdat_s = lapply(split(sdat, sdat$id), function(ZZ) {
    nY = sum(ZZ$Y)
    nn = nrow(ZZ)
    ZZ$ss = (nY==0)*1+(nY>0 & nY<nn)*2+(nY==nn)*3
    ZZ
  })
  sdat = do.call(rbind, sdat_s)
  attr(sdat, 'true') = ests
  sdat
}

idsXss = function(DAT) {
  # DAT = list of subject-level data (output from gen_data function)
  
  # Output: out (list of ids by each sampling stratum) and attr(out,'summ') which contains key frequencies
  #         idXss will be used when performing sampling to keep track of remaining ids
  
  # Generate ids by each sampling stratum
  N = length(DAT)
  DAT_ss = split(sapply(DAT, function(ZZ) ZZ$id[1]), sapply(DAT, function(ZZ) ZZ$ss[1]))
  
  # Calculate frequencies in each sampling statrum and frequencies of sum(Y) and mean(Y)
  ss_sum = sapply(DAT_ss, length)
  sy_sum = table(sapply(DAT, function(ZZ) sum(ZZ$Y)))
  sy_mean = table(sapply(DAT, function(ZZ) mean(ZZ$Y)))
  
  out = DAT_ss
  attr(out, 'freq') <- list('ss'=ss_sum, 'sy'=sy_sum, 'sm'=sy_mean)
  out
}

get_sample = function(DESIGN, IDSXSS, RSEXACT, REPLACE=FALSE) {
  # DESIGN  = vector of expected sample counts
  # IDSXSS  = list of ids by sampling stratum (length=3)
  # RSEXACT = indicator to force exact random sampling
  # REPLACE = indicator if sampling with replacement should be performed
  
  # Output: list that contains sampling probs, and lists of sampled ids by ss and remaining ids by ss
  
  tmp_design = as.list(DESIGN)
  blah_samp = vector('list', length(tmp_design))
  
  if (RSEXACT) {
    blah = do.call(rbind, mapply(function(AA, BB) list(cbind(AA, BB)), AA=as.list(as.numeric(names(IDSXSS))), BB=IDSXSS))
    blah_tmp = data.frame(blah[sort(sample(x=seq(nrow(blah)), round(sum(DESIGN), 0), replace=FALSE)),])
    if (ncol(blah_tmp) == 1) blah_tmp = data.frame(t(blah_tmp))
    blah_samp = blah_tmp
    blah_samp[,1] = factor(blah_samp[,1], levels=1:3)
    blah_samp = split(blah_samp$BB, blah_samp$AA)
  }
  
  tmp_samp = mapply(function(AA, BB, CC) {
    nn = length(AA)
    if (nn < BB) BB = nn
    sp = ifelse(nn==0, 0, BB/nn)
    sampled = which(rbinom(nn, size=1, prob=sp) == 1)
    if (length(sampled) > 0 | !is.null(CC)) {
      if (is.null(CC)) {
        sampled = AA[sampled]
      } else {
        sampled = CC
      }
    }
    attr(sampled,'sp') = sp
    list(sampled)}, AA=IDSXSS, BB=tmp_design, CC=blah_samp)
  
  tmp_sp = sapply(tmp_samp, function(ZZ) attr(ZZ, 'sp'))
  
  tmp_remain = mapply(function(AA, BB) {
    out = AA
    if (length(BB) > 0)  out = AA[!(AA%in%BB)]
    out
  }, AA=IDSXSS, BB=tmp_samp)
  
  list('sp'=tmp_sp, 'idsXss_samp'=tmp_samp, 'idsXss_remain'=tmp_remain)
}

prep_data = function(mean.formula, lv.formula, t.formula, id, data, inits, samp.probs, samp.probi, offset) {
  
  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  
  y  = model.response(mean.f, 'numeric') 
  uy = unique(y)
  
  x  = model.matrix(mean.formula, mean.f)
  
  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if (!is.null(t.formula)) {
    x.t  = model.matrix(t.formula, model.frame(t.formula, data))
  }
  if (!is.null(lv.formula)) {
    x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))
  }
  
  if (!is.matrix(samp.probs)) {
    samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)
  }
  if (is.null(samp.probi)) {
    samp.probi = matrix(1, nrow=length(y), ncol=1)
  }
  if (!is.null(inits)) {
    inits = unlist(inits)
  }
  
  if (is.null(inits)) {
    inits = c(glm(mean.formula, family='binomial',data=data)$coef, rep(1, ncol(x.t)+ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x = x[,-c(omit_dup_col)]
      inits = inits[-c(omit_dup_col)]
    }
  }
  
  if (is.null(offset)) {
    offset = rep(0, length(y))
  }
  
  x.t = cbind(x.t, NULL)
  x.lv = cbind(x.lv, NULL)
  
  nms.inits = names(inits)
  ix.mean = seq(length(nms.inits))
  ix.t = grep('gamma', nms.inits)
  ix.lv = grep('sigma', nms.inits)
  ix.mean = ix.mean[-c(which(ix.mean %in% c(ix.t, ix.lv)))]
  
  if (length(ix.lv) == 0) {
    inits.lv = -200 #### JS changed this 
  } else {
    inits.lv = inits[ix.lv]
  }
  
  if (length(ix.t) == 0) {
    inits.t = 0 
  } else {
    inits.t = inits[ix.t]
  }
  inits = c(inits[ix.mean], inits.t ,inits.lv)
  
  if (length(inits) != ncol(x) + ncol(x.t) + ncol(x.lv)) { 
    stop("Parameter length incongruous with X, Xgam, and Xsig")
  }
  paramlengths = c(ncol(x), ncol(x.t), ncol(x.lv))
  
  ## Lagged response
  Ylag = rep(0, length(y))
  for(i in 2:length(y)) {
    if (id[i] == id[i-1]) {
      Ylag[i] = y[i-1] 
    }
  }
  
  id.tmp = split(id, id)
  X.tmp = split(x, id)
  Y.tmp = split(y, id)
  Ylag.tmp = split(Ylag, id)
  Xgam.tmp = split(x.t, id)
  Xsig.tmp = split(x.lv, id)
  SampProbi.tmp = split(samp.probi, id)
  SampProbs.tmp = split(samp.probs, id)
  offset.tmp = split(offset, id)
  
  subjectData = list()
  uid = as.character(unique(id))
  for (j in seq(along=uid)) {
    i = uid[j]
    subjectData[[j]] = list(id=as.character(unique(id.tmp[[i]])), 
                            X=matrix(X.tmp[[i]], ncol=ncol(x)), 
                            Y=as.double(Y.tmp[[i]]),
                            Ylag=as.double(Ylag.tmp[[i]]),
                            Xgam=matrix(Xgam.tmp[[i]], ncol=ncol(x.t)),
                            Xsig=matrix(Xsig.tmp[[i]], ncol=ncol(x.lv)),
                            SampProbi=unique(SampProbi.tmp[[i]]),
                            SampProbs=matrix(SampProbs.tmp[[i]], ncol=ncol(samp.probs))[1,],
                            Offset=as.double(offset.tmp[[i]]))
  }
  names(subjectData) = uid
  
  beta_nms = colnames(x)
  alpha_nms = c(if (!is.null(t.formula)) { paste('gamma', colnames(x.t), sep=':') }, 
                if (!is.null(lv.formula)) {paste('log(sigma)', colnames(x.lv), sep=':')}) 
  
  list(params=inits, paramlengths=paramlengths, subjectData=subjectData, samp.probi=tapply(samp.probi, id, unique),
       nms=list(beta_nms, alpha_nms))
}

Odds2Prob = function(odds) {
  odds/(1+odds)
}

ImputeData = function(FIT, DAT, M=30, ROB=FALSE, MARG.EXP.FORMULA, verbose=FALSE, EXPVAR=EXPVAR, Sampled=Sampled) {
  # FIT - MMLongit object
  # DAT - complete data set with sampled flag (Sampled)
  # M - # of imputations
  # ROB - indicator, if FALSE use model-based covariance
  # MARG.EXP.FORMULA - marginal exposure model formula: [Xe | Xo]
  # EXPVAR - exposure / expensive variable: Xe
  # Sampled - should be set to 1 if included in the sample and 0 otherwise
  
  Sampled  = DAT$Sampled = DAT[,as.character(Sampled)]
  imp_data = vector('list', M)
  
  for (ixx in seq(M)) { 
    if (verbose) cat("\r","Imputing dataset number:", ixx)
    # Step 0: Update date files
    # Create versions of DAT with updated Xe values
    DAT_Xe1 = DAT_Xe0 = DAT
    
    DAT_Xe1[,EXPVAR] = 1
    DAT_Xe1[,"Xe.i"] = DAT_Xe1[,"time"]
    
    DAT_Xe0[,EXPVAR] = 0
    DAT_Xe0[,"Xe.i"] = 0

    # Non-sampled subjects
    DAT_S0 = DAT[Sampled==0,]
    DAT_S0_Xe1 = DAT_Xe1[Sampled==0, ]
    DAT_S0_Xe0 = DAT_Xe0[Sampled==0, ]
    dup = duplicated(DAT_S0$id)
    DAT_S0_v0 = DAT_S0[!dup, ]
    
    # Sampled subjects
    DAT_S1     <- DAT[Sampled==1,]
    DAT_S1_Xe1 <- DAT_Xe1[Sampled==1, ]
    DAT_S1_Xe0 <- DAT_Xe0[Sampled==1, ]
    dup        <- duplicated(DAT_S1$id)
    DAT_S1_v0  <- DAT_S1[!dup,]
    
    # Step 1: extract estimates from FIT and draw m^th theta from MVN
    
    theta = c(FIT$beta, FIT$alpha)
    if (ROB) { 
      cov_theta = FIT$rob.cov
    } else { 
      cov_theta = FIT$mod.cov
    }
    
    theta.m = mvrnorm(n=1, mu=theta, Sigma=cov_theta)
    
    nms.inits = names(theta.m)
    ix.mean = seq(length(nms.inits))
    ix.t = grep('gamma', nms.inits)
    ix.lv = grep('sigma', nms.inits)
    ix.mean = ix.mean[-c(which(ix.mean %in% c(ix.t, ix.lv)))]
    
    if (length(ix.lv) == 0) {
      inits.lv = -200
    } else {
      inits.lv = theta.m[ix.lv]
    }
    
    if (length(ix.t) == 0) {
      inits.t = 0
    } else {
      inits.t = theta.m[ix.t]
    }
    theta.m0 = c(theta.m[ix.mean], inits.t ,inits.lv)
    
    # Step 2: Using the estimated theta=(beta, alpha), calculate likelihood contributions for sampled subjects
    # for both Xe=1 and Xe=0.
    prep_input = attr(FIT, 'args')
    
    tmp_S1_Xe1 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, t.formula=prep_input$t.formula, 
                           id=DAT_S1_Xe1[,prep_input$id], data=DAT_S1_Xe1, inits=theta.m, samp.probs=prep_input$samp.probs, 
                           samp.probi=prep_input$samp.probi, offset=prep_input$offset)$subjectData
    tmp_S1_Xe0 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, t.formula=prep_input$t.formula, 
                           id=DAT_S1_Xe0[,prep_input$id], data=DAT_S1_Xe0, inits=theta.m, samp.probs=prep_input$samp.probs, 
                           samp.probi=prep_input$samp.probi, offset=prep_input$offset)$subjectData
    
    # These functions calculate the conditional likelihood contributions by and the ascertainment correction
    # for sampled subjects.  The difference of ascertainment corrections (i.e. log of the ratio 
    # log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo))= log(pr(S=1 | Xe=1, Xo))-log(pr(S=1 | Xe=0, Xo))) provides
    # the offset for the marginal exposure model 
    
    LLSC_1 = MMLB:::LogLScoreCalc(params=theta.m0, subjectData=tmp_S1_Xe1, Q=FIT$LLSC_args$Q, W=FIT$LLSC_args$W, Z=FIT$LLSC_args$Z, 
                                  ParamLengths=FIT$LLSC_args$ParamLengths, CondLike=FIT$control['cond.like']==1, EmpiricalCheeseCalc=ROB) 
    LLSC_0 = MMLB:::LogLScoreCalc(params=theta.m0, subjectData=tmp_S1_Xe0, Q=FIT$LLSC_args$Q, W=FIT$LLSC_args$W, Z=FIT$LLSC_args$Z, 
                                  ParamLengths=FIT$LLSC_args$ParamLengths, CondLike=FIT$control['cond.like']==1, EmpiricalCheeseCalc=ROB)
    
    # Note: The ascertainment corrections used for the marginal exposure model offset 
    # are already on the log scale, so log(a/b) = loga-logb
    logAC_1 = attr(LLSC_1, "ACSubj") 
    logAC_0 = attr(LLSC_0, "ACSubj")
    DAT_S1_v0$offset = logAC_1-logAC_0
    
    # Perform offsetted logistic regression [Note, only using first visit info i.e. no time-varying information.
    fit.exp = glm(MARG.EXP.FORMULA, offset=offset, data=DAT_S1_v0, family=binomial())
    
    # Extract alpha estimates and sample from the sampling distribution
    alpha = coef(fit.exp)
    cov_alpha = vcov(fit.exp)
    
    alpha.m = mvrnorm(n=1, mu=alpha, Sigma=cov_alpha)
    
    # Step 5: Apply the model and theta to calculate the ascertainment corrections for unsampled people, 
    
    nr = nrow(DAT_S0_Xe1)
    
    tmp_S0_Xe1 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, t.formula=prep_input$t.formula, 
                           id=DAT_S0_Xe1[,prep_input$id], data=DAT_S0_Xe1, inits=theta.m, samp.probs=prep_input$samp.probs[rep(1,nr),], 
                           samp.probi=prep_input$samp.probi[rep(1,nr)], offset=prep_input$offset[rep(1,nr)])$subjectData
    tmp_S0_Xe0 = prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, t.formula=prep_input$t.formula, 
                           id=DAT_S0_Xe0[,prep_input$id], data=DAT_S0_Xe0, inits=theta.m, samp.probs=prep_input$samp.probs[rep(1,nr),], 
                           samp.probi=prep_input$samp.probi[rep(1,nr)], offset=prep_input$offset[rep(1,nr)])$subjectData
    
    LLSC_1 = MMLB:::LogLScoreCalc(params=theta.m0, subjectData=tmp_S0_Xe1, Q=FIT$LLSC_args$Q, W=FIT$LLSC_args$W, Z=FIT$LLSC_args$Z, 
                                  ParamLengths=FIT$LLSC_args$ParamLengths, CondLike=FIT$control['cond.like']==1, EmpiricalCheeseCalc=ROB) 
    LLSC_0 = MMLB:::LogLScoreCalc(params=theta.m0, subjectData=tmp_S0_Xe0, Q=FIT$LLSC_args$Q, W=FIT$LLSC_args$W, Z=FIT$LLSC_args$Z, 
                                  ParamLengths=FIT$LLSC_args$ParamLengths, CondLike=FIT$control['cond.like']==1, EmpiricalCheeseCalc=ROB)
    
    ## Applied to the non-sampled subjects: log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo)).
    logAC_1 = attr(LLSC_1, "ACSubj") # Ascertainment corrections already log transformed
    logAC_0 = attr(LLSC_0, "ACSubj")
    offset_log_S0 = logAC_1-logAC_0
    
    # Create a temporary dataframe for non-sampled subjects in order to create the
    # model.matrix that is used to calculate marginal exposure odds
    tmp = DAT_S0_v0
    tmp[,EXPVAR] = sample(c(0,1), size=nrow(tmp), replace=TRUE)
    mm_tmp = model.matrix(as.formula(MARG.EXP.FORMULA), tmp) # Added as.formula
    
    # Compute the conditional exposure odds in non-sampled subjects
    # Notice that all conditional odds are equal: [Xe | Y, Xo, S=1] = [Xe | Y, Xo, S=0] = [Xe | Y, Xo] since we sampled based on Y
    # So, we multiply the odds (conditional on being sampled), pr[Xe=1 | Xo, S=1]/pr[Xe=0 | Xo, S=1] 
    # by the likelihood ratio (conditional on being sampled), pr[Y | Xe=1, Xo, S=1]/ pr[Y | Xe=1, Xo, S=1]
    # in the unsampled subjects in order to obtain conditional odds in unsampled: pr[Xe=1 | Y, Xo, S=1]/pr[Xe=0 | Y, Xo, S=1]
    marg.odds.Xe.Xo.S1 = exp(mm_tmp %*% alpha.m)*exp(offset_log_S0)  
    LR.Y.Xe.Xo.S1 = exp(attr(LLSC_1, "LogLikeSubj")-attr(LLSC_0, "LogLikeSubj"))
    cond.odds.Xe.Y.Xo.S0 = LR.Y.Xe.Xo.S1*marg.odds.Xe.Xo.S1  
    
    # Convert odds to probability, log(odds) = log(1/[1-p]) ---> p = odds/(1+odds)
    # Impute Xe for non-sampled subjects
    prXe1 = Odds2Prob(cond.odds.Xe.Y.Xo.S0)
    XeS0 = rbinom(nrow(tmp), size=1, prob=prXe1)
    
    # Add XeS0 info to DAT_S0_v0
    DAT_S0_v0[,EXPVAR] = XeS0
    DAT_S0 = DAT_S0[,-which(colnames(DAT_S0) == EXPVAR)]
    DAT_S0 = merge(DAT_S0, DAT_S0_v0[,c('id', EXPVAR)], by='id')
    DAT_S0$Xe.i = DAT_S0$Xe*DAT_S0$time
    
    DAT_imp = rbind(DAT_S1, DAT_S0[,colnames(DAT_S1)])
    
    # Store imputation data
    imp_data[[ixx]] = DAT_imp
  } 
  imp_data
}

get_design = function(N, M, INC=5) {
  # N   = vector of population sample sizes (A,B,C)
  # M   = target sample size
  # INC = increment size
  #### This part is rewritten by RT ####
  N1_seq = seq(0, min(N[1]-N[1]%%INC, M), INC)
  N2_seq = seq(0, min(N[2]-N[2]%%INC, M), INC)
  
  grids = expand.grid(N1_seq, N2_seq)
  grids = grids[rowSums(grids) <= M,]
  grids = cbind(grids, M-rowSums(grids))
  colnames(grids)[3] = "var3"
  grids = grids[grids[,3] <= N[3],]
  #### This part is rewritten by RT ####
  rownames(grids) = seq(nrow(grids))
  grids
}

get_opt = function(DAT) {
  if (class(DAT) == 'MMLongit') {
    var_gt = vcov(DAT)$beta['Xe.i','Xe.i']  
    var_g = vcov(DAT)$beta['Xe','Xe']
    #### This part is rewritten by RT ####
    d_opt = det(vcov(DAT)$beta[c('Xe', 'Xe.i'),c('Xe', 'Xe.i')])
    #### This part is rewritten by RT ####
    
  } else {
    # class = MIresult
    var_gt = vcov(DAT)['Xe.i','Xe.i']  
    var_g = vcov(DAT)['Xe','Xe']
    #### This part is rewritten by RT ####
    d_opt = det(vcov(DAT)[c('Xe', 'Xe.i'),c('Xe', 'Xe.i')])
    #### This part is rewritten by RT ####
  }
  c(var_gt, var_g, d_opt)
}

idea = function(iFIT, iDAT, iS1INFO, REMAINING, MARG.EXP.FORMULA, N2FIXED, METHOD, WAVED2, INC, 
                MODEL, LVMODEL, TMODEL, ITERLIM, MIREP, RNGp, VERBOSE, TARGETVAR, TARGET, BREP, Q, GRIDINC) {
    # Perform IDEA, use FA data to estimate Xe using data from previous stages
    if (VERBOSE) cat('\n\r',date(),'\n\r','Performing interim design evaluation analysis\n')
    
    #### This part is rewritten by RT ####
    # if (METHOD == 'ACML') MIREP = 1 # Overwrite MIREP if performing ACML analysis
    #### This part is rewritten by RT ####
    
    # Impute data (MIREP reps)
    #### This part is rewritten by RT ####
    ## change MIREP to BREP
    tmpXe = ImputeData(FIT=iFIT, DAT=iDAT, M=BREP, MARG.EXP.FORMULA=MARG.EXP.FORMULA, 
                     Sampled='samplek', verbose=VERBOSE, EXPVAR='Xe')
    #### This part is rewritten by RT ####
    
    # For each imputed data set, find best ss or design
    best_des_mat = matrix(NA, ncol=3, nrow=BREP)
    
    for(mxx in seq(BREP)) {
    
        if (N2FIXED) {
            # Adaptive design
            if(VERBOSE) cat('\n\r',date(),'\n\r','For fixed N2, calculating possible designs\n')
            
            #### This part is rewritten by RT #### 
            # Define candidate designs
            poss_des2 = get_design(N=REMAINING, M=WAVED2, INC=INC)
            # poss_des2 = poss_des2[poss_des2[,2]>0,]
            poss_des2 = data.frame(poss_des2)
            colnames(poss_des2) = c('n20','n21','n22')
            #### This part is rewritten by RT #### 
        } else {
            poss_vals_rng = ceiling(quantile(seq(1, REMAINING[2]), probs=c(RNGp[1], RNGp[2])))
            poss_vals = ceiling(seq(poss_vals_rng[1], poss_vals_rng[2], length.out=GRIDINC))
            poss_des2 = data.frame(n20 = 0, n21 = poss_vals, n22 = 0)
        }
    
        # Evaluate all possible candidate designs and store optimality summary
        optimal_mat = matrix(NA, nrow=nrow(poss_des2), ncol=3)
          
        for(opt_des in seq(nrow(poss_des2))) {
            adat2_tmp2 = tmpXe[[mxx]]
            s2d = poss_des2[opt_des,]
            
            #### This part is rewritten by RT ####
            adat2_tmp2_s0 = adat2_tmp2[which(adat2_tmp2$samplek == 0),]
            sdat2_tmp2_s0 = split(adat2_tmp2_s0, adat2_tmp2_s0$id)
            tmp2_s0_ss = idsXss(DAT=sdat2_tmp2_s0)
            s2info = get_sample(DESIGN=s2d, IDSXSS=tmp2_s0_ss, RSEXACT=FALSE)
            #### This part is rewritten by RT ####
            
            dat_s = split(adat2_tmp2, adat2_tmp2$id)
            dat_ix = which(adat2_tmp2$id %in% unlist(s2info$idsXss_samp))
            adat2_tmp2$samplek[dat_ix] = 2 
            adat2_tmp2[dat_ix,c('sprob1', 'sprob2', 'sprob3')] = matrix(s2info$sp, nrow=1)[rep(1, length(dat_ix)),]
            
            adat2_tmp3 = adat2_tmp2
            adat2_tmp2 = adat2_tmp2[adat2_tmp2$samplek>0,]
            
            # Estimate Stage K sampling probabilities
            prob_sampled = do.call(rbind, lapply(split(adat2_tmp2, adat2_tmp2$samplek), function(ZZ) ZZ[1,c('sprob1', 'sprob2', 'sprob3')]))
            prob_not_sampled = apply(1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
            prob_not_sampled = rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
            prob_sampled = prob_sampled*prob_not_sampled 
            prob_sampled$samplek = as.numeric(rownames(prob_sampled))
            colnames(prob_sampled)[1:3] = c('sp1', 'sp2', 'sp3') 
            
            adat2_tmp2 = merge(adat2_tmp2, prob_sampled, by='samplek', all.x=TRUE)
            adat2_tmp2 = adat2_tmp2[order(adat2_tmp2$id, adat2_tmp2$time),]
            
            fit2 = mm(as.formula(MODEL), dat=adat2_tmp2, id=id, lv.formula=LVMODEL, t.formula=TMODEL, cond.like=TRUE, 
                      samp.probs=as.matrix(adat2_tmp2[,c('sp1','sp2','sp3')]), iter.lim=ITERLIM, q=Q, return_args=TRUE)
            
            if (METHOD == 'MI') {
                # Perform MI analysis
                adat2_tmp3 = adat2_tmp3[order(adat2_tmp3$id, adat2_tmp3$time),]
                adat2_tmp3$Xe[adat2_tmp3$samplek == 0] = NA      # ignore Xe, Xe.i for non-sampled subjects
                adat2_tmp3$Xe.i[adat2_tmp3$samplek == 0] = NA
                adat2_tmp3$samplek[adat2_tmp3$samplek != 0] = 1    # ImputeData requires all sampled subjects = 1 (not 1 and 2)
                
                impdat2 = ImputeData(FIT=fit2, DAT=adat2_tmp3, M=MIREP, MARG.EXP.FORMULA=MARG.EXP.FORMULA, Sampled='samplek', 
                                   verbose=VERBOSE, EXPVAR='Xe')
                
                betas = vars = vector('list', length(impdat2))
                
                for(zzx in seq_along(impdat2)) {
                fits2 =  mm(as.formula(MODEL), dat=impdat2[[zzx]], id=id, lv.formula=LVMODEL, t.formula=TMODEL, iter.lim=ITERLIM, q=Q)
                betas[[zzx]] = c(fits2$beta, fits2$alpha)
                vars[[zzx]] = fits2$mod.cov
                }
                fit2 = MIcombine(betas, vars)
            }
            #### This part is rewritten by RT #### 
            optimal_mat[opt_des,] = get_opt(DAT=fit2)^(-1)
            #### This part is rewritten by RT ####
        }
          
        #### This part is rewritten by RT ####
        poss_des2$opt = NA
        poss_des2$opt = optimal_mat[,match(TARGET, c('Xei', 'Xe', 'XeXei'))]
        
        if (N2FIXED) {
            des_fit = lm(opt~(n20+n22+I(n20^2)+I(n22^2))^2, data=poss_des2)
            new_poss_des2 = expand.grid(n20=seq(0, WAVED2), n21=seq(0, WAVED2), n22=seq(0, WAVED2))
            new_poss_des2 = new_poss_des2[rowSums(new_poss_des2) == WAVED2,]
            for (i in 1:3) {
                new_poss_des2 = new_poss_des2[which(new_poss_des2[,i] < REMAINING[i]),]
            }
            preds = predict(des_fit, newdata=new_poss_des2)
            
            best_des = new_poss_des2[which(preds == max(preds)),]
            # best_des = poss_des2[which.max(poss_des2$opt),c("n20", "n21", "n22")]            
        } else {
            des_fit = lm(opt~(n21+I(n21^2))^2, data=poss_des2)
            new_poss_des2 = data.frame(n20 = 0, n21 = seq(min(poss_des2$n21), max(poss_des2$n21)), n22 = 0)
            preds = predict(des_fit, newdata=new_poss_des2)
            best_des = new_poss_des2[which(preds >= TARGETVAR^(-1))[1],]
        }

        #### This part is rewritten by RT ####
        
        best_des_mat[mxx,] = as.numeric(best_des)
    }
    best_des = ceiling(apply(best_des_mat, 2, mean, na.rm=TRUE))
    attr(best_des, 'hx') = best_des_mat
    best_des
}

twowave_aODS = function(DATA=NULL, SEEDS=1, YSEED=NULL, NLOW=4, NHIGH=6, N=5000, PXE=0.25,
                        BETAZ=data.frame('Z'=1), GAMMA=2, SIGMA=NULL, OMEGAZ.Xe=1, 
                        BETAM=data.frame('Int'=-1.50, 'time'=-0.25, 'Xe'=1, 'Xe.i'=.25), 
                        LVMODEL=NULL, TMODEL=~1, MODEL='Y~time+Xe+Xe.i+Z',
                        MARG.EXP.FORMULA='Xe~Z', METHOD='ACML',
                        N2FIXED=TRUE, WAVED1=c(25, 50, 25), WAVED2=400, 
                        TARGET='Xei', TARGETVAR=NULL, 
                        MIREP=1, GRIDREP=1, GRIDINC=20, INC=50, RNGp=c(0.05, 0.95),
                        Q=2, ITERLIM=500, VERBOSE=FALSE, BREP=1, MIREP.final=50) {
  
  # Arguments
  # DATA = if NULL, then simulate data using:
  #   SEEDS, YSEED = seed (YSEED != NULL if we want a force fixing outcome)
  #   N, (NLOW, NHIGH) = sample size, and range for number of visits per subject  
  #   PXE = probab ility of exposure
  #   BETAM, BETAZ = marginal parameter values (M), confounder (Z)
  #   OMEGAZ.Xe = mean difference in Z between for Xe=1, Xe=0
  #   LVMODEL, TMODEL = dependence model specification: latent variable (LV), transition (T) models
  #   MODEL = marginal model specification
  #   MARG.EXP.FORMULA = marginal exposure model specification
  # METHOD = 'ACML' or 'MI'
  # N2FIXED = if TRUE, perform adaptive design ODS; else perform adaptive sample size ODS
  # WAVED1 = ODS design for wave one
  # WAVED2 = overall size of wave two ODS (n_20 + n_21 + n_22)
  # TARGET = object of inference (either 'Xei','Xe','XeXei')
  # TARGETVAR = if !NULL, then target variance for Var(Xe.i)
  # MIREP = number of imputed data sets for the MI approach in interim design evaluation (see ImputeData) 
  # GRIDREP =  number of replications when computing N2 via uniform sampling (i.e., GRIDREP = 2 will replicate each guess twice)
  # INC = increment size of ODS design (i.e., 10 implies D[10,50,10], D[20,40,10], etc)
  # RNGp = precentiles used to restrict design space for adaptive sample size designs (default=0.05, 0.95 which 
  #        means we will restrict the grid search to n2 central 90 percent of possible values)
  # Q = number of quadrature points for numerical integration
  # ITERLIM = iteration limit
  # VERBOSE = if TRUE, print fitting details
  # BREP = Number of complet cohort datasets generated in the interim design evaluation
  # MIREP.final = number of imputed data sets for the MI approach in the final data analysis (see ImputeData)
  
  tA = Sys.time()
  
  # Generate data
  set.seed(SEEDS)
  
  if(is.null(DATA)) { 
    if(VERBOSE) {
      cat('\n\r',date(),'\n\r','Generating data')
    }
    
    XI = gen_XI(SEED=SEEDS, NLOW, NHIGH, N, PXE, OMEGAZ.Xe)
    if(is.null(YSEED)) {
      YSEED = SEEDS
    }
    sdat = gen_data(SEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, BETAZ=BETAZ,BETAM=BETAM, GAMMA=GAMMA, 
                    SIGMA=SIGMA, Q=Q, XI=XI, MEANFORMULA=MODEL, LVFORMULA=LVMODEL, TFORMULA=TMODEL)
    sdat = split(sdat, sdat$id)
    fcdat = data.frame(do.call(rbind, lapply(sdat, as.matrix)))
    
    fc_ss_dist = table(sapply(sdat, function(ZZ) ZZ$ss[1]))
    
    adat = fcdat
  } else {
    adat = fcdat = DATA
    sdat = split(adat, adat$id)
  }
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing full cohort analysis')
  # Perform full cohort analysis
  fc_fit = mm(as.formula(MODEL), data=adat, lv.formula=LVMODEL, t.formula=TMODEL, id=id, q=Q, iter.lim=ITERLIM)
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing wave 1 analysis')
  # Perform wave 1 analysis
  tmp_ss = idsXss(DAT=sdat)
  
  # Frequencies by sampling stratum
  n_ss = sapply(tmp_ss, length)
  
  # Random sample
  adat = adat[order(adat$id, adat$time),]
  adat_rs = adat
  ids = unique(adat$id)
  samp = sum(WAVED1)+WAVED2
  ids_rs = sample(ids, samp, replace=FALSE)
  rs_sample = (seq(nrow(adat)) %in% which(adat$id %in% ids_rs))*1
  adat_rs$sample = rs_sample
  
  rs.tlv = mm(as.formula(MODEL), dat=adat_rs[adat_rs$sample==1,], id=id, q=Q, lv.formula=LVMODEL, t.formula=TMODEL, iter.lim=ITERLIM) 
  rs_ml = list('coef'=coef(rs.tlv), 'vcov'=vcov(rs.tlv), 'convergence'=rs.tlv$control['convergence_code'])
  
  #### This part is rewritten by RT ####
  ### Single wave fixed design start ###
  adat_single = adat
  
  adat_single$Xe.o = adat_single$Xe
  adat_single$Xe.i.o = adat_single$Xe.i
  adat_single$Xe = adat_single$Xe.i = NA
  adat_single$samplek = 0
  adat_single$sprob3 = adat_single$sprob2 = adat_single$sprob1 = NA
  
  wave1.p = WAVED1/sum(WAVED1)
  wave.single = ceiling((sum(WAVED1)+WAVED2)*wave1.p)
  
  d_single = apply(cbind(wave.single, attr(tmp_ss,'freq')$ss), 1, min)
  s_single_info = get_sample(DESIGN=d_single, IDSXSS=tmp_ss, RSEXACT=FALSE)
  
  s_single_id = sort(unique(unlist(s_single_info$idsXss_samp)))
  s_single_idX = which(adat_single$id%in%s_single_id)
  adat_single$samplek[s_single_idX] = 1
  adat_single[s_single_idX, c('sprob1','sprob2','sprob3')] = matrix(s_single_info$sp, nrow=1)[rep(1, length(s_single_idX)),]
  adat_single$Xe[s_single_idX] = adat_single$Xe.o[s_single_idX]
  adat_single$Xe.i[s_single_idX] = adat_single$Xe.i.o[s_single_idX]
  
  adat_single_acml = adat_single[which(adat_single$samplek == 1),]
  # ACML
  fit_single = mm(as.formula(MODEL), dat=adat_single_acml, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                  cond.like=TRUE, samp.probs=s_single_info$sp, iter.lim=ITERLIM, q=Q, return_args=TRUE)
  fit_single_res = list('coef'=coef(fit_single), 'vcov'=vcov(fit_single), 'convergence'=fit_single$control['convergence_code'])
  
  # MI
  if(METHOD=='MI') {
      tmpXe_mi_single = ImputeData(FIT=fit_single, DAT=adat_single, M=MIREP.final, MARG.EXP.FORMULA=MARG.EXP.FORMULA, Sampled='samplek', verbose=VERBOSE, EXPVAR='Xe')
      
      fits_single = lapply(tmpXe_mi_single, 
                           function (ZZ) {
                               mm(as.formula(MODEL), dat=ZZ, id=id, lv.formula=LVMODEL, t.formula=TMODEL, iter.lim=ITERLIM, q=Q)
                               }
                           )
      
      betas_single = lapply(fits_single, function(ZZ) c(ZZ$beta, ZZ$alpha))
      vars_single = lapply(fits_single, function(ZZ) ZZ$mod.cov)
      
      mi_vals_single = MIcombine(betas_single, vars_single)
      alpha_ests_single = c(grep('gamma', names(coef(mi_vals_single))), grep('sigma', names(coef(mi_vals_single))))
      beta_ests_single = seq(length(coef(mi_vals_single)))[-alpha_ests_single]
      mi_converge_single = table(c(1:5, sapply(fits_single, function(ZZ) ZZ$control['convergence_code'])))-1
      fit_single_res = list('coef'=list('beta'=coef(mi_vals_single)[beta_ests_single], 'alpha'=coef(mi_vals_single)[alpha_ests_single]),
                            'vcov'=list('beta'=mi_vals_single$variance[beta_ests_single,beta_ests_single],
                                        'alpha'=mi_vals_single$variance[alpha_ests_single,alpha_ests_single]),
                            'convergence'=mi_converge_single)
  }
  ### Single wave fixed design end ###
  #### This part is rewritten by RT ####
  
  # Create analysis data set for final analysis (fcdata augmented with sampling probabilities)
  adat$Xe.o = adat$Xe                               # Copy of true (original) Xe value
  adat$Xe.i.o = adat$Xe.i                             # Copy of original Xe.i value
  adat$Xe = adat$Xe.i = NA                       # Replace Xe, Xe.i with NA 
  adat$samplek = 0                                     # numeric value noting stage k value when sampled
  adat$sprob3 = adat$sprob2 = adat$sprob1 = NA      # sampling probabilities for strata 1-3
  
  # Create wave 1 sampling info and perform analysis
  w1d = apply(cbind(WAVED1, attr(tmp_ss,'freq')$ss), 1, min)
  s1info = get_sample(DESIGN=w1d, IDSXSS=tmp_ss, RSEXACT=FALSE)
  
  # Update adat with wave 1 information
  s1id = sort(unique(unlist(s1info$idsXss_samp)))
  s1idX = which(adat$id%in%s1id)
  adat$samplek[s1idX] = 1
  adat[s1idX, c('sprob1','sprob2','sprob3')] = matrix(s1info$sp, nrow=1)[rep(1, length(s1idX)),]
  adat$Xe[s1idX] = adat$Xe.o[s1idX]
  adat$Xe.i[s1idX] = adat$Xe.i.o[s1idX]
  
  # Perform wave 1 analysis
  adat1 = adat[(adat$samplek==1), ] #is.na(adat$samplek) 
  fit1 = mm(as.formula(MODEL), dat=adat1, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
            cond.like=TRUE, samp.probs=s1info$sp, iter.lim=ITERLIM, q=Q, return_args=TRUE)
  fit1_res = list('coef'=coef(fit1), 'vcov'=vcov(fit1), 'convergence'=fit1$control['convergence_code'])
  
  # Create temporary data set
  adat2 = adat
  
  tmp_ss_k = s1info$idsXss_remain
  tmp_ss_kn = sapply(tmp_ss_k, length)
  prob_sampled = matrix(s1info$sp, nrow=1)
  colnames(prob_sampled) = c('sp1', 'sp2', 'sp3')
  
  # Perform IDEA to find wave 2 design or sample size
  idea_res = idea(iFIT=fit1, iDAT=adat, iS1INFO=s1info, REMAINING=tmp_ss_kn, MARG.EXP.FORMULA=MARG.EXP.FORMULA, N2FIXED=N2FIXED, 
                  METHOD=METHOD, WAVED2=WAVED2, INC=INC, MODEL=MODEL, LVMODEL=LVMODEL, TMODEL=TMODEL, ITERLIM=ITERLIM, 
                  MIREP=MIREP, TARGETVAR=TARGETVAR, TARGET=TARGET, VERBOSE=VERBOSE, RNGp=RNGp, BREP=BREP, Q=Q, GRIDINC=GRIDINC)
  
  design_info = list('WAVED1'=WAVED1, 'WAVED2'=WAVED2, 'best'=idea_res) 

  # Perform second wave design
  adat2_tmp = adat2
  adat2_tmp$Xe = adat2_tmp$Xe.o
  adat2_tmp$Xe.i = adat2_tmp$Xe.i.o
    
  s2d = idea_res
  s2info = get_sample(DESIGN=s2d, IDSXSS=tmp_ss_k, RSEXACT=FALSE)
  dat_s = split(adat2_tmp, adat2_tmp$id)
  dat_ix = which(adat2_tmp$id %in% unlist(s2info$idsXss_samp))
  adat2_tmp$samplek[dat_ix] = 2
  adat2_tmp[dat_ix,c('sprob1', 'sprob2', 'sprob3')] = matrix(s2info$sp, nrow=1)[rep(1, length(dat_ix)),]
    
  # Estimate Stage K sampling probabilities
  prob_sampled = do.call(rbind, lapply(split(adat2_tmp, adat2_tmp$samplek), function(ZZ) ZZ[1,c('sprob1', 'sprob2', 'sprob3')]))
  if (any(rownames(prob_sampled) == 0)) prob_sampled = prob_sampled[-which(rownames(prob_sampled) == 0),]
  prob_not_sampled = apply(1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
  prob_not_sampled = rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])
  prob_sampled = prob_sampled*prob_not_sampled
  prob_sampled$samplek = as.numeric(rownames(prob_sampled))
  colnames(prob_sampled)[1:3] = c('sp1', 'sp2', 'sp3')
    
  adat2_tmp = merge(adat2_tmp, prob_sampled, by='samplek', all.x=TRUE)
  adat2_tmp_mi = adat2_tmp
  adat2_tmp_mi[,c('sprob1', 'sprob2', 'sprob3')] = adat2_tmp_mi[,c('sp1','sp2','sp3')] 
  adat2_tmp_mi$samplek[adat2_tmp_mi$samplek>0] = 1
    
  adat2_tmp = adat2_tmp[adat2_tmp$samplek>0,]
  adat2_tmp = adat2_tmp[order(adat2_tmp$id, adat2_tmp$time),]
    
  adat2_tmp_mi = adat2_tmp_mi[order(adat2_tmp_mi$id, adat2_tmp_mi$time),]
  
  fit2 = mm(as.formula(MODEL), dat=adat2_tmp, id=id, lv.formula=LVMODEL, t.formula=TMODEL,
            cond.like=TRUE, samp.probs=as.matrix(adat2_tmp[,c('sp1','sp2','sp3')]),
            iter.lim=ITERLIM, q=Q, return_args=TRUE)
  
  best_ests = list('coef'=coef(fit2), 'vcov'=vcov(fit2), 'convergence'=fit2$control['convergence_code'])
  
  # MI
  if(METHOD=='MI') {
    # Update Xe, Xe.i
    adat2_tmp_mi$Xe[adat2_tmp_mi$samplek == 0] = NA
    adat2_tmp_mi$Xe.i[adat2_tmp_mi$samplek == 0] = NA
    
    tmpXe_mi = ImputeData(FIT=fit2, DAT=adat2_tmp_mi, M=MIREP.final, MARG.EXP.FORMULA=MARG.EXP.FORMULA, Sampled='samplek', verbose=VERBOSE, EXPVAR='Xe')
    
    fits = lapply(tmpXe_mi, function(ZZ) { 
      mm(as.formula(MODEL), dat=ZZ, id=id, lv.formula=LVMODEL, t.formula=TMODEL, iter.lim=ITERLIM, q=Q)
    })
    
    betas = lapply(fits, function(ZZ) c(ZZ$beta, ZZ$alpha))
    vars = lapply(fits, function(ZZ) ZZ$mod.cov)
    
    mi_vals = MIcombine(betas, vars)
    alpha_ests = c(grep('gamma', names(coef(mi_vals))), grep('sigma', names(coef(mi_vals))))
    beta_ests = seq(length(coef(mi_vals)))[-alpha_ests]
    mi_converge = table(c(1:5, sapply(fits, function(ZZ) ZZ$control['convergence_code'])))-1
    best_ests = list('coef'=list('beta'=coef(mi_vals)[beta_ests], 'alpha'=coef(mi_vals)[alpha_ests]),
                     'vcov'=list('beta'=mi_vals$variance[beta_ests,beta_ests], 'alpha'=mi_vals$variance[alpha_ests,alpha_ests]),
                     'convergence'=mi_converge)
  }

  tB = Sys.time()
  
  outs = list('seed'=SEEDS, 'time'=difftime(tB, tA, units='mins'), 'ss'=attr(tmp_ss, 'freq')$ss, 
              'fc'=list('coef'=coef(fc_fit), 'vcov'=vcov(fc_fit), 'modcov'=fc_fit$mod.cov),
              'rs'=rs_ml, 'acml_1'=fit1_res, 'ests'=best_ests,
              'design_info'=design_info, 'method'=METHOD, 'single'=fit_single_res)
  outs
}
