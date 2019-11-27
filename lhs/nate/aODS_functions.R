
cluster.summary <- function( id, x, fun ){ 
  xlist <- split( x, id )
  nj    <- unlist( lapply( xlist, length ) )
  xj    <- unlist( lapply( xlist, fun) )
  xsummary <- rep( xj, nj )
  xsummary
}
trace <- function(ZZ) sum(diag(ZZ))
get_opt <- function(DAT) {
  var_gt <- vcov(DAT)$beta['Xe.i','Xe.i']  
  var_g  <- vcov(DAT)$beta['Xe','Xe']
  d_opt  <- sqrt(det( solve(vcov(DAT)$beta)[c('Xe','Xe.i'),c('Xe','Xe.i')] ))
  c(var_gt, var_g, d_opt)
}
idsXss     <- function(DAT) {
  # DAT = list of subject-level data (output from gen_data function)
  
  # Output: out (list of ids by each sampling stratum) and attr(out,'summ') which contains key frequencies
  #         idXss will be used when performing sampling to keep track of remaining ids
  
  # Generate ids by each sampling stratum
  N      <- length(DAT)
  #DAT_ss <- split(seq(N), sapply(DAT, function(ZZ) ZZ$ss[1]))
  DAT_ss <- split(sapply(DAT, function(ZZ) ZZ$id[1]), sapply(DAT, function(ZZ) ZZ$ss[1]))
  
  # Calculate frequencies in each sampling statrum and frequencies of sum(Y) and mean(Y)
  ss_sum  <- sapply(DAT_ss, length)
  sy_sum  <- table(sapply(DAT, function(ZZ) sum(ZZ$Y)))
  sy_mean <- table(sapply(DAT, function(ZZ) mean(ZZ$Y)))
  
  out <- DAT_ss
  attr(out, 'freq') <- list('ss'=ss_sum, 'sy'=sy_sum, 'sm'=sy_mean)
  out
}
get_sample <- function(DESIGN, IDSXSS, RSEXACT, REPLACE=FALSE) {
  # DESIGN  = vector of expected sample counts
  # IDSXSS  = list of ids by sampling stratum (length=3)
  # RSEXACT = indicator to force exact random sampling
  # REPLACE = indicator if sampling with replacement should be performed
  
  # Output: list that contains sampling probs, and lists of sampled ids by ss and remaining ids by ss
  
  tmp_design  <- as.list(DESIGN)
  blah_samp   <- vector('list',length(tmp_design))
  
  if( RSEXACT ) {
    blah <- do.call(rbind, mapply( function(AA,BB) list(cbind(AA,BB)), 
                                   AA=as.list( as.numeric(names(IDSXSS)) ),BB=IDSXSS))
    blah_tmp <-  data.frame( blah[ sort( sample(x=seq(nrow(blah)),round( sum(DESIGN),0) , replace=FALSE)) ,] )
    if(ncol(blah_tmp)==1) blah_tmp <- data.frame(t(blah_tmp))
    blah_samp     <- blah_tmp
    blah_samp[,1] <- factor(blah_samp[,1], levels=1:3)
    blah_samp     <- split(blah_samp$BB, blah_samp$AA)
  }
  
  tmp_samp    <- mapply( function(AA,BB,CC) {
    nn <- length(AA)
    if(nn < BB) BB <- nn
    sp <- ifelse(nn==0,0, BB/nn)
    sampled <- which( rbinom( nn, size=1, prob=sp )== 1 )
    if(length(sampled)>0 | !is.null(CC)) {
      if(is.null(CC)) {
        sampled <- AA[sampled]
      } else {
        sampled <- CC
      }
    }
    attr(sampled,'sp') <- sp
    list(sampled) }, AA=IDSXSS ,BB=tmp_design, CC=blah_samp)
  
  tmp_sp <- sapply( tmp_samp, function(ZZ) attr(ZZ,'sp') )
  
  tmp_remain <- mapply( function(AA,BB) {
    out <- AA
    if(length(BB)>0)  out <- AA[!(AA%in%BB)]
    out
  }, AA=IDSXSS, BB=tmp_samp)
  
  list('sp'=tmp_sp,'idsXss_samp'=tmp_samp,'idsXss_remain'=tmp_remain)
}
Odds2Prob  <- function(odds) odds/(1+odds) 
ods        <- function(SEED, DAT, DES1, DES2, MOD, STEPMAX,TS=FALSE,VERBOSE=FALSE,INITS=NULL,RETURNDAT=FALSE) {
  
  adat2 <- lapply( split(DAT, DAT$id), function(ZZ) { 
    ZZ$nobs   <- nrow(ZZ)
    ZZ$sumY   <- sum(ZZ$Y)
    ZZ$ss     <- ifelse(ZZ$sumY==ZZ$nobs, 3, ifelse(ZZ$sumY==0, 1, 2))
    ZZ$sample <- NA
    ZZ})
  
  tmp_ss <- idsXss(DAT=adat2) 
  DAT$sample <- NA
  
  # Obtain phase 1 ODS sample
  set.seed(SEED)
  s1d      <- apply( cbind( DES1, attr(tmp_ss,'freq')$ss ), 1, min )
  s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss)
  
  s1id     <- sort( unique(unlist( s1info$idsXss_samp )) )
  s1idX    <- which(DAT$id%in%s1id)
  DAT$sample[s1idX] <- 1
  
  # Obtain phase 2 ODS sample
  tmp_ss_2      <- s1info$idsXss_remain
  tmp_ss_2n     <- sapply(tmp_ss_2,length)
  prob_sampled  <- matrix( s1info$sp,nrow=1 )
  colnames(prob_sampled) <- c('sp1','sp2','sp3')
  # print(prob_sampled)
  if(!TS) { 
    #    ts.fit.lv <- mm(as.formula(MOD), dat=DAT[DAT$sample==1,], id=id, lv.formula=~1, t.formula=NULL, iter.lim=500,verbose=VERBOSE,
    #                    cond.like=TRUE, samp.probs=as.numeric(s1info$sp),step.max=STEPMAX)
    #    ts.fit.t <- mm(as.formula(MOD), dat=DAT[DAT$sample==1,], id=id, lv.formula=NULL, t.formula=~1, iter.lim=500,verbose=VERBOSE,
    #                   cond.like=TRUE, samp.probs=as.numeric(s1info$sp),step.max=STEPMAX)
    ts.fit.tlv <- mm(as.formula(MOD), dat=DAT[DAT$sample==1,], id=id, lv.formula=~1, t.formula=~1, iter.lim=500,verbose=VERBOSE,
                     cond.like=TRUE, samp.probs=as.numeric(s1info$sp),step.max=STEPMAX,inits=INITS)
    blah <- ts.fit.tlv #list(ts.fit.lv, ts.fit.t, ts.fit.tlv)
    if(RETURNDAT) { 
      attr(blah,'data') <- DAT
      attr(blah,'sp')   <- s1info$sp
    }
  }
  if(TS){
    s2d      <- apply( cbind( DES2, sapply(tmp_ss_2,length)), 1, min )
    s2info   <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_2) 
    s2idX    <- which(DAT$id %in% unlist(s2info$idsXss_samp))
    DAT$sample[s2idX] <- 2 
    
    # Estimate phase 2 sampling probabilities
    prob_sampled        <- data.frame( rbind(s1info$sp,s2info$sp) )
    prob_not_sampled    <- apply( 1-prob_sampled, 2, cumprod)                  # cumulative prob of not being sampled
    prob_not_sampled    <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
    prob_sampled        <- prob_sampled * prob_not_sampled 
    prob_sampled$sample <- c(1,2)
    colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3')
    
    # Update newdata with sampling probabilities
    DAT2  <- merge(DAT, prob_sampled, by='sample',all.x=TRUE)
    DAT2  <- DAT2[ order(DAT2$id, DAT2$visit1), ]
    
    # Perform two-stage fixed ODS design
    DAT3  <- DAT2[!is.na(DAT2$sample), ]
    
    # ts.fit.lv <- mm(as.formula(MOD), dat=DAT3, id=id, lv.formula=~1, iter.lim=500,verbose=VERBOSE,
    #                  cond.like=TRUE, samp.probs=as.matrix(DAT3[,c('sp1','sp2','sp3')]),step.max=STEPMAX)
    #  ts.fit.t <- mm(as.formula(MOD), dat=DAT3, id=id, t.formula=~1, iter.lim=500,verbose=VERBOSE,
    #                 cond.like=TRUE, samp.probs=as.matrix(DAT3[,c('sp1','sp2','sp3')]),step.max=STEPMAX)
    ts.fit.tlv <- mm(as.formula(MOD), dat=DAT3, id=id, lv.formula=~1, t.formula=~1, iter.lim=500,verbose=VERBOSE,
                     cond.like=TRUE, samp.probs=as.matrix(DAT3[,c('sp1','sp2','sp3')]),step.max=STEPMAX,inits=INITS)
    blah <- ts.fit.tlv#list(ts.fit.lv, ts.fit.t, ts.fit.tlv)
    if(RETURNDAT) attr(blah,'data') <- DAT2
  }
  blah
}
prep_data <- function(mean.formula, lv.formula, t.formula, id, data, inits,
                      samp.probs, samp.probi, offset) {
  
  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  
  y  = model.response(mean.f,'numeric') 
  uy = unique(y)
  
  x  = model.matrix(mean.formula,mean.f)
  
  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data)) 
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))
  
  if(!is.matrix(samp.probs)) samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)
  if(is.null(samp.probi)) samp.probi = matrix(1,nrow=length(y),ncol=1) 
  if(!is.null(inits)) inits = unlist(inits)
  
  if(is.null(inits)) {
    inits = c(glm(mean.formula,family='binomial',data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x            = x[,-c(omit_dup_col)]
      inits        = inits[-c(omit_dup_col)]
    }
  }
  
  if(is.null(offset)) offset = rep(0, length(y))
  
  x.t = cbind(x.t, NULL)
  x.lv = cbind(x.lv, NULL)
  
  if ( length(inits) != ncol(x) + ncol(x.t) + ncol(x.lv) ) { stop("Parameter length incongruous with X, Xgam, and Xsig")}
  paramlengths <- c(ncol(x), ncol(x.t), ncol(x.lv))
  
  ## Lagged response
  Ylag = rep(0, length(y))
  for(i in 2:length(y)) { if(id[i]==id[i-1]) Ylag[i]<-y[i-1] }
  
  id.tmp        = split(id, id)
  X.tmp         = split(x,id)
  Y.tmp         = split(y,id)
  Ylag.tmp      = split(Ylag,id)
  Xgam.tmp      = split(x.t,id)
  Xsig.tmp      = split(x.lv,id)
  SampProbi.tmp = split(samp.probi,id)
  SampProbs.tmp = split(samp.probs,id)
  offset.tmp    = split(offset,id)
  
  subjectData <- vector('list', length=length(unique(id)))
  subjectData <- list()
  uid <- as.character(unique(id))
  for(j in seq(along=uid)){
    i <- uid[j]
    subjectData[[j]] <- list(id=as.character(unique(id.tmp[[i]])), 
                             X=matrix(X.tmp[[i]], ncol=ncol(x)), 
                             Y=as.double(Y.tmp[[i]]),
                             Ylag=as.double(Ylag.tmp[[i]]),
                             Xgam=matrix(Xgam.tmp[[i]], ncol=ncol(x.t)),
                             Xsig=matrix(Xsig.tmp[[i]], ncol=ncol(x.lv)),
                             SampProbi=unique(SampProbi.tmp[[i]]),
                             SampProbs=matrix(SampProbs.tmp[[i]], ncol=ncol(samp.probs))[1,],
                             Offset=as.double(offset.tmp[[i]]))
  }
  names(subjectData) <- uid
  
  beta_nms  <- colnames(x)
  alpha_nms <- c( if(!is.null(t.formula)){paste('gamma',colnames(x.t),sep=':')}, 
                  if(!is.null(lv.formula)){paste('log(sigma)',colnames(x.lv),sep=':')}) 
  
  list(params=inits, paramlengths=paramlengths, subjectData=subjectData, samp.probi=tapply(samp.probi, id, unique),
       nms=list(beta_nms, alpha_nms))
}
indir_imp_data  <- function(FIT, DAT, M=30, ROB=FALSE, MARG.EXP.FORMULA, verbose=FALSE,EXPVAR=EXPVAR, Sampled=Sampled, samp.probs=c(1,1,1)) {
  # FIT - MMLongit object
  # DAT - complete data set with sampled flag (Sampled)
  # M   - # of imputations
  # ROB - indicator, if FALSE use model-based covariance
  # FORM - mean model formula
  # MARG.EXP.FORMULA - marginal exposure model formula
  # EXPVAR - exposure / expensive variable
  # Sampled - set to 1 if included in the sample and 0 otherwise
  imp_data <- vector('list', M)
  
  for(ix in seq(M)) {
    # Step 0: Update date files
    # Create versions of DAT with updated Xe values
    # All 0 -> fev.rs6657868!="CC" (TT or TC), All 1 -> fev.rs6657868=="CC"
    DAT_Xe1 <- DAT_Xe0 <- DAT
    DAT_Xe1[,EXPVAR] <- 1
    DAT_Xe0[,EXPVAR] <- 0
    
    # Non-sampled subjects, remove Xe
    #DAT_S0    <- DAT[DAT$Sampled==0,]
    #DAT_S0[,EXPVAR] <- NA
    #DAT_S0_v0 <- DAT_S0[DAT_S0$year==0, ]
    # DAT_S0_Xe1 <- DAT_Xe1[DAT_Xe1$S==0, ]
    #DAT_S0_Xe0 <- DAT_Xe0[DAT_Xe0$S==0, ]
    
    DAT_S0     <- DAT[DAT[,Sampled]==0,]
    DAT_S0_Xe1 <- DAT_Xe1[DAT_Xe1[,Sampled]==0, ]
    DAT_S0_Xe0 <- DAT_Xe0[DAT_Xe0[,Sampled]==0, ]
    dup        <- duplicated(DAT_S0$id)
    DAT_S0_v0  <- DAT_S0[!dup, ]
    
    # Sampled subjects
    #DAT_S1    <- DAT[DAT$S==1,]
    #DAT_S1_v0 <- DAT[DAT$S==1 & DAT$year==0,]
    #DAT_S1_Xe1 <- DAT_Xe1[DAT_Xe1$S==1, ]
    #DAT_S1_Xe0 <- DAT_Xe0[DAT_Xe0$S==1, ]
    
    DAT_S1     <- DAT[DAT[,Sampled]==1,]
    DAT_S1_Xe1 <- DAT_Xe1[DAT_Xe1[,Sampled]==1, ]
    DAT_S1_Xe0 <- DAT_Xe0[DAT_Xe0[,Sampled]==1, ]
    dup        <- duplicated(DAT_S1$id)
    DAT_S1_v0  <- DAT_S1[!dup,]
    
    # Step 1: extract estimates from FIT and draw m^th theta from MVN
    theta      <- c(FIT$beta, FIT$alpha)
    if(ROB) {
      cov_theta <- FIT$rob.cov
    } else {
      cov_theta <- FIT$mod.cov
    }
    
    theta.m <- mvrnorm(n = 1, mu=theta, Sigma=cov_theta)
    ##
    ##
    ## allow this to work for a random intercept only model
    ##
    ##
    #bbb <- length(theta.m)
    #if (is.na(match("gamma:(Intercept)", names(FIT$alpha)))) theta.m <- c(theta.m[1:(bbb-1)], 0, theta.m[bbb])
    
    # Step 2: Using theta and calculate liac for sampled subjects!
    # liac_S1_Xe0 <- mm_LogLScore(mean.formula=FORM, iter.lim=500, step.max=10,
    #                             lv.formula= ~1, t.formula = ~1, id=id, data=DAT_S1_Xe0, 
    #                             cond.like=TRUE,samp.probs=samp.probs, inits=theta.m) 
    # liac_S1_Xe1 <- mm_LogLScore(mean.formula=FORM, iter.lim=500, step.max=10,
    #                             lv.formula= ~1, t.formula = ~1, id=id, data=DAT_S1_Xe1, 
    #                             cond.like=TRUE,samp.probs=samp.probs, inits=theta.m) 
    prep_input <- attr(FIT,'args')
    tmp_S1_Xe1 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, 
                            t.formula=prep_input$t.formula, 
                            id=DAT_S1_Xe1[,prep_input$id], data=DAT_S1_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi, 
                            offset=prep_input$offset)$subjectData
    tmp_S1_Xe0 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula, 
                            id=DAT_S1_Xe0[,prep_input$id], data=DAT_S1_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi, 
                            offset=prep_input$offset)$subjectData
    #print("blah")
    # Step 3: Perform offsetted logistic regression on sampled subset 
    # Note: liac are already on the log scale, so log(a/b) = loga-logb
    
    LLSC_1 <- MMLB:::LogLScoreCalc(params = theta.m, subjectData = tmp_S1_Xe1, 
                                   Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z = FIT$LLSC_args$Z, 
                                   ParamLengths = FIT$LLSC_args$ParamLengths,
                                   CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB) 
    LLSC_0 <- MMLB:::LogLScoreCalc(params = theta.m, subjectData = tmp_S1_Xe0, 
                                   Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z =  FIT$LLSC_args$Z, 
                                   ParamLengths = FIT$LLSC_args$ParamLengths,
                                   CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    #offset_log <- attr(liac_S1_Xe1,'liac')[,2]-attr(liac_S1_Xe0,'liac')[,2]
    
    logAC_1 <- attr(LLSC_1, "ACSubj") # Ascertainment corrections already log transformed
    logAC_0 <- attr(LLSC_0, "ACSubj")
    DAT_S1_v0$offset <- logAC_1-logAC_0
    #print("blah")
    # Perform offsetted logistic regression [Note, only using first visit info - time-invariant exposure]
    fit.exp <- glm(MARG.EXP.FORMULA, offset=offset, data=DAT_S1_v0,family=binomial())
    
    # lp  <- predict(fit.exp)
    # lp2 <- model.matrix(fit.exp) %*% cbind(coef(fit.exp),NULL)+offset_log
    # predict.glm does account for offset term!
    
    # Step 4: Extract alpha estimates and sampled m^th alpha from MVN
    alpha      <- coef(fit.exp)
    cov_alpha  <- vcov(fit.exp)
    
    alpha.m <- mvrnorm(n = 1, mu=alpha, Sigma=cov_alpha)
    
    # Step 5: Compute predicted probabilies using alpha.m and the design matrix of non-sampled subjects
    # liac_S0_Xe0 <- mm_LogLScore(mean.formula=FORM, iter.lim=500, step.max=10,
    #                             lv.formula= ~1, t.formula = ~1, id=id, data=DAT_S0_Xe0, 
    #                             cond.like=TRUE,samp.probs=samp.probs, inits=theta.m) 
    # liac_S0_Xe1 <- mm_LogLScore(mean.formula=FORM, iter.lim=500, step.max=10,
    #                             lv.formula= ~1, t.formula = ~1, id=id, data=DAT_S0_Xe1, 
    #                             cond.like=TRUE,samp.probs=samp.probs, inits=theta.m) 
    
    nr <- nrow(DAT_S0_Xe1)
    
    tmp_S0_Xe1 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, 
                            t.formula=prep_input$t.formula, 
                            id=DAT_S0_Xe1[,prep_input$id], data=DAT_S0_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),], 
                            samp.probi=prep_input$samp.probi[rep(1,nr)], 
                            offset=prep_input$offset[rep(1,nr)])$subjectData
    tmp_S0_Xe0 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula, 
                            id=DAT_S0_Xe0[,prep_input$id], data=DAT_S0_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),], samp.probi=prep_input$samp.probi[rep(1,nr)], 
                            offset=prep_input$offset[rep(1,nr)])$subjectData
    
    LLSC_1 <- MMLB:::LogLScoreCalc(params = theta.m, subjectData = tmp_S0_Xe1, 
                                   Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z = FIT$LLSC_args$Z, 
                                   ParamLengths = FIT$LLSC_args$ParamLengths,
                                   CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB) 
    LLSC_0 <- MMLB:::LogLScoreCalc(params = theta.m, subjectData = tmp_S0_Xe0, 
                                   Q = FIT$LLSC_args$Q, W = FIT$LLSC_args$W, Z =  FIT$LLSC_args$Z, 
                                   ParamLengths = FIT$LLSC_args$ParamLengths,
                                   CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    
    logAC_1 <- attr(LLSC_1, "ACSubj") # Ascertainment corrections already log transformed
    logAC_0 <- attr(LLSC_0, "ACSubj")
    
    offset_log_S0 <- logAC_1-logAC_0
    
    # # Create 'junk' model.matrix for non-sampled subjects
    # tmp <- DAT_S0_v0
    # tmp[,EXPVAR] <- sample(c(0,1),size=nrow(tmp),replace=TRUE)
    # fit.exp.tmp <- glm(FORM2, data=tmp, family=binomial())
    # mm_tmp <- model.matrix(fit.exp.tmp)
    
    # Need to create a temporary data for non-sampled subjects in order to create the
    # model.matrix that is used to calculate marginal exposure odds
    tmp          <- DAT_S0_v0
    tmp[,EXPVAR] <- sample(c(0,1),size=nrow(tmp),replace=TRUE)
    mm_tmp       <- model.matrix(MARG.EXP.FORMULA, tmp)
    
    # Compute the conditional exposure odds in non-sampled subjects
    # Notice that all conditional odds are equal: [Xe | Y, Xo, S=1] = [Xe | Y, Xo, S=0] = [Xe | Y, Xo] since we sampled based on Y
    marg.odds.Xe.Xo.S1   <- exp(mm_tmp %*% alpha.m)*exp(offset_log_S0)  
    LR.Y.Xe.Xo.S1        <- exp( attr(LLSC_1, "LogLikeSubj")-attr(LLSC_0, "LogLikeSubj"))
    cond.odds.Xe.Y.Xo.S0 <-  LR.Y.Xe.Xo.S1*marg.odds.Xe.Xo.S1  
    
    # Impute Xe for non-sampled subjects 
    # Convert odds to probability, log(odds) = log(1/[1-p]) ---> p = odds/(1+odds)
    prXe1     <- Odds2Prob(cond.odds.Xe.Y.Xo.S0)
    XeS0      <- rbinom(nrow(tmp), size=1, prob=prXe1)
    
    # Add XeS0 info to DAT_S0_v0
    # DAT_S0_v0$snp <- XeS0
    DAT_S0_v0[,EXPVAR] <- XeS0
    DAT_S0 <- DAT_S0[,-which(colnames(DAT_S0)==EXPVAR)]
    DAT_S0 <- merge(DAT_S0, DAT_S0_v0[,c('id',EXPVAR)],by='id')
    DAT_imp <- rbind(DAT_S1, DAT_S0[,colnames(DAT_S1)])
    
    # Store imputation data
    imp_data[[ix]] <- DAT_imp
    if(verbose) cat('\r',ix)
  } # end of for ix loop
  
  imp_data
}
do.iMI   <- function(FIT, DAT, M, MARG.EXP.FORMULA,VERBOSE,EXPVAR,SAMPLED) {
  imp_dat <- indir_imp_data(FIT=FIT, DAT=DAT, M=M, MARG.EXP.FORMULA=MARG.EXP.FORMULA,verbose=VERBOSE, 
                            EXPVAR=EXPVAR,Sampled=SAMPLED)
  imp_res <- lapply( imp_dat, function(ZZ) {
    imp_fit <- mm(as.formula(Y_mod), dat=ZZ, id=id, lv.formula=~1, t.formula=~1, iter.lim=1000)
    list(ests=unlist(coef(imp_fit)),vcovs=imp_fit$mod.cov)
  } )
  
  mi_combine <- MIcombine(lapply(imp_res, function(ZZ)ZZ$ests), lapply(imp_res, function(ZZ)ZZ$vcovs))
  #mi_est     <- mi_combine$coefficients
  #mi_ses     <- sqrt(diag(mi_combine$variance))
  #cbind('est'=mi_est, 'ses'=mi_ses)
  list(mi_combine$coefficients,mi_combine$variance)
}
expit    <- function(Z) {
  expZ <- exp(Z)
  expZ/(1+expZ)
}
gen_XI   <- function(SEED, NLOW, NHIGH, N, PXE, BETAZ) {
  # Arguments
  #
  # SEED  = seed for rng
  # NLOW  = minumum cluster size
  # NHIGH = maximum cluster size
  # N     = number of subjects to sample
  # PXE   = proprotion of Xe in population
  # BETAZ = vector of coefficients to create binary Z
  
  # Output: design matrix
  
  # Set seed
  set.seed(SEED)
  
  # Generate number of clusters per subjects
  nclust <- sample( seq(NLOW,NHIGH), N, replace=TRUE)
  if(NLOW==NHIGH) nclust <- rep(NLOW, N)
  
  # Generate key exposure (Xe) and a possibly related binary variable Z 
  Xe     <- rbinom(N,size=1,prob=PXE)
  # Z      <- rbinom(N,size=1,prob=expit(BETAZ[2]*Xe + BETAZ[1]) )
  Z_Xe0  <- round( rnorm(sum(Xe==0),0,1), 4 )
  Z_Xe1  <- round( rnorm(sum(Xe==1),1,1), 4 )
  Z      <- rep(NA, length(Xe))
  Z[Xe==0] <- Z_Xe0
  Z[Xe==1] <- Z_Xe1
  
  XI     <- vector('list', N)
  for ( i in seq(N) ){
    XI[[i]] <- cbind(i, rep(1, nclust[i]), seq(nclust[i])-1 , Xe[i], ( seq(nclust[i])-1 )*Xe[i], Z[i])
  }
  XI <- data.frame( do.call(rbind, XI) )
  colnames(XI) <- c('id','int','time','Xe','Xe.i','Z')
  XI
}
gen_data <- function(SEED, NLOW, NHIGH, N, PXE, MEANFORMULA, LVFORMULA, TFORMULA, BETAZ, BETAM, GAMMA, SIGMA, Q, XI=NULL) {
  # Set seed
  set.seed(SEED)
  
  if(is.null(XI)) stop('Specify XI')
  nclust <- as.numeric(table(XI$id))
  N      <- length(nclust)
  XI$Y   <- 0 # Just add junk column for Y, this will be rewritten later!
  
  if(is.character(MEANFORMULA)) MEANFORMULA <- as.formula(MEANFORMULA)
  if(is.character(LVFORMULA))   LVFORMULA   <- as.formula(LVFORMULA)  
  if(is.character(TFORMULA))    TFORMULA    <- as.formula(TFORMULA)
  
  # Generate Y
  ests <- BETAM
  if(!is.null(GAMMA)) ests <- c(ests, GAMMA)
  if(!is.null(SIGMA)) ests <- c(ests, log(SIGMA))
  
  if(is.data.frame(BETAM) | is.data.frame(BETAZ) ) {
    beta.mz <- c( as.numeric(BETAM), as.numeric(BETAZ) )
  } else {
    beta.mz <- c( BETAM, BETAZ )
  }
  
  sdat <- GenBinaryY(mean.formula=MEANFORMULA, lv.formula=LVFORMULA, t.formula=TFORMULA, 
                     beta = beta.mz, sigma = SIGMA, gamma = GAMMA, id=id, data=XI,q=Q, Yname = "Y") 
  sdat_s <- lapply( split(sdat, sdat$id), function(ZZ) {
    nY <- sum(ZZ$Y)
    nn <- nrow(ZZ)
    ZZ$ss <- (nY==0)*1+(nY>0 & nY<nn)*2+(nY==nn)*3
    ZZ
  })
  sdat <- do.call(rbind, sdat_s)
  attr(sdat,'true') <- ests
  sdat
}
prep_data      <- function(mean.formula, lv.formula, t.formula, id, data, inits,
                           samp.probs, samp.probi, offset) {
  
  mean.f = model.frame(mean.formula, data)
  mean.t = attr(mean.f, "terms")
  
  y  = model.response(mean.f,'numeric') 
  uy = unique(y)
  #print(summary(y))
  #print(summary(id))
  x  = model.matrix(mean.formula,mean.f)
  
  x.t = x.lv = matrix(0, ncol=1, nrow=length(y))
  if(!is.null(t.formula))   x.t  = model.matrix(t.formula,model.frame(t.formula, data)) 
  if(!is.null(lv.formula))  x.lv = model.matrix(lv.formula, model.frame(lv.formula, data))
  
  if(!is.matrix(samp.probs)) samp.probs = matrix(samp.probs, nrow = length(y), ncol = 3, byrow = TRUE)
  if(is.null(samp.probi)) samp.probi = matrix(1,nrow=length(y),ncol=1) 
  if(!is.null(inits)) inits = unlist(inits)
  
  if(is.null(inits)) {
    inits = c(glm(mean.formula,family='binomial',data=data)$coef, rep(1, ncol(x.t) + ncol(x.lv)))
    if(any(is.na(inits))) {
      omit_dup_col = which(is.na(inits))
      x            = x[,-c(omit_dup_col)]
      inits        = inits[-c(omit_dup_col)]
    }
  }
  
  if(is.null(offset)) offset = rep(0, length(y))
  
  x.t = cbind(x.t, NULL)
  x.lv = cbind(x.lv, NULL)
  
  nms.inits <- names(inits)
  ix.mean   <- seq(length(nms.inits))
  ix.t      <- grep('gamma',nms.inits)
  ix.lv     <- grep('sigma',nms.inits)
  ix.mean   <- ix.mean[-c(which(ix.mean %in% c(ix.t, ix.lv)))]
  
  if(length(ix.lv)==0) {
    inits.lv <- -5 
  } else {
    inits.lv <- inits[ix.lv]
  }
  
  if(length(ix.t)==0) {
    inits.t <- 0 
  } else {
    inits.t <- inits[ix.t]
  }
  inits <- c( inits[ix.mean], inits.t ,inits.lv )
  
  if ( length(inits) != ncol(x) + ncol(x.t) + ncol(x.lv) ) { stop("Parameter length incongruous with X, Xgam, and Xsig")}
  paramlengths <- c(ncol(x), ncol(x.t), ncol(x.lv))
  
  #print(length(id))
  #print(length(y))
  ## Lagged response
  Ylag = rep(0, length(y))
  for(i in 2:length(y)) { if(id[i]==id[i-1]) Ylag[i]<-y[i-1] }
  
  id.tmp        = split(id, id)
  X.tmp         = split(x,id)
  Y.tmp         = split(y,id)
  Ylag.tmp      = split(Ylag,id)
  Xgam.tmp      = split(x.t,id)
  Xsig.tmp      = split(x.lv,id)
  SampProbi.tmp = split(samp.probi,id)
  SampProbs.tmp = split(samp.probs,id)
  offset.tmp    = split(offset,id)
  
  subjectData <- vector('list', length=length(unique(id)))
  subjectData <- list()
  uid <- as.character(unique(id))
  for(j in seq(along=uid)){
    i <- uid[j]
    subjectData[[j]] <- list(id=as.character(unique(id.tmp[[i]])), 
                             X=matrix(X.tmp[[i]], ncol=ncol(x)), 
                             Y=as.double(Y.tmp[[i]]),
                             Ylag=as.double(Ylag.tmp[[i]]),
                             Xgam=matrix(Xgam.tmp[[i]], ncol=ncol(x.t)),
                             Xsig=matrix(Xsig.tmp[[i]], ncol=ncol(x.lv)),
                             SampProbi=unique(SampProbi.tmp[[i]]),
                             SampProbs=matrix(SampProbs.tmp[[i]], ncol=ncol(samp.probs))[1,],
                             Offset=as.double(offset.tmp[[i]]))
  }
  names(subjectData) <- uid
  
  beta_nms  <- colnames(x)
  alpha_nms <- c( if(!is.null(t.formula)){paste('gamma',colnames(x.t),sep=':')}, 
                  if(!is.null(lv.formula)){paste('log(sigma)',colnames(x.lv),sep=':')}) 
  
  list(params=inits, paramlengths=paramlengths, subjectData=subjectData, samp.probi=tapply(samp.probi, id, unique),
       nms=list(beta_nms, alpha_nms))
}
ImputeData     <- function(FIT, DAT, M=30, ROB=FALSE, MARG.EXP.FORMULA, verbose=FALSE,EXPVAR=EXPVAR, Sampled=Sampled) {
  # FIT - MMLongit object
  # DAT - complete data set with sampled flag (Sampled)
  # M   - # of imputations
  # ROB - indicator, if FALSE use model-based covariance
  # MARG.EXP.FORMULA - marginal exposure model formula: [Xe | Xo]
  # EXPVAR - exposure / expensive variable: Xe
  # Sampled - should be set to 1 if included in the sample and 0 otherwise
  #EXPVAR   = DAT$EXPVAR  = DAT[, as.character(substitute(EXPVAR))]
  Sampled  = DAT$Sampled = DAT[, as.character(substitute(Sampled))]
  
  imp_data <- vector('list', M)
  
  for(ixx in seq(M)) { if (verbose) cat("\r","Imputing dataset number:", ixx)
    # Step 0: Update date files
    # Create versions of DAT with updated Xe values
    DAT_Xe1 <- DAT_Xe0 <- DAT
    #DAT_Xe1$EXPVAR <- ifelse (!is.na(DAT_Xe1$EXPVAR), 1, 1 )
    #DAT_Xe0$EXPVAR <- ifelse (!is.na(DAT_Xe0$EXPVAR), 0, 0 )
    
    #print(summary(DAT_Xe1$EXPVAR))
    #print(summary(DAT_Xe0))
    DAT_Xe1[,EXPVAR][is.na(DAT_Xe1[,EXPVAR])] <- 1
    DAT_Xe0[,EXPVAR][is.na(DAT_Xe0[,EXPVAR])] <- 0
    DAT_Xe1[,EXPVAR] <- 1
    DAT_Xe0[,EXPVAR] <- 0
    #print(summary(DAT_Xe0))
    #DAT_Xe1[,EXPVAR & !is.na(EXPVAR)] <- 1
    #DAT_Xe0[,EXPVAR & !is.na(EXPVAR)] <- 0
    
    # Non-sampled subjects
    DAT_S0     <- DAT[Sampled==0,]
    DAT_S0_Xe1 <- DAT_Xe1[Sampled==0, ]
    DAT_S0_Xe1$Xe.i <- DAT_Xe1[Sampled==0, 'time']
    DAT_S0_Xe0 <- DAT_Xe0[Sampled==0, ]
    DAT_S0_Xe0$Xe.i <- 0
    dup        <- duplicated(DAT_S0$id)
    DAT_S0_v0  <- DAT_S0[!dup, ]
    
    # Sampled subjects
    DAT_S1     <- DAT[Sampled==1,]
    DAT_S1_Xe1 <- DAT_Xe1[Sampled==1, ]
    DAT_S1_Xe0 <- DAT_Xe0[Sampled==1, ]
    dup        <- duplicated(DAT_S1$id)
    DAT_S1_v0  <- DAT_S1[!dup,]
    
    #print(summary(DAT_S0_Xe1))
    #print(summary(DAT_S1_Xe1))
    # Step 1: extract estimates from FIT and draw m^th theta from MVN
    
    theta      <- c(FIT$beta, FIT$alpha)
    if(ROB) { cov_theta <- FIT$rob.cov
    } else { cov_theta <- FIT$mod.cov}
    
    theta.m <- mvrnorm(n = 1, mu=theta, Sigma=cov_theta)
    
    # Update theta.m to account for PROFILEVAR
    theta.m2 <- c( coef(attr(FIT,'FIT'))$beta, coef(attr(FIT,'FIT'))$alpha)
    theta.m2[match(names(theta.m),names(theta.m2))] <- theta.m
    theta.m2[match( attr(FIT, 'PROFILEVAR') ,names(theta.m2))] <- as.numeric(attr(FIT, 'BETAM')[match(attr(FIT, 'PROFILEVAR') ,names(theta.m2))]) 
    
    theta.m   <- theta.m2
    nms.inits <- names(theta.m)
    ix.mean   <- seq(length(nms.inits))
    ix.t      <- grep('gamma',nms.inits)
    ix.lv     <- grep('sigma',nms.inits)
    ix.mean   <- ix.mean[-c(which(ix.mean %in% c(ix.t, ix.lv)))]
    
    if(length(ix.lv)==0) {
      inits.lv <- -5 
    } else {
      inits.lv <- theta.m[ix.lv]
    }
    
    if(length(ix.t)==0) {
      inits.t <- 0 
    } else {
      inits.t <- theta.m[ix.t]
    }
    theta.m0 <- c( theta.m[ix.mean], inits.t ,inits.lv )
    
    
    # Step 2: Using the estimated theta=(beta,alpha), calculate likelihood contributions for sampled subjects
    # for both Xe=1 and Xe=0.  
    #print(names(attr(FIT,'args')))
    prep_input <- attr( attr(FIT,'FIT'), 'args' )
    
    tmp_S1_Xe1 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, 
                            t.formula=prep_input$t.formula, 
                            id=DAT_S1_Xe1[,prep_input$id], data=DAT_S1_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi, 
                            offset=prep_input$offset)$subjectData
    tmp_S1_Xe0 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula, 
                            id=DAT_S1_Xe0[,prep_input$id], data=DAT_S1_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs, samp.probi=prep_input$samp.probi, 
                            offset=prep_input$offset)$subjectData
    
    # These functions calculate the conditional likelihood contributions by and the ascertainment correction
    # for sampled subjects.  The difference of ascertainment corrections (i.e. log of the ratio 
    # log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo))= log(pr(S=1 | Xe=1, Xo))-log(pr(S=1 | Xe=0, Xo))) provides
    # the offset for the marginal exposure model 
    
    LLSC_1 <- MMLB:::LogLScoreCalc(params = theta.m0, subjectData = tmp_S1_Xe1, 
                                   Q = attr(FIT,'FIT')$LLSC_args$Q, W = attr(FIT,'FIT')$LLSC_args$W,
                                   Z = attr(FIT,'FIT')$LLSC_args$Z, 
                                   ParamLengths = attr(FIT,'FIT')$LLSC_args$ParamLengths,
                                   CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB) 
    LLSC_0 <- MMLB:::LogLScoreCalc(params = theta.m0, subjectData = tmp_S1_Xe0, 
                                   Q = attr(FIT,'FIT')$LLSC_args$Q, W = attr(FIT,'FIT')$LLSC_args$W, 
                                   Z = attr(FIT,'FIT')$LLSC_args$Z, 
                                   ParamLengths = attr(FIT,'FIT')$LLSC_args$ParamLengths,
                                   CondLike = FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    
    # Note: The ascertainment corrections used for the marginal exposure model offset 
    #       are already on the log scale, so log(a/b) = loga-logb
    logAC_1 <- attr(LLSC_1, "ACSubj") 
    logAC_0 <- attr(LLSC_0, "ACSubj")
    DAT_S1_v0$offset <- logAC_1-logAC_0
    
    # Perform offsetted logistic regression [Note, only using first visit info i.e. no time-varying information.
    fit.exp <- glm(MARG.EXP.FORMULA, offset=offset, data=DAT_S1_v0,family=binomial())
    
    # Extract alpha estimates and sample from the sampling distribution
    alpha      <- coef(fit.exp)
    cov_alpha  <- vcov(fit.exp)
    
    alpha.m <- mvrnorm(n = 1, mu=alpha, Sigma=cov_alpha)
    
    # Step 5: Apply the model and theta to calculate the ascertainment corrections for unsampled people, 
    
    nr <- nrow(DAT_S0_Xe1)
    
    tmp_S0_Xe1 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula, 
                            t.formula=prep_input$t.formula, 
                            id=DAT_S0_Xe1[,prep_input$id], data=DAT_S0_Xe1, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),], 
                            samp.probi=prep_input$samp.probi[rep(1,nr)], 
                            offset=prep_input$offset[rep(1,nr)])$subjectData
    tmp_S0_Xe0 <- prep_data(mean.formula=prep_input$mean.formula, lv.formula=prep_input$lv.formula,
                            t.formula=prep_input$t.formula, 
                            id=DAT_S0_Xe0[,prep_input$id], data=DAT_S0_Xe0, inits=theta.m,
                            samp.probs=prep_input$samp.probs[rep(1,nr),], 
                            samp.probi=prep_input$samp.probi[rep(1,nr)], 
                            offset=prep_input$offset[rep(1,nr)])$subjectData
    
    LLSC_1 <- MMLB:::LogLScoreCalc(params = theta.m0, subjectData = tmp_S0_Xe1, 
                                   Q = attr(FIT,'FIT')$LLSC_args$Q, 
                                   W = attr(FIT,'FIT')$LLSC_args$W, 
                                   Z = attr(FIT,'FIT')$LLSC_args$Z, 
                                   ParamLengths = attr(FIT,'FIT')$LLSC_args$ParamLengths,
                                   CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB) 
    LLSC_0 <- MMLB:::LogLScoreCalc(params = theta.m0, subjectData = tmp_S0_Xe0, 
                                   Q = attr(FIT,'FIT')$LLSC_args$Q, 
                                   W = attr(FIT,'FIT')$LLSC_args$W, 
                                   Z =  attr(FIT,'FIT')$LLSC_args$Z, 
                                   ParamLengths = attr(FIT,'FIT')$LLSC_args$ParamLengths,
                                   CondLike =  FIT$control['cond.like']==1, EmpiricalCheeseCalc = ROB)
    
    ## Applied to the non-sampled subjects: log(pr(S=1 | Xe=1, Xo)/pr(S=1 | Xe=1, Xo)).
    logAC_1 <- attr(LLSC_1, "ACSubj") # Ascertainment corrections already log transformed
    logAC_0 <- attr(LLSC_0, "ACSubj")
    offset_log_S0 <- logAC_1-logAC_0
    
    # Create a temporary dataframe for non-sampled subjects in order to create the
    # model.matrix that is used to calculate marginal exposure odds
    tmp          <- DAT_S0_v0
    tmp[,EXPVAR] <- sample(c(0,1),size=nrow(tmp),replace=TRUE)
    mm_tmp       <- model.matrix(as.formula(MARG.EXP.FORMULA), tmp) # Added as.formula
    
    # Compute the conditional exposure odds in non-sampled subjects
    # Notice that all conditional odds are equal: [Xe | Y, Xo, S=1] = [Xe | Y, Xo, S=0] = [Xe | Y, Xo] since we sampled based on Y
    # So, we multiply the odds (conditional on being sampled), pr[Xe=1 | Xo, S=1]/pr[Xe=0 | Xo, S=1] 
    # by the likelihood ratio (conditional on being sampled), pr[Y | Xe=1, Xo, S=1]/ pr[Y | Xe=1, Xo, S=1]
    # in the unsampled subjects in order to obtain conditional odds in unsampled: pr[Xe=1 | Y, Xo, S=1]/pr[Xe=0 | Y, Xo, S=1]
    marg.odds.Xe.Xo.S1   <- exp(mm_tmp %*% alpha.m)*exp(offset_log_S0)  
    LR.Y.Xe.Xo.S1        <- exp( attr(LLSC_1, "LogLikeSubj")-attr(LLSC_0, "LogLikeSubj"))
    cond.odds.Xe.Y.Xo.S0 <-  LR.Y.Xe.Xo.S1*marg.odds.Xe.Xo.S1  
    
    # Convert odds to probability, log(odds) = log(1/[1-p]) ---> p = odds/(1+odds)
    # Impute Xe for non-sampled subjects
    prXe1     <- Odds2Prob(cond.odds.Xe.Y.Xo.S0)
    XeS0      <- rbinom(nrow(tmp), size=1, prob=prXe1)
    
    # Add XeS0 info to DAT_S0_v0
    DAT_S0_v0[,EXPVAR] <- XeS0
    DAT_S0 <- DAT_S0[,-which(colnames(DAT_S0)==EXPVAR)]
    DAT_S0 <- merge(DAT_S0, DAT_S0_v0[,c('id',EXPVAR)],by='id')
    DAT_S0$Xe.i <- DAT_S0$Xe*DAT_S0$time
    
    DAT_imp <- rbind(DAT_S1, DAT_S0[,colnames(DAT_S1)])
    
    # Store imputation data
    imp_data[[ixx]] <- DAT_imp
    #if(verbose) cat('\r',ix)
  } # end of for ix loop
  
  imp_data
}
find_n2 <- function(MODEL, LVMODEL, TMODEL, Q, VNM, EVXEI, DAT_WXE, TMPSS, NREPS2=50, 
                    ITERLIM, VERBOSE,NSTART=10, QRPROB=.5,FIT1) {
  # VNM     = vnm argument (name of model parameter, 'Xe.i')
  # EVXEI   = target expected variance for Xe.i
  # DAT_WXE = complete data after feasibility analysis or mi
  
  # remove sp1 sp2 sp3 columns, if present
  blah <- colnames(DAT_WXE) %in% c('sp1','sp2','sp3')
  if(any(blah)) DAT_WXE <- DAT_WXE[,-which(blah)]
  
  # If we use D[0,X,0], we need to find N2.
  #fc_fit_wXe  <- mm(as.formula(MODEL),data=DAT_WXE[!is.na(DAT_WXE$samplek),], lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  fc_fit_wXe  <- mm(as.formula(MODEL),data=DAT_WXE, lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  #  EV          <- t(as.matrix(PROFILELC)) %*% vcov(fc_fit_wXe)$beta %*% as.matrix(PROFILELC)
  EV          <- vcov(fc_fit_wXe)$beta['Xe.i','Xe.i']
  N_FC        <- table(DAT_WXE[!duplicated(DAT_WXE$id),'ss'])[2]
  #  fit0        <- mm(as.formula(MODEL),data=DAT_WXE[!is.na(DAT_WXE$samplek),], lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  #  fit0        <- mm(as.formula(MODEL),data=DAT_WXE[DAT_WXE$samplek==1,], lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  fit0        <- FIT1 
  #  EV0         <- t(as.matrix(PROFILELC)) %*% vcov(fit0)$beta %*% as.matrix(PROFILELC)
  EV0          <- vcov(fit0)$beta['Xe.i','Xe.i']  
  #  N_0         <- table(factor( DAT_WXE[!duplicated(DAT_WXE$id) & !is.na(DAT_WXE$samplek),'ss'], levels=c(1,2,3)) )[2]
  N_0         <- table(factor( DAT_WXE[!duplicated(DAT_WXE$id) & DAT_WXE$samplek==1,'ss'], levels=c(1,2,3)) )[2]
  PLOTDAT     <- cbind(c(N_FC,N_0), log(c(EV,EV0)))
  colnames(PLOTDAT) <- c('N','logV')
  #NEWDATA     <- data.frame(N=seq(N_FC))
  
  # Sampling frequencies DAT_WXE (same as original data set since its based on Y)
  n_ss_n2   <- sapply(TMPSS, length)
  n1        <- length(unique(DAT_WXE$id))
  
  # Eligible data for stage 2 sampling
  #  DAT_WXE2   <- DAT_WXE[is.na(DAT_WXE$samplek), ] # data that can be sampled during stage 2: use this or S1D
  DAT_WXE2   <- DAT_WXE[DAT_WXE$samplek==0, ] # data that can be sampled during stage 2: use this or S1D
  DAT_WXE2_s <- split(DAT_WXE2, DAT_WXE2$id)
  tmp_ss_n2  <- idsXss(DAT_WXE2_s)
  n2         <- sapply(tmp_ss_n2,length)
  
  # idsXss tabulates row numbers (not ids), we will need to adjust these according to the following key
  id_key     <- data.frame('ix'=seq(length(DAT_WXE2_s)), 'id'=sapply(DAT_WXE2_s, function(ZZ)ZZ$id[1]), 
                           'ss'=sapply(DAT_WXE2_s, function(ZZ)ZZ$ss[1]))
  # See the following example
  # DAT_WXE2b <- DAT_WXE2
  # DAT_WXE2b$id <- DAT_WXE2b$id + 10000000
  # DAT_WXE2b_s <- split(DAT_WXE2b, DAT_WXE2b$id)
  # tmp_ss_n2b  <- idsXss(DAT_WXE2b_s) # notice that the ix are returned not the ids!
  
  if(EV>EVXEI) {
    # Store results
    outs <-  list('med_n2'=NULL,'mean_n2'=NULL, 'data'=NULL, 'threshold_met'=FALSE)
  } else {
    
    cont   <- TRUE
    ii     <- 0
    # NSTART <- 20# max(NSTART, 25)
    qs <- rep( seq(0.1, 0.9, 0.1),each=2)
    NSTART  <- length(qs)      # estimate variance for each qs value twice
    nreps2  <- NREPS2 + NSTART # number of reps AFTER initial NSTART 
    PLOTDAT <- rbind(PLOTDAT, matrix(NA, ncol=2, nrow=nreps2))
    tmp_out <- matrix(NA, ncol=2, nrow=nreps2)
    #PLOTDAT2 <- data.frame(matrix(NA, ncol=2, nrow=NREPS2))
    #colnames(PLOTDAT2) <- c('N','logV')
    
    while(cont & ii<nreps2) {
      ii      <- ii+1
      tmpdat  <- DAT_WXE 
      if(ii <= NSTART) {
        n2q <- ceiling(quantile( seq(n2[2]) , prob=qs[ii]))
      } else {
        
        # Find estimated log variance of Xei that is closest to EVEXI
        # plot(PLOTDAT[1:ii,1], PLOTDAT[1:ii,2])
        # abline(h=log(EVEXI))
        
        PLOTDAT <- PLOTDAT[order(PLOTDAT[,1]),]
        signs   <- ifelse(ii%%2==0,1,-1)
        
        find_sign_change <- max( which( !duplicated(sign(PLOTDAT[1:ii,2] - log(EVXEI+signs*EVXEI/5))) ) )
        n2_seq <- seq(PLOTDAT[find_sign_change-1,'N'], PLOTDAT[find_sign_change,'N'])
        n2q <- sample(n2_seq,1)
      }
      
      # if(ii<NSTART+1) {
      #   if(ii>1) PLOTDAT <- rbind(PLOTDAT, PLOTDAT2[ii-1,])
      #   fit_logV_N <- lm(logV~N,data=PLOTDAT)
      #   n2q <- ceiling( which( exp(predict(fit_logV_N, newdata=NEWDATA) )<=EVXEI )[1])
      # } else {
      #   if(ii==NSTART+1) PLOTDAT <- PLOTDAT[-seq(10),]
      #   PLOTDAT <- rbind(PLOTDAT, PLOTDAT2[ii-1,])
      #   #n2q <- ceiling( which( exp(predict(rq(logV~N,tau=QRPROB,data=PLOTDAT),newdata=NEWDATA) )<=EVXEI )[1])
      #   fit_logV_N <- lm(logV~N,data=PLOTDAT)
      #   n2q <- ceiling( which( exp(predict(fit_logV_N, newdata=NEWDATA) )<=EVXEI )[1])
      #   
      #   # kk       <- which(!is.na(tmp_out[,2]))
      #   # kk_dat   <- data.frame(cbind( yy=tmp_out[kk,2],xx=tmp_out[kk,1]) )
      #   # kk_dat2  <- data.frame(xx=seq(n2[2]))
      #   # 
      #   # preds_nl <- predict( with(kk_dat, loess(yy~xx,control=loess.control(surface="direct"))), newdata=kk_dat2 )
      #   # n2q      <- which(preds_nl <=EVXEI)[1]
      #   # if(is.na(n2q)) n2q <- max(kk_dat2)
      #   # if(n2q==max(kk_dat2) & ii>NSTART+2) cont<-FALSE
      # }
      
      # Perform extreme sampling design
      tmp_design <- c(0,n2q,0)
      
      # Create stage 2 sampling info and perform analysis
      s2d      <- apply( cbind( tmp_design, attr(tmp_ss_n2,'freq')$ss ), 1, min )
      s2info   <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_n2, RSEXACT=FALSE)
      
      if(length(s2info$idsXss_samp[[2]])==0) {
        # Force to sample at least 1
        pick_one <- sample(s2info$idsXss_remain[[2]],1, replace=FALSE)
        s2info$idsXss_remain[[2]] <- s2info$idsXss_remain[[2]][s2info$idsXss_remain[[2]]!=pick_one]
        s2info$idsXss_samp[[2]]   <- pick_one
      }
      
      # Update DAT_WXE2 with stage 2 information
      #s2id  <- id_key[sort( unique(unlist( s2info$idsXss_samp )) ),'id']
      s2id  <- id_key[match( sort( unique(unlist( s2info$idsXss_samp )) ), id_key[,'id'] ),'id']
      s2idX <- which(tmpdat$id%in%s2id)
      tmpdat$samplek[s2idX] <- 2
      tmpdat[s2idX,c('sprob1','sprob2','sprob3')] <- matrix( s2info$sp, nrow=1 )[rep(1, length(s2idX)),]
      
      # Estimate Stage K sampling probabilities
      prob_sampled     <- do.call(rbind, lapply( split(tmpdat, tmpdat$samplek), function(ZZ) ZZ[1,c('sprob1','sprob2','sprob3')] ))
      prob_sampled     <- prob_sampled[-1,] # omit 0 row in table
      prob_not_sampled <- apply( 1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
      prob_not_sampled <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
      prob_sampled     <- prob_sampled * prob_not_sampled 
      prob_sampled$samplek <- as.numeric( rownames(prob_sampled) )
      colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3')
      
      tmpdat2 <- merge(tmpdat, prob_sampled, by='samplek',all.x=TRUE)
      tmpdat3 <- tmpdat2[tmpdat2$samplek>0, ]
      tmpdat3 <- tmpdat3[ order(tmpdat3$id, tmpdat3$time), ]
      #tmpdat3 <- tmpdat2[!is.na(tmpdat2$samplek),]
      
      fit2 <- mm(as.formula(MODEL), dat=tmpdat3, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                 cond.like=TRUE, samp.probs=as.matrix( tmpdat3[,c('sp1','sp2','sp3')] ),
                 iter.lim=ITERLIM,q=Q)
      
      #      tmp_out[ii, 1:2] <- c(n2q, t(as.matrix(PROFILELC)) %*% vcov(fit2)$beta %*% as.matrix(PROFILELC))
      #      PLOTDAT2[ii,1:2] <- c(n2q,log(t(as.matrix(PROFILELC)) %*% vcov(fit2)$beta %*% as.matrix(PROFILELC)))
      tmp_out[ii, 1:2] <- c(n2q, vcov(fit2)$beta['Xe.i','Xe.i'] )
      #      PLOTDAT2[ii,1:2] <- c(n2q,log( vcov(fit2)$beta['Xe.i','Xe.i']))
      PLOTDAT[ii+2,1:2] <- c(n2q,log( vcov(fit2)$beta['Xe.i','Xe.i']))
      if(VERBOSE) cat('\n\n',tmp_out[ii,])
      #if(VERBOSE) cat('\r',ii)
    }
    
    n_dat  <- data.frame( tmp_out )
    colnames(n_dat) <- c('N','V')
    
    if(any(is.na(n_dat))) n_dat <- n_dat[!is.na(n_dat[,1]), ]
    n_dat$keep <- (n_dat$N<=quantile(n_dat$N,.90) & n_dat$N>quantile(n_dat$N,.1))*1
    n_dat2 <- n_dat[n_dat$keep==1, ]
    
    n_dat3 <- data.frame(N=seq(min(n_dat2[,1]),max(n_dat2[,1])))
    V      <- predict( loess(V~N,data=n_dat2,control=loess.control(surface="direct")), newdata=n_dat3 )
    n_dat3 <- cbind(n_dat3, V)
    
    # Store results
    meetEVXEI <- TRUE
    n2    <- n_dat3[ which(V <=EVXEI)[1], 1]
    n2var <- n_dat3[ which(V <=EVXEI)[1], 2]
    if(is.na(n2)) { 
      meetEVXEI <- FALSE
      n2 <- n2var <- NA
    }
    
    outs <- list('med_n2'=ceiling(median(tmp_out[-seq(NSTART),1])),
                 'mean_n2'=ceiling(mean(tmp_out[-seq(NSTART),1])), 
                 'data'=n_dat, 'threshold_met'=any(n_dat[,2]<=EVXEI),
                 'n2_lowess'=n2,'n2_lowess_var'=n2var)
  }
  outs
}

adaptive_ss <- function(DATA=NULL,SEEDS=1, YSEED=NULL, NLOW=5, NHIGH=5, N=5000, PXE=.25,
                        BETAZ=data.frame('Z'=1), BETAM=c(-1.50,-0.25,1,.25), GAMMA=2, SIGMA=NULL, 
                        LVMODEL=NULL, TMODEL=~1, MODEL='Y~time+Xe+Xe.i+Z',
                        PMODEL='Y~time', Q=2, ITERLIM=500, PROFILEVAR=c('Xe','Xe.i'),
                        VERBOSE=TRUE, S1D =c(25,50,25), BREP=250, PROFILEEV=0.01,MARG.EXP.FORMULA) {
  
  tA <- Sys.time()
  
  if(is.null(DATA)) { 
    if(VERBOSE) cat('\n\r',date(),'\n\r','Generating data')
    
    # Generate data
    set.seed(SEEDS)
    
    XI     <- gen_XI(SEED=35, NLOW, NHIGH, N, PXE, BETAZ)
    if(is.null(YSEED)) YSEED <- SEEDS
    sdat   <- gen_data(SEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, BETAZ=BETAZ,BETAM=BETAM, GAMMA=GAMMA, 
                       SIGMA=SIGMA, Q=Q, XI=XI,MEANFORMULA=MODEL, LVFORMULA=LVMODEL, TFORMULA=TMODEL)
    sdat   <- split(sdat, sdat$id)
    fcdat  <- data.frame(do.call(rbind, lapply(sdat, as.matrix)))
    
    fc_ss_dist <- table(sapply(sdat, function(ZZ) ZZ$ss[1]))
    
    adat <- fcdat
  } else {
    adat <- fcdat <- DATA
    sdat <- split(adat, adat$id)
  }
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing full cohort analysis')
  # Perform full cohort analysis
  fc_fit   <- mm(as.formula(MODEL),data=fcdat, lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing stage 1 analysis')
  # Perform stage 1 analysis
  tmp_ss <- idsXss(DAT=sdat)
  
  # Frequencies by sampling stratum
  n_ss   <- sapply(tmp_ss, length)
  
  # Create analysis data set for final analysis (fcdata augmented with sampling probabilities)
  adat <- fcdat
  adat$Xe.o    <- adat$Xe                               # Copy of true (original) Xe value
  adat$Xe.i.o  <- adat$Xe.i                             # Copy of original Xe.i value
  adat$Xe      <- adat$Xe.i <- NA                       # Replace Xe, Xe.i with NA 
  adat$samplek <- 0                                     # numeric value noting stage k value when sampled
  adat$sprob3  <- adat$sprob2 <- adat$sprob1 <- NA      # sampling probabilities for strata 1-3
  
  # Create stage 1 sampling info and perform analysis
  s1d      <- apply( cbind( S1D, attr(tmp_ss,'freq')$ss ), 1, min )
  s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss, RSEXACT=FALSE)
  
  # Update adat with stage 1 information
  s1id  <- sort( unique(unlist( s1info$idsXss_samp )) )
  s1idX <- which(adat$id%in%s1id)
  adat$samplek[s1idX] <- 1
  adat[s1idX,c('sprob1','sprob2','sprob3')] <- matrix( s1info$sp, nrow=1 )[rep(1, length(s1idX)),]
  adat$Xe[s1idX]   <- adat$Xe.o[s1idX]
  adat$Xe.i[s1idX] <- adat$Xe.i.o[s1idX]
  
  # Perform Stage 1 analysis
  #adat1  <- adat[!is.na(adat$samplek), ]
  adat1  <- adat[adat$samplek>0, ]
  fit1   <- mm(as.formula(MODEL), dat=adat1, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
               cond.like=TRUE, samp.probs=s1info$sp, iter.lim=ITERLIM, q=Q,return_args=TRUE)
  
  fit1_res <- list('coef'=coef(fit1),'vcov'=vcov(fit1))
  
  # Fix PROFILEVAR (treat as offset using profile likelihood - compare MODEL and PMODEL)
  dat1  <- do.call(rbind, sdat[as.numeric(unlist(s1info$idsXss_samp))])
  
  mf    <- model.matrix(as.formula(MODEL), data=dat1)
  pmf   <- mf
  pmf[, !( colnames(mf) %in% PROFILEVAR ) ] <- 0
  
  # pmf   <- model.matrix(as.formula(PMODEL), data=dat1)
  
  # betam <- coef NATEHERE #c(BETAM,BETAZ)
  # names(betam) <- colnames(mf)
  
  betam <- coef(fc_fit)$beta     
  names(betam) <- colnames(mf)
  
  #dat1$offset <- as.matrix( mf[,offset_vals] ) %*% as.matrix( betam[offset_vals] )
  dat1$offset <- as.matrix( pmf ) %*% as.matrix( betam )
  dat2 <- dat1
  dat2[,colnames(dat2) %in% PROFILEVAR] <- 0
  
  fit1b   <- mm(as.formula(MODEL), dat=dat2, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                cond.like=TRUE, samp.probs=s1info$sp, offset=dat2$offset, iter.lim=ITERLIM, q=Q,return_args = TRUE)
  
  attr(fit1b, 'PROFILEVAR') <- PROFILEVAR
  attr(fit1b, 'BETAM')      <- coef(fc_fit)$beta # BETAM
  attr(fit1b, 'FIT')        <- fit1
  
  # # Perform feasibility analysis, use FA data to estimate Xe using data from previous stages
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing interim design analysis')
  # # Create temporary data set
  adat2      <- adat

  tmp_ss_k      <- s1info$idsXss_remain
  tmp_ss_kn     <- sapply(tmp_ss_k,length)
  prob_sampled  <- matrix( s1info$sp,nrow=1 )
  colnames(prob_sampled) <- c('sp1','sp2','sp3')
  
  tmpXe <- ImputeData(FIT=fit1b, DAT=adat, M=BREP, MARG.EXP.FORMULA=MARG.EXP.FORMULA, Sampled='samplek', 
                      verbose=VERBOSE, EXPVAR='Xe')
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Estimating distribution of n2 of D[0,n2,0]')
  n2_fa    <- rep(NA, BREP)
  n2_fa_hx <- vector('list',BREP)
  n2_thres <- rep(NA, BREP)
  
  for(rx in seq(BREP)) {
    # FA version
    # adat2_tmp      <- adat2
    # tmpXe2         <- tmpXe
    # tmpXe2$Xe_em   <- rbinom(nrow(tmpXe), 1, prob=tmpXe2$pXe) 
    # adat2_tmp      <- merge(adat2_tmp, tmpXe2, by='id',all.x=TRUE)
    # adat2_tmp$Xe   <- adat2_tmp$Xe_em                
    # adat2_tmp$Xe.i <- adat2_tmp$Xe_em*adat2_tmp$time
    
    adat2_tmp <- tmpXe[[rx]]
    
    N2est <- find_n2(MODEL, LVMODEL, TMODEL, Q, VNM=PROFILEVAR, EVXEI=PROFILEEV, DAT_WXE=adat2_tmp, TMPSS=tmp_ss,
                  NREPS2=100, ITERLIM=ITERLIM, VERBOSE=VERBOSE, NSTART=30, QRPROB=.5,
                  FIT1=fit1)
    
    n2_fa[rx]      <- N2est$n2_lowess
    n2_fa_hx[[rx]] <- N2est$data
    n2_thres[rx]   <- N2est$threshold_met*1
    
    if(FALSE){
      pdf(file.path('~','Dropbox','Nate','Research','dissertation','paper2','draft','est_n2_example.pdf'),height=8,width=11)
      par(las=1)
      plot_dat <- N2est$data
      plot_dat <- plot_dat[-seq(NSTART),]
      plthres <- mean(plot_dat[,2]<=PROFILEEV)
      
      plot(plot_dat, typ='n',ylab=expression(paste("Var(", beta[et], ")")), 
           xlab=expression(paste(N[0*','*n[i]]^s,' Sample Size')),
           sub=substitute(paste('P(',widehat(V)<=thres,')=',p),list(thres=PROFILEEV,p=round(plthres,3))),
           cex.axis=.75)
      abline(h=PROFILEEV,col='lightgrey')
      text(plot_dat, labels=seq(nrow(N2est$data))[-seq(NSTART)], cex=0.5)
      aa   <- density(plot_dat[,1])
      aa$y <- aa$y/(max(aa$y)*(1/((par('usr')[3]+(par('usr')[4]-par('usr')[3])/10) -par('usr')[3])))
      aa$y <- aa$y+par('usr')[3]
      lines(aa$x,aa$y)
      #norm_dat <- dnorm(x=seq(min(N2$data[,1]),max(N2$data[,1])),mean=mean(N2$data[,1]), sd(N2$data[,1]))
      #norm_dat <- norm_dat/(max(norm_dat)*(1/((par('usr')[3]+(par('usr')[4]-par('usr')[3])/10) -par('usr')[3])))
      #  norm_dat <- norm_dat+par('usr')[3]
      #  lines(seq(min(N2$data[,1]),max(N2$data[,1])), norm_dat, col='red')
      dev.off()
    }
    if(VERBOSE) cat('\r',paste('Estimating distribution of n2 of D[0,n2,0]: ', rx, ' of ', BREP,sep=''))
  } # end of rx loop
  
  # Pick best design and perform analysis (and save n2_fa)
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing stage 2 analysis')
  
  # Empirical distribution of n2
  emp_n2 <- data.frame( table(n2_fa),as.numeric(prop.table(table(n2_fa)) ))
  colnames(emp_n2) <- c('n2','freq','prop')
  best_n2 <- ceiling(mean(n2_fa,na.rm=TRUE))
  
  # n2_dat    <- do.call(rbind, n2_fa_hx )
  # huh       <- exp(predict( lm(log(V)~N,data=n2_dat), newdata=data.frame('N'=seq(1,tmp_ss_kn[2]))))
  # n2_via_hx <- which(huh<=PROFILEEV)[1]
  
  d2 <- list( c(0,best_n2,0) )
  s2d    <- lapply( d2, function(ZZ) apply(cbind( ZZ, sapply(tmp_ss_k,length)), 1, min ))
  sKinfo <- lapply( s2d, function(ZZ) get_sample(DESIGN=ZZ, IDSXSS=tmp_ss_k, RSEXACT=FALSE) )
  
  #fit2_fa <- NULL
  
  for(dx in seq(length(d2))) {
    if(dx==1) {
      #dat_s  <- split( adat2, adat2$id )
      dat_ix <- which(adat2$id %in% unlist(sKinfo[[1]]$idsXss_samp))
      adat2$samplek[dat_ix] <- 2 
      adat2[dat_ix,c('sprob1','sprob2','sprob3')] <- matrix( sKinfo[[1]]$sp, nrow=1 )[rep(1, length(dat_ix)),]
      
      # Estimate Stage K sampling probabilities
      prob_sampled     <- do.call(rbind, lapply( split(adat2, adat2$samplek), function(ZZ) ZZ[1,c('sprob1','sprob2','sprob3')] ))
      prob_sampled     <- prob_sampled[-1,] # remove first row for 0
      prob_not_sampled <- apply( 1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
      prob_not_sampled <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
      prob_sampled     <- prob_sampled * prob_not_sampled 
      prob_sampled$samplek <- as.numeric( rownames(prob_sampled) )
      colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3') 
      
      adat3 <- merge(adat2, prob_sampled, by='samplek',all.x=TRUE)
      adat3 <- adat3[ order(adat3$id, adat3$time), ]
      adat3 <- adat3[adat3$samplek>0,]
      adat3$Xe <- adat3$Xe.o
      adat3$Xe.i <- adat3$Xe.i.o
      
      fit2 <- mm(as.formula(MODEL), dat=adat3, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                 cond.like=TRUE, samp.probs=as.matrix( adat3[,c('sp1','sp2','sp3')] ),
                 iter.lim=ITERLIM,q=Q)
    } 
  } # end of dx for loop
  
  tB <- Sys.time()
  outs <- list('seed'=SEEDS,'time'=difftime(tB,tA,units='mins'),
               'fc'= list('coef'=coef(fc_fit), 'vcov'=vcov(fc_fit), 'modcov'=fc_fit$mod.cov),
               'designs'=rbind(s1d, s2d[[1]]),'n2_dat'=n2_fa_hx,'n2_dist'=emp_n2,'n2_thres'=mean(n2_thres),
               'fit2'=list('coef'=coef(fit2), 'vcov'=vcov(fit2), 'modcov'=fit2$mod.cov), 
               'target'=list(PROFILEVAR,PROFILEEV),'fit1'=list('coef'=coef(fit1), 'vcov'=vcov(fit1)))
  outs
}

adaptive_des <- function(DATA=NULL,SEEDS=1, YSEED=NULL, NLOW=5, NHIGH=5, N=5000, PXE=.25,
                         BETAZ=data.frame('Z'=1), BETAM=data.frame('Int'=-1.50,'time'=-0.25,'Xe'=1,'Xe.i'=.25), 
                         GAMMA=2, SIGMA=NULL, 
                         LVMODEL=NULL, TMODEL=~1, MODEL='Y~time+Xe+Xe.i+Z',
                         Q=2, ITERLIM=500, PROFILEVAR=c('Xe','Xe.i'),VERBOSE=TRUE,
                         S1D =c(25,50,25), BREP=250,N2=400, MOD=10,MARG.EXP.FORMULA) {
  
  tA <- Sys.time()
  
  if(is.null(DATA)) { 
    if(VERBOSE) cat('\n\r',date(),'\n\r','Generating data')
    
    # Generate data
    set.seed(SEEDS)
    
    XI     <- gen_XI(SEED=35, NLOW, NHIGH, N, PXE, BETAZ)
    if(is.null(YSEED)) YSEED <- SEEDS
    sdat   <- gen_data(SEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, BETAZ=BETAZ,BETAM=BETAM, GAMMA=GAMMA, 
                       SIGMA=SIGMA, Q=Q, XI=XI,MEANFORMULA=MODEL, LVFORMULA=LVMODEL, TFORMULA=TMODEL)
    sdat   <- split(sdat, sdat$id)
    fcdat  <- data.frame(do.call(rbind, lapply(sdat, as.matrix)))
    
    fc_ss_dist <- table(sapply(sdat, function(ZZ) ZZ$ss[1]))
    
    adat <- fcdat
  } else {
    adat <- fcdat <- DATA
    sdat <- split(adat, adat$id)
  }
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing full cohort analysis')
  # Perform full cohort analysis
  fc_fit   <- mm(as.formula(MODEL),data=adat, lv.formula=LVMODEL,t.formula=TMODEL, id=id,q=Q)
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing stage 1 analysis')
  # Perform stage 1 analysis
  tmp_ss <- idsXss(DAT=sdat)
  
  # Frequencies by sampling stratum
  n_ss   <- sapply(tmp_ss, length)
  
  # Random sample
  adat      <- adat[order(adat$id, adat$time), ]
  adat_rs   <- adat
  ids       <- unique(adat$id)
  samp      <- sum(S1D)+N2
  ids_rs    <- sample(ids, samp, replace=FALSE)
  rs_sample <- (seq(nrow(adat)) %in% which(adat$id %in% ids_rs))*1
  adat_rs$sample <- rs_sample
  
  rs.tlv <- mm(as.formula(MODEL), dat=adat_rs[adat_rs$sample==1,], id=id,
               lv.formula=LVMODEL,  t.formula=TMODEL, verbose=VERBOSE, iter.lim=ITERLIM) 
  rs_ml  <- list('coef'=coef(rs.tlv),'vcov'=vcov(rs.tlv))
  
  # Create analysis data set for final analysis (fcdata augmented with sampling probabilities)
  adat$Xe.o    <- adat$Xe                               # Copy of true (original) Xe value
  adat$Xe.i.o  <- adat$Xe.i                             # Copy of original Xe.i value
  adat$Xe      <- adat$Xe.i <- NA                       # Replace Xe, Xe.i with NA 
  adat$samplek <- 0                                     # numeric value noting stage k value when sampled
  adat$sprob3  <- adat$sprob2 <- adat$sprob1 <- NA      # sampling probabilities for strata 1-3
  
  # Create stage 1 sampling info and perform analysis
  s1d      <- apply( cbind( S1D, attr(tmp_ss,'freq')$ss ), 1, min )
  s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss, RSEXACT=FALSE)
  
  # Update adat with stage 1 information
  s1id  <- sort( unique(unlist( s1info$idsXss_samp )) )
  s1idX <- which(adat$id%in%s1id)
  adat$samplek[s1idX] <- 1
  adat[s1idX,c('sprob1','sprob2','sprob3')] <- matrix( s1info$sp, nrow=1 )[rep(1, length(s1idX)),]
  adat$Xe[s1idX]   <- adat$Xe.o[s1idX]
  adat$Xe.i[s1idX] <- adat$Xe.i.o[s1idX]
  
  # Perform Stage 1 analysis
  adat1  <- adat[(adat$samplek==1), ] #is.na(adat$samplek) | 
  fit1   <- mm(as.formula(MODEL), dat=adat1, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
               cond.like=TRUE, samp.probs=s1info$sp, iter.lim=ITERLIM, q=Q,return_args=TRUE)
  
  fit1_res <- list('coef'=coef(fit1),'vcov'=vcov(fit1))
  
  # Fix PROFILEVAR (treat as offset using profile likelihood - compare MODEL and PMODEL)
  dat1  <- do.call(rbind, sdat[as.numeric(unlist(s1info$idsXss_samp))])
  
  mf    <- model.matrix(as.formula(MODEL), data=dat1)
  pmf   <- mf
  pmf[, !( colnames(mf) %in% PROFILEVAR ) ] <- 0
  
  # pmf   <- model.matrix(as.formula(PMODEL), data=dat1)
  #betam <- c(BETAM,BETAZ)
  betam <- coef(fc_fit)$beta     
  names(betam) <- colnames(mf)
  # Update betam with profilevals
  # See arguments PVAL_XeXei, VAL_Xe and PVAL_Xei
  
  #dat1$offset <- as.matrix( mf[,offset_vals] ) %*% as.matrix( betam[offset_vals] )
  dat1$offset <- as.matrix( pmf ) %*% as.matrix( betam )
  dat2 <- dat1
  dat2[,colnames(dat2) %in% PROFILEVAR] <- 0
  
  fit1b   <- mm(as.formula(MODEL), dat=dat2, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                cond.like=TRUE, samp.probs=s1info$sp, offset=dat2$offset, iter.lim=ITERLIM, q=Q,return_args = TRUE)
  
  attr(fit1b, 'PROFILEVAR') <- PROFILEVAR
  attr(fit1b, 'BETAM')      <- coef(fc_fit)$beta # BETAM
  attr(fit1b, 'FIT')        <- fit1
  
  # # Perform IDE, use FA data to estimate Xe using data from previous stages
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing interim design analysis')
  # # Create temporary data set
  adat2      <- adat

  tmp_ss_k      <- s1info$idsXss_remain
  tmp_ss_kn     <- sapply(tmp_ss_k,length)
  prob_sampled  <- matrix( s1info$sp,nrow=1 )
  colnames(prob_sampled) <- c('sp1','sp2','sp3')
  
  tmpXe <- ImputeData(FIT=fit1b, DAT=adat, M=BREP, MARG.EXP.FORMULA=MARG.EXP.FORMULA, Sampled='samplek', verbose=VERBOSE, EXPVAR='Xe')
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','For fixed N2, calculating possible designs')
  poss_des2 <- cbind(seq(0, min(tmp_ss_kn[c(1,3)]),by=MOD), 0, seq(0, min(tmp_ss_kn[c(1,3)]),by=MOD))
  poss_des2[,2] <- N2-rowSums( poss_des2 ) 
  poss_des2     <- poss_des2[ !(poss_des2[,2] < poss_des2[,1]), ]
  
  # Compute results for each possible design 
  poss_des2_list  <- split(poss_des2, row(poss_des2))
  nr_poss_des     <- nrow(poss_des2)
  true_results    <- vector('list', nr_poss_des)  # analyzing data as a two-stage fixed design
  true_results2   <- vector('list', nr_poss_des)  # analyzing data as a singe-stage fixed design (same data as two-stage)
  ss_w_des1_des2  <- vector('list', nr_poss_des)  # analyzing data as a singe-stage fixed design (diff data)
  
  for(nnx in seq(nr_poss_des)) {
    adat2_tmp      <- adat2
    adat2_tmp$Xe   <- adat2_tmp$Xe.o
    adat2_tmp$Xe.i <- adat2_tmp$Xe.i.o
    
    s2d    <- apply(cbind( poss_des2_list[[nnx]], sapply(tmp_ss_k,length)), 1, min )
    s2info <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_k, RSEXACT=FALSE)
    dat_s  <- split( adat2_tmp, adat2_tmp$id )
    dat_ix <- which(adat2_tmp$id %in% unlist(s2info$idsXss_samp))
    adat2_tmp$samplek[dat_ix] <- 2
    adat2_tmp[dat_ix,c('sprob1','sprob2','sprob3')] <- matrix( s2info$sp, nrow=1 )[rep(1, length(dat_ix)),]
    
    adat2_tmp <- adat2_tmp[adat2_tmp$samplek>0,]
    
    # Estimate Stage K sampling probabilities
    prob_sampled     <- do.call(rbind, lapply( split(adat2_tmp, adat2_tmp$samplek), function(ZZ) ZZ[1,c('sprob1','sprob2','sprob3')] ))
    prob_not_sampled <- apply( 1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
    prob_not_sampled <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])
    prob_sampled     <- prob_sampled * prob_not_sampled
    prob_sampled$samplek <- as.numeric( rownames(prob_sampled) )
    colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3')
    
    adat2_tmp <- merge(adat2_tmp, prob_sampled, by='samplek',all.x=TRUE)
    adat2_tmp <- adat2_tmp[ order(adat2_tmp$id, adat2_tmp$time), ]
    adat2_tmp2 <- adat2_tmp[!is.na(adat2_tmp$samplek),]
    
    fit2 <- mm(as.formula(MODEL), dat=adat2_tmp2, id=id, lv.formula=LVMODEL, t.formula=TMODEL,
               cond.like=TRUE, samp.probs=as.matrix( adat2_tmp2[,c('sp1','sp2','sp3')] ),
               iter.lim=ITERLIM,q=Q)
    true_results[[nnx]] <- list('coef'=coef(fit2),'vcov'=vcov(fit2))
    
    # Analyze data as if from a single stage
    adat_junk         <- fcdat
    adat_junk$Xe.o    <- adat_junk$Xe                               # Copy of true (original) Xe value
    adat_junk$Xe.i.o  <- adat_junk$Xe.i                             # Copy of original Xe.i value
    adat_junk$Xe      <- adat_junk$Xe.i <- NA                       # Replace Xe, Xe.i with NA 
    adat_junk$samplek <- 0                                     # numeric value noting stage k value when sampled
    adat_junk$sprob3  <- adat_junk$sprob2 <- adat_junk$sprob1 <- NA      # sampling probabilities for strata 1-3
    
    s12info <- get_sample(DESIGN=s1d + s2d, IDSXSS=tmp_ss, RSEXACT=FALSE)
    s12id   <- sort( unique(unlist( s12info$idsXss_samp )) )
    s12idX  <- which(adat_junk$id%in%s12id)
    adat_junk$samplek[s12idX] <- 1
    adat_junk[s12idX,c('sprob1','sprob2','sprob3')] <- matrix( s12info$sp, nrow=1 )[rep(1, length(s12idX)),]
    adat_junk$Xe[s12idX]   <- adat$Xe.o[s12idX]
    adat_junk$Xe.i[s12idX] <- adat$Xe.i.o[s12idX]
    
    fit12_junk <- mm(as.formula(MODEL), dat=adat_junk[!is.na(adat_junk$samplek), ], id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                     cond.like=TRUE, samp.probs=s12info$sp, iter.lim=ITERLIM, q=Q)
    ss_w_des1_des2[[nnx]] <- list('coef'=coef(fit12_junk),'vcov'=vcov(fit12_junk))
    
    # Very similar fit to the previous model, but this is using the actually observed data in the 2-stage experiment
    # not new data using different sampling probabilities.
    fit12 <- mm(as.formula(MODEL), dat=adat2_tmp2, id=id, lv.formula=LVMODEL, t.formula=TMODEL,
                cond.like=TRUE, samp.probs=s12info$sp,iter.lim=ITERLIM,q=Q)
    true_results2[[nnx]] <- list('coef'=coef(fit12),'vcov'=vcov(fit12))
  }
  attr(true_results,'poss_des') <- poss_des2
  
  # Estimate optimality criteria for all candidate designs 

  optimal_mat <- rep( list(matrix(NA, nrow=BREP, ncol=3)), nr_poss_des )
  
  for(bx in seq(BREP)) {
    for(nnx in seq(nr_poss_des)) {
      adat2_tmp2 <- tmpXe[[bx]]
      s2d        <- apply(cbind( poss_des2_list[[nnx]], sapply(tmp_ss_k,length)), 1, min )
      s2info     <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_k, RSEXACT=FALSE) 
      dat_s      <- split( adat2_tmp2, adat2_tmp2$id )
      dat_ix     <- which(adat2_tmp2$id %in% unlist(s2info$idsXss_samp))
      adat2_tmp2$samplek[dat_ix] <- 2 
      adat2_tmp2[dat_ix,c('sprob1','sprob2','sprob3')] <- matrix( s2info$sp, nrow=1 )[rep(1, length(dat_ix)),]
      
      adat2_tmp2 <- adat2_tmp2[adat2_tmp2$samplek>0,]
      
      # Estimate Stage K sampling probabilities
      prob_sampled     <- do.call(rbind, lapply( split(adat2_tmp2, adat2_tmp2$samplek), function(ZZ) ZZ[1,c('sprob1','sprob2','sprob3')] ))
      prob_not_sampled <- apply( 1-prob_sampled, 2, cumprod)                      # cumulative prob of not being sampled
      prob_not_sampled <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
      prob_sampled     <- prob_sampled * prob_not_sampled 
      prob_sampled$samplek <- as.numeric( rownames(prob_sampled) )
      colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3') 
      
      adat2_tmp2 <- merge(adat2_tmp2, prob_sampled, by='samplek',all.x=TRUE)
      adat2_tmp2 <- adat2_tmp2[ order(adat2_tmp2$id, adat2_tmp2$time), ]
      #adat2_tmp3 <- adat2_tmp2[!is.na(adat2_tmp2$samplek),]
      
      fit2 <- mm(as.formula(MODEL), dat=adat2_tmp2, id=id, lv.formula=LVMODEL, t.formula=TMODEL, 
                 cond.like=TRUE, samp.probs=as.matrix( adat2_tmp2[,c('sp1','sp2','sp3')] ),
                 iter.lim=ITERLIM,q=Q)
      
      optimal_mat[[nnx]][bx,] <- get_opt(DAT=fit2)
    }
  }
  
  optimal_mat <- lapply(optimal_mat, function(ZZ) {
    ZZ <- data.frame(ZZ) 
    colnames(ZZ) <- c('var_gt', 'var_g', 'd_opt')
    ZZ})
  
  tB <- Sys.time()
  
  outs <- list('seed'=SEEDS,'time'=difftime(tB,tA,units='mins'),
               'fc'= list('coef'=coef(fc_fit), 'vcov'=vcov(fc_fit), 'modcov'=fc_fit$mod.cov),
               'rs'=rs_ml,
               'optimal_mat'=optimal_mat, 'true_results'=true_results, 'fit1'=fit1_res,
               'true_results2'=true_results2,'ss_w_des1_des2'=ss_w_des1_des2)
  outs
}

fixed2stage <- function(DATA=NULL, SEEDS=1, YSEED=NULL, NLOW=5, NHIGH=5, N=5000, PXE=.25,
                        BETAZ=data.frame('Z'=1), BETAM=data.frame('Int'=-1.50,'time'=-0.25,'Xe'=1,'Xe.i'=.25), 
                        GAMMA=2, SIGMA=NULL, LVMODEL=NULL, TMODEL=~1, MODEL='Y~time+Xe+Xe.i+Z',
                        Q=2, ITERLIM=500,VERBOSE=TRUE, SAMP=500, SS_DES, TS_DES1, TS_DES2) {
  tA <- Sys.time()
  
  if(is.null(DATA)) { 
    if(VERBOSE) cat('\n\r',date(),'\n\r','Generating data')
    
    # Generate data
    set.seed(SEEDS)
    
    XI     <- gen_XI(SEED=35, NLOW, NHIGH, N, PXE, BETAZ)
    if(is.null(YSEED)) YSEED <- SEEDS
    sdat   <- gen_data(SEED=YSEED, NLOW=NLOW, NHIGH=NHIGH, N=N, PXE=PXE, BETAZ=BETAZ,BETAM=BETAM, GAMMA=GAMMA, 
                       SIGMA=SIGMA, Q=Q, XI=XI,MEANFORMULA=MODEL, LVFORMULA=LVMODEL, TFORMULA=TMODEL)
    sdat   <- split(sdat, sdat$id)
    fcdat  <- data.frame(do.call(rbind, lapply(sdat, as.matrix)))
    
    fc_ss_dist <- table(sapply(sdat, function(ZZ) ZZ$ss[1]))
    
    adat <- fcdat
  } else {
    adat <- DATA
  }
  
  if(VERBOSE) cat('\n\r',date(),'\n\r','Performing full cohort analysis')
  
  # Full cohort analysis
  fc.tlv <- mm(as.formula(MODEL), dat=adat, id=id, lv.formula=LVMODEL,  t.formula=TMODEL, verbose=FALSE, iter.lim=ITERLIM)
  fc_ml  <- cbind( unlist(coef(fc.tlv)),unlist( lapply( vcov(fc.tlv), function(ZZ) sqrt(diag(as.matrix(ZZ))) )) )
  initial.vals <- c(fc.tlv$beta,fc.tlv$alpha)
  
   
  if(VERBOSE) cat('\n\r',date(),'\n\r',paste('Performing RS of size ',SAMP,sep=''))
  # Random sample
  adat_rs   <- adat
  ids       <- unique(adat$id)
  ids_rs    <- sample(ids, SAMP, replace=FALSE)
  rs_sample <- (seq(nrow(adat)) %in% which(adat$id %in% ids_rs))*1
  adat_rs$sample <- rs_sample
  
  rs.tlv <- mm(as.formula(MODEL), dat=adat_rs[adat_rs$sample==1,], id=id, lv.formula=LVMODEL,  t.formula=TMODEL, verbose=FALSE, 
               iter.lim=ITERLIM,return_args=TRUE) # return_args for CD+iMI
  rs_ml  <- list(unlist(coef(rs.tlv)), rs.tlv$mod.cov)#cbind( unlist(coef(rs.tlv)),unlist( lapply( vcov(rs.tlv), function(ZZ) sqrt(diag(as.matrix(ZZ))) )) )
  
  if(VERBOSE) cat('\n\r',date(),'\n\r',paste('Performing single stage ODS designs'))
  # Single-stage ODS
  n_ss    <- length(SS_DES)
  
  ss_acml <- vector('list',n_ss)
  ss_samp_probs <- vector('list', n_ss)
  
  for(sx in seq(n_ss)) {
    adat_ss  <- adat
    tmp_ss   <- idsXss(DAT=split(adat_ss,adat_ss$id)) 
    s1d      <- apply( cbind( SS_DES[[sx]], attr(tmp_ss,'freq')$ss ), 1, min )
    s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss,RSEXACT=FALSE)
    s1id     <- sort( unique(unlist( s1info$idsXss_samp )) )
    adat_ss$sample <- (adat_ss$id%in%s1id)*1
    adat_ss  <- adat_ss[ order(adat_ss$id, adat_ss$time), ]
    ss_sp    <- s1info$sp
    
    ss_dat            <- adat_ss[adat_ss$sample==1,]
    ss_dat$samp.probi <- ss_sp[ss_dat$ss]
    ss.tlv <- mm(as.formula(MODEL), dat=ss_dat, id=id, lv.formula=LVMODEL,  t.formula=TMODEL, verbose=FALSE, 
                 iter.lim=ITERLIM,return_args=TRUE, cond.like=TRUE, samp.probs=as.numeric(ss_sp)) # return_args for CD+iMI
    ss_acml[[sx]]  <- list(unlist(coef(ss.tlv)), ss.tlv$mod.cov)
    
    ss_samp_probs[[sx]] <- ss_sp
  }
  
  if(VERBOSE) cat('\n\r',date(),'\n\r',paste('Performing two stage fixed ODS designs'))
  # Two-stage ODS
  tsd  <- expand.grid(seq(length(TS_DES1)),seq(length(TS_DES2)))
  n_ts <- nrow(tsd)
  
  ts_acml <- vector('list',n_ts)
  ts_samp_probs <- vector('list', n_ts)
  
  for(tx in seq(n_ts)) {
    adat_ts  <- adat
    
    tsx <- tsd[tx,]
    
    s1d      <- apply( cbind( TS_DES1[[tsx[1,1]]], attr(tmp_ss,'freq')$ss ), 1, min )
    s1info   <- get_sample(DESIGN=s1d, IDSXSS=tmp_ss, RSEXACT=FALSE)
    s1id     <- sort( unique(unlist( s1info$idsXss_samp )) )
    adat_ts$sample <- (adat_ts$id%in%s1id)*1
    
    # Stage 2
    tmp_ss_2   <- s1info$idsXss_remain
    tmp_ss_2n  <- sapply(tmp_ss_2,length)
    s2d        <- apply( cbind( TS_DES2[[tsx[1,2]]], sapply(tmp_ss_2,length)), 1, min )
    s2info     <- get_sample(DESIGN=s2d, IDSXSS=tmp_ss_2, RSEXACT=FALSE) 
    adat_ts$sample <- adat_ts$sample + (adat_ts$id %in% unlist(s2info$idsXss_samp))*2
    
    # Estimate phase 2 sampling probabilities
    prob_sampled        <- data.frame( rbind(s1info$sp,s2info$sp) )
    prob_not_sampled    <- apply( 1-prob_sampled, 2, cumprod)                  # cumulative prob of not being sampled
    prob_not_sampled    <- rbind(1, prob_not_sampled[-nrow(prob_not_sampled),])   
    prob_sampled        <- prob_sampled * prob_not_sampled 
    prob_sampled$sample <- c(1,2)
    colnames(prob_sampled)[1:3] <- c('sp1','sp2','sp3')
    
    adat_ts <- merge(adat_ts, prob_sampled, by='sample',all.x=TRUE)
    adat_ts <- adat_ts[ order(adat_ts$id, adat_ts$time), ]
    
    ts_dat <- adat_ts[adat_ts$sample>0,]
    ts_sp  <- ts_dat[,c('sp1','sp2','sp3')]
    ts_dat$samp.probi <- ts_sp[cbind(seq_along(nrow(ts_dat)),ts_dat$ss)]
    ts_dat <- ts_dat[order(ts_dat$id, ts_dat$time), ]
    
    ts.tlv <- mm(as.formula(MODEL), dat=ts_dat, id=id, lv.formula=LVMODEL,  t.formula=TMODEL, verbose=FALSE, 
                 iter.lim=ITERLIM,return_args=TRUE, cond.like=TRUE, samp.probs=as.matrix(ts_sp)) # return_args for CD+iMI
    ts_acml[[tx]]  <-list(unlist(coef(ts.tlv)), ts.tlv$mod.cov) 
    
    ts_samp_probs[[tx]] <- as.numeric( t(prob_sampled[,-4]) )
  }
  
  
  outs <- list(fc_ml, rs_ml,ss_acml,ts_acml)
  attr(outs, 'ss_sp') <- ss_samp_probs
  attr(outs, 'ts_sp') <- ts_samp_probs
  attr(outs, 'tsd')   <- tsd
  outs
}