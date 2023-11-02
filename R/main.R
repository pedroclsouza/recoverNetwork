#' Recover Network
#' 
#' Package from paper "Identifying Network Ties from Panel Data: Theory and an Application to Tax Competition"
#' @param data Panel data. See vignette for format and restriction.
#' @param lambda Penalization parameters in the format (p1,p1*,p2)
#' @param exoeffects Restricts gamma=0 if exoeffects=0 (optional, default is exoeffects=1)
#' @param timeweights Introduces time weights whcih is either 1 (default) or a vector with the same length as the number of time periods. This allows for reweighting the observations and estimate the version of the codes with the dynamic Wt.
#' @param eta Used in the computation of the initial conditions as explained in Step 4 in Appendix B.2. Default is 0.05.
#' @param Wder The derivative of the objective function with respect to the elements of Wij. This can be user-supplied to avoid the cost of computing the derivative. Computational gains for large N can be substantial. Default is NA.
#' @param adaestmethod Selection of estimation method for the adaptive step. Only supported version (default) is cvxr.
#' @param dopostols Equal to 1 if post-OLS is to be performed. Default is not using post-OLS, so dopostols=0, and final estimates are drawn from the unpenalized GMM.
#' @param Wfixed Used if part of the network is known. In this case, supply a NxN matrix with NA in the (i,j) elements that that are not known (and thus to be estimated), and with the known network intensity in the elements that are known. Supply Wfixed=-1 if network is entirely unknown and to be estimated.
#' @return List with initial conditions, elastic net, adaptive elastic net, and unpenalized GMM steps. See vignette for format.

recoverNetwork <- function(data,lambda,exoeffects=1,docv=0,timeweights=1,eta=0.05,Wder=NA,adaestmethod="cvxr",dopostols=0,Wfixed=-1){
  
  # settings ==================================================================
  
  param <- list()
  
  param$opt$outermeth='L-BFGS-B'
  param$opt$outeriter=100
  param$opt$en$estmethod='glmnet-coupled'
  param$opt$en$rsn=1
  param$opt$adaen$estmethod=adaestmethod
  param$opt$adaen$rsn=1
  param$r <- 0
  
  param$Wder=Wder
  param$Wfixed=Wfixed
  param$Wfixed0=Wfixed
  
  param$exoeffects=exoeffects
  param$standardise=1
  
  param$kappa=2.5
  param$eta=eta 
  param$verbose=1
  
  param$N <- length(unique(data$id))
  param$T <- length(unique(data$time))
  param$K <- length(grep("[x]",colnames(data)))
  
  param$lambdaL1 <- lambda[1]
  param$lambdaL2 <- lambda[2]
  param$lambdaL1Star <- lambda[3]
  
  param$fe <- c("id","time")
  
  param$timeweights <- timeweights
  
  # force main diagonal of Wfixed to be 0
  if (length(param$Wfixed)==1) {
    param$Wfixed <- matrix(NA,nrow=param$N,ncol=param$N)
  }
  diag(param$Wfixed) <- 0
  
  param$dopostpop <- 0
  param$dopostols <- dopostols
  param$doexternaliv <- 0
  
  # organise data ==============================================================
  
  datalist <- organiseData(data,param)
  
  data <- datalist$data
  datasq <- datalist$datasq
  
  result <- list()
  result$sdx <- datalist$sdx
  result$sdxnorm <- datalist$sdxnorm
  
  # initial conditions =========================================================
  
  if (param$verbose==1){message('Initial conditions')}
  
  # coefficients for the non-network model
  
  formula <- paste0("y~",paste(colnames(data)[grep("[x]",colnames(data))],collapse = "+"))
  betahat <- stats::lm(formula,data=data)
  betahat <- betahat$coefficients[-1]
  
  if (param$exoeffects==1) {
    p = c(0.5,betahat,rep(0,length(betahat)))
  } else{
    p = c(0.5,betahat)
  }
  
  # pre-filtering of elements of W
  # note: the function computeDerivatives is not computationally efficient
  
  if (length(param$Wder)==1) {
    param$Wder <- computeDerivatives(data,datasq,param)
    result$Wder <- param$Wder
  } else {
    message('skip computation of derivative')
  }
  
  WfixedWder<- (-param$Wder > param$eta)+0
  iX_Wfixed <- which(WfixedWder==0)
  param$Wfixed[iX_Wfixed] <- 0
  
  if (length(param$Wfixed0)>1) {
    iiX <- which(!is.na(param$Wfixed0))
    param$Wfixed[iiX] <- param$Wfixed0[iiX] 
  }
  
  result$ini$p <- p
  result$ini$Wfixed <- param$Wfixed
  
  # elastic net ================================================================
  
  if (param$verbose==1){message('Elastic Net')}
  
  # penalty
  
  penalty <- param$Wfixed
  penalty[is.na(penalty)] <- 1
  penalty[penalty==0] <- Inf
  result$en$penalty <- penalty
  
  # run optimisation
  
  o1 = outerOptim(p,lambdaL1=param$lambdaL1,lambdaL2=param$lambdaL2,penalty=penalty,estmethod=param$opt$en$estmethod,rsn=param$opt$en$rsn,datasq,param)
  W1 = f(as.matrix(o1[1:(2*param$K+1)]),lambdaL1=param$lambdaL1,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$en$estmethod,rsn=param$opt$en$rsn,datasq,param,returnW=1)
  p1 = o1[1:(2*param$K+1)]
  colnames(p1)<-NULL
  rownames(p1)<-NULL
  
  # save results
  
  result$en$optim$value <- o1$value
  result$en$W <- W1
  result$en$p <- p1
  
  result$en$rho = p1[1][[1]]
  result$en$beta = p1[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$en$gamma = p1[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$en$gamma = 0*result$en$beta 
  }
  
  # adaptive elastic net =======================================================
  
  if (param$verbose==1){message('Adaptive Elastic Net')}
  
  # penalty
  
  penalty <- result$en$W
  penalty[penalty <= 0.05] = 0.05
  penalty <- penalty^(-param$kappa)
  penalty[which(param$Wfixed==0)] <- Inf
  result$adaen$penalty <- penalty
  
  # run optimisation
  
  o2 = outerOptim(p,lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty=penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param)
  
  if (param$exoeffects==1){
    W2 = f(as.matrix(o2[1:(2*param$K+1)]),lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param,returnW=1)
    p2 = o2[1:(2*param$K+1)]
  } else {
    W2 = f(as.matrix(o2[1:(param$K+1)]),lambdaL1=param$lambdaL1Star,lambdaL2=param$lambdaL2,penalty,estmethod=param$opt$adaen$estmethod,rsn=param$opt$adaen$rsn,datasq,param,returnW=1)
    p2 = o2[1:(param$K+1)]
  }
  
  colnames(p2)<-NULL
  rownames(p2)<-NULL
  
  # save results
  
  result$adaen$optim$value <- o2$value
  result$adaen$W <- W2
  result$adaen$p <- p2
  
  result$adaen$rho = p2[1][[1]]
  result$adaen$beta = p2[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$adaen$gamma = p2[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$adaen$gamma = 0*result$adaen$beta
  }
  
  # unpenalised gmm ============================================================
  
  if (param$verbose==1){message('Unpenalised GMM')}
  
  WfixedSaved <- param$Wfixed
  param$Wfixed <- result$adaen$W
  
  # penalty
  penalty <- param$Wfixed
  penalty[is.na(penalty)] <- 1
  penalty[penalty==0] <- Inf
  
  # run optimisation
  
  o3 = outerOptim(unlist(p2),lambdaL1=0,lambdaL2=0,penalty=penalty,estmethod="givenw",rsn=param$opt$adaen$rsn,datasq,param)
  if (param$exoeffects==1){
    p3 = o3[1:(2*param$K+1)]
  } else {
    p3 = o3[1:(param$K+1)]
  }
  colnames(p3)<-NULL
  rownames(p3)<-NULL
  
  # save results
  result$unpenalisedgmm$optim$value <- o3$value
  result$unpenalisedgmm$p <- p3
  result$unpenalisedgmm$W <- result$adaen$W
  
  result$unpenalisedgmm$rho = p3[1][[1]]
  result$unpenalisedgmm$beta = p3[2:(param$K+1)][[1]] * result$sdxnorm
  if (param$exoeffects==1){
    result$unpenalisedgmm$gamma = p3[(param$K+2):(2*param$K+1)][[1]] * result$sdxnorm
  } else {
    result$unpenalisedgmm$gamma = 0*result$unpenalisedgmm$beta
  }
  
  # asymtotics (work in progress) ==============================================
  
  rho = result$unpenalisedgmm$rho
  beta = unlist(result$unpenalisedgmm$p[2:(param$K+1)])
  if (param$exoeffects==1){
    gamma = unlist(result$unpenalisedgmm$p[(param$K+2):(2*param$K+1)])
  } else {
    gamma = 0*beta
  }
  
  G_rho <- (momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="rho",1))
  G_beta <- sapply(1:param$K,function(i) momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="beta",i))
  G_gamma <- sapply(1:param$K,function(i) momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="gamma",i))
  iX <- which(is.na(WfixedSaved))
  G_W <- sapply(iX,function(i) (momentDer_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param,s="W",iX)))
  
  if (param$exoeffects==1){ 
    
    G <- (cbind(G_rho,G_beta,G_gamma, G_W))
    
  } else {
    
    G <- (cbind(G_rho,G_beta, G_W))
    
  }
  
  Omega <- moment_asyvar(rho=rho,beta=beta,gamma=gamma,W=result$unpenalisedgmm$W,datasq,param)
  
  invGtG <- ginv(t(G)%*%G)
  GtG <- t(G)%*%G
  vtheta <- ginv(GtG) %*% t(G)%*%Omega%*%G %*% ginv(GtG)
  vtheta <- vtheta/(param$T)
  
  vtheta <- diag(vtheta)
  
  result$unpenalisedgmm$asyvar <- vtheta
  
  # post peers-of-peers IV =====================================================
  
  if (param$dopostpop==1) {
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    if (param$exoeffects==1){
      
      textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("WSqx",names(dataIV))],collapse="+"))
      
    } else {
      
      textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("WSqx",names(dataIV))],collapse="+"))
      
    }
    
    panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
    panelest_coeftable <- panelest$coeftable
    panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
    panelest_1ststage <- fixest::fitstat(panelest,"ivwald")
    
    result$pop$panelest <- panelest
    result$pop$panelest_coeftable <- panelest_coeftable
    result$pop$panelest_etable <- panelest_etable
    result$pop$panelest_1ststage <- panelest_1ststage
    
    result$pop$rho <- panelest$coefficients["fit_Wy"]
    result$pop$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
    
    if (param$exoeffects==1){
      
      result$pop$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
      
    } else {
      
      result$pop$gamma <- 0*result$pop$beta
      
    }
    
  }
  
  # external IV ================================================================
  
  if (param$doexternaliv==1) {
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    sqnz <- names(data)[grep("[z]",names(data))]
    
    if (length(sqnz)>0){
      
      if (param$exoeffects==1){
        
        textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("Wz",names(dataIV))],collapse="+"))
        
      } else {
        
        textformula <- paste0("y ~ -1+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time|Wy~", paste(names(dataIV)[grep("Wz",names(dataIV))],collapse="+"))
        
      }
      
      panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
      panelest_coeftable <- panelest$coeftable
      panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
      panelest_1ststage <- fixest::fitstat(panelest,"ivwald")
      
      result$extiv$panelest_coeftable <- panelest_coeftable
      result$extiv$panelest_etable <- panelest_etable
      result$extiv$panelest_1ststage <- panelest_1ststage
      
      result$extiv$rho <- panelest$coefficients["fit_Wy"]
      result$extiv$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
      
      if (param$exoeffects==1){
        
        result$extiv$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
        
      } else {
        
        result$extiv$gamma <- 0*result$extiv$beta
        
      }
      
    }
    
  }
  
  # post OLS ===================================================================
  
  if (param$dopostols==1){
    
    dataIV = popIV(result$unpenalisedgmm$W,data,datasq,param)
    
    if (param$exoeffects==1){
      
      textformula <- paste0("y ~ -1+Wy+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"+", paste(names(dataIV)[grep("Wx",names(dataIV))],collapse="+"),"|id+time")
      
    } else {
      
      textformula <- paste0("y ~ -1+Wy+",paste0(names(dataIV)[grep("x",names(dataIV))][1:param$K],collapse="+"),"|id+time")
      
    }
    
    panelest <- fixest::feols(as.formula(textformula),data=dataIV,vcov="hetero")
    panelest_coeftable <- panelest$coeftable
    panelest_etable <- fixest::etable(panelest, signifCode = c("***"=0.01, "**"=0.05, "*"=0.10))
    
    result$ols$panelest_coeftable <- panelest_coeftable
    result$ols$panelest_etable <- panelest_etable
    
    result$ols$rho <- panelest$coefficients["Wy"]
    result$ols$beta <- panelest$coefficients[names(dataIV)[grep("x",names(dataIV))][1:param$K]] * result$sdxnorm
    
    if (param$exoeffects==1){
      
      result$ols$gamma  <- panelest$coefficients[names(dataIV)[grep("Wx",names(dataIV))]] * result$sdxnorm
      
    } else {
      
      result$ols$gamma <- 0*result$ols$beta
      
    }
    
  }
  
  # BIC ========================================================================
  
  value <- result$adaen$optim$value
  W <- result$adaen$W
  penalty <- result$adaen$penalty
  value_nopen <- value - (param$lambdaL1Star*sum(abs(W)*penalty,na.rm=TRUE) + param$lambdaL2*sum(W^2))
  
  BIC <- log(value_nopen) + sum(abs(result$adaen$W)>1e-5)*log(param$T)/param$T
  
  result$BIC <- BIC
  
  # cross-validation ===========================================================
  
  if (docv==1) {
    
    crossvalidation <- cv(data,param)
    result$crossvalidation <- crossvalidation
    
  }
  
  # return =====================================================================
  
  return(result)
  
}

# optimisation =================================================================

outerOptim = function(p,lambdaL1,lambdaL2,penalty,estmethod,rsn=0,datasq,param){
  
  # bounds
  
  if (param$exoeffects==1) {
    
    lower <- c(0,rep(-Inf,param$K),rep(-Inf,param$K))
    upper <- c(.95,rep(Inf,param$K),rep(Inf,param$K))
    
  } else {
    
    lower <- c(0,rep(-Inf,param$K))
    upper <- c(.95,rep(Inf,param$K))
    
  }
  
  if (param$opt$outermeth == 'L-BFGS-B') {
    
    o <- optimx::optimx(
      p,
      f,
      lower=lower,
      upper=upper,
      control=list(trace=0,maxit=param$maxit,kkt=FALSE,save.failures=TRUE,pgtol=1e-5),
      hessian=FALSE,
      method=param$opt$outermeth,
      itnmax=param$opt$outeriter,
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  if (param$opt$outermeth == 'newuoa') {
    
    o <- optimx::optimx(
      p,
      f,
      control=list(trace=0,kkt=FALSE,save.failures=TRUE),
      hessian=FALSE,
      method=param$opt$outermeth,
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  if (param$opt$outermeth == 'testing') {
    
    o <- optimx::optimx(
      p,
      f,
      method=c('Nelder-Mead','BFGS','L-BFGS-B','CG','nlm','nlminb','spg','ucminf','newuoa','bobyqa','nmkb','hjkb','Rcgmin'),
      lambdaL1=lambdaL1,
      lambdaL2=lambdaL2,
      penalty=penalty,
      estmethod=estmethod,
      rsn=rsn,
      datasq=datasq,
      param=param,
      returnW=0)
    
  }
  
  return(o)
  
}

f <- function(p,lambdaL1,lambdaL2,penalty,estmethod,rsn=0,datasq,param,returnW=0){
  
  # download parameters
  
  rho = p[1]
  beta = p[2:(param$K+1)]
  
  if (param$exoeffects==1){
    gamma = p[(param$K+2):(2*param$K+1)]
  } else {
    gamma = 0*beta
  }
  
  # penalisation parameters
  
  lambda <- lambdaL1+lambdaL2
  if (lambda==0){lambda <- 1}
  alpha <- lambdaL1/lambda
  
  # prepare data
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  # Y <- y-x*beta
  texteval <- paste0("Y<-datasq$sqy-",paste0("datasq$",sqn,"*beta[",1:length(sqn),"]",collapse="-"))
  eval(parse(text=texteval))
  
  # X <- t(rho*datasq$Y + gamma*datasq$X)
  evaltext <- paste0("X<-t(rho*datasq$sqy+",paste0("datasq$",sqn,"*gamma[",1:length(sqn),"]",collapse="+"),")")
  eval(parse(text=evaltext))
  
  # Z <- t(datasq$X)
  evaltext <- paste0("Z<-cbind(",paste(paste0("t(datasq$",sqn,")"),collapse=","),")")
  eval(parse(text=evaltext))
  
  # allocating space for W and the error term
  
  W <- param$Wfixed
  err <- matrix(ncol=0,nrow=1)
  
  # GLMNET - DECOUPLED
  
  if (estmethod=='glmnet'){
    
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      if (param$r>0) {
        
        ZYia <- rbind(ZYi,param$r)
        ZXia <- rbind(ZXi,param$r)
        
      } else {
        
        ZYia <- ZYi
        ZXia <- ZXi
        
      }
      
      if (length(jX)>=2){
        
        fit <- glmnet::glmnet(ZXia,ZYia,
                              lambda=lambda,
                              alpha=alpha,
                              penalty.factor=peni,
                              lower.limits=0,
                              upper.limits=1,
                              intercept=FALSE,
                              standardize=FALSE)
        
        fit <- as.numeric(fit$beta)
        
        if (rsn==1){fit <- fit / sum(fit)}
        
        W[i,jX] <- fit
        
        err = c(err,ZYi-ZXi%*%fit)
        
      }
      
      if (length(jX)==1){
        
        W[i,jX] <- 1
        
        err = c(err,ZYi-ZXi*1)
        
      }
      
      if (length(jX)==0){
        
        err = c(err,ZYi)
        
      }
      
    }
    
  }
  
  # GLMNET - COUPLED
  
  if (estmethod=='glmnet-coupled'){
    
    # allocate space for ZY and ZX
    
    ZY <- matrix(ncol=1,nrow=0)
    ZX <- matrix(ncol=0,nrow=0)
    pen <- matrix(ncol=1,nrow=0)
    
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      if (param$r > 0) {
        ZYia <- rbind(ZYi,param$r)
        ZXia <- rbind(ZXi,param$r)
      } else {
        ZYia <- ZYi
        ZXia <- ZXi
      }
      
      ZY <- rbind(ZY,ZYia) ## (N^2x1)
      ZX <- Matrix::bdiag(ZX,ZXia) ## (N^{2}XN*)
      pen <- rbind(pen,peni)
      
    }
    
    fit <- glmnet::glmnet(ZX,ZY,lambda=lambda,alpha=alpha,penalty.factor=pen,lower.limits=0,upper.limits=1,intercept=FALSE)
    
    fit = as.numeric(fit$beta)
    
    fit = (1+param$lambdaL2/(param$T^2))*fit ## Debiasing Step
    
    # put coefficients back in W
    W <- t(param$Wfixed)
    W[which(is.na(W))] <- as.numeric(fit)
    W <- t(W)
    
    ## row-sum normalize
    
    if (rsn==1) {
      
      rsnW = rowSums(W)
      
      for (i in 1:param$N){
        
        if (rsnW[i] > 0){W[i,] = W[i,]/rsnW[i]}else{W[i,] = W[i,]/1}
        
      }
      
    }
    
    # compute error and prepare return
    
    err = ZY-ZX%*%fit
    
    f_obj = sum(err^2) + lambdaL1*sum(abs(W)*penalty,na.rm=TRUE) + lambdaL2*sum(W^2)
    
    # prepare return
    
    if (returnW==0){
      
      return(f_obj)
      
    } else {
      
      return(W)
      
    }
    
  }
  
  # CVXR
  
  if (estmethod=='cvxr'){
    
    pen <- matrix(ncol=1,nrow=0)    
    ZY <- matrix(ncol=1,nrow=0)
    ZX <- matrix(ncol=0,nrow=0)
    
    theta = CVXR::Variable(sum(is.na(param$Wfixed)))
    
    length_optim_Wfixed = rowSums(is.na(param$Wfixed))
    stop = c(1,cumsum(length_optim_Wfixed))
    
    constraints = list()
    
    # theta_pos = rowSums(is.na(param$Wfixed))
    # theta_pos = c(0,cumsum(theta_pos))
    
    for (i in 1:param$N){
      
      Yi <- as.matrix(Y[i,])
      
      jX <- which(is.na(param$Wfixed[i,]))
      Xi <- X[,jX]
      peni <- as.matrix(penalty[i,jX])
      pen <- rbind(pen,peni)
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      ZY <- rbind(ZY,ZYi) ## (N^2x1)
      ZX <- Matrix::bdiag(ZX,ZXi) ## (N^{2}XN*)
      
      # st <- theta_pos[i]+1
      # en <- theta_pos[i+1]
      # 
      # constraints <- append(constraints,sum(theta[st:en])==1)
      
      if (i == 1){constraints = append(constraints,(sum(theta[stop[i]:stop[i+1]]) == 1))}
      if (i > 1){constraints = append(constraints,(sum(theta[(stop[i]+1):stop[i+1]]) == 1))}
      
    }
    
    constraints = append(constraints,(theta >= rep(0,sum(is.na(param$Wfixed)))))
    
    loss <- sum((ZY - ZX %*% theta)^2) / (2 * dim(ZX)[2])
    
    obj <- loss + elastic_reg(theta,lambda=lambda,alpha=alpha,pen=pen)
    
    prob <- CVXR::Problem(CVXR::Minimize(obj),constraints=constraints)
    
    result <- CVXR::psolve(prob,solver="ECOS")
    
    fit = result$getValue(theta)
    
    # put coefficients back in W
    
    W <- t(param$Wfixed)
    W[which(is.na(W))] <- as.numeric(fit)
    W <- t(W)
    
    # objective function
    
    err = ZY-ZX%*%fit
    
  }
  
  # NO ESTIMATION METHOD - GIVEN W
  
  if (estmethod=='givenw'){
    
    for (i in 1:param$N){
      
      Yi <- Y[i,]
      
      jX <- setdiff(seq(1:param$N),i)
      
      Xi <- X[,jX]
      
      ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
      ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
      
      fit <- W[i,jX]
      
      err = c(err,ZYi-ZXi%*%fit)
      
    }
    
  }
  
  # compute the final objective function
  
  f_obj = sum(err^2) + lambdaL1*sum(abs(W)*penalty,na.rm=TRUE) + lambdaL2*sum(W^2)
  
  # prepare return
  
  if (returnW==0){
    
    return(f_obj)
    
  } else {
    
    return(W)
    
  }
  
}

elastic_reg <- function(theta,lambda,alpha,pen=pen) {
  
  ridge <- (1 - alpha) / 2 * sum(theta^2)
  lasso <- alpha * CVXR::p_norm(pen*theta, 1)
  l <- lambda * (lasso + ridge)
  return(l)
  
}

bic <- function(p,W,datasq,param){
  
  # download parameters
  
  rho = p[1]
  beta = p[2:(param$K+1)]
  
  if (param$exoeffects==1){
    gamma = p[(param$K+2):(2*param$K+1)]
  } else {
    gamma = 0*beta
  }
  
  # prepare data
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  # Y <- y-x*beta
  texteval <- paste0("Y<-datasq$sqy-",paste0("datasq$",sqn,"*beta[",1:length(sqn),"]",collapse="-"))
  eval(parse(text=texteval))
  
  # X <- t(rho*datasq$Y + gamma*datasq$X)
  evaltext <- paste0("X<-t(rho*datasq$sqy+",paste0("datasq$",sqn,"*gamma[",1:length(sqn),"]",collapse="+"),")")
  eval(parse(text=evaltext))
  
  # Z <- t(datasq$X)
  evaltext <- paste0("Z<-cbind(",paste(paste0("t(datasq$",sqn,")"),collapse=","),")")
  eval(parse(text=evaltext))
  
  # allocate space for ZY and ZX
  
  ZY <- matrix(ncol=1,nrow=0)
  ZX <- matrix(ncol=0,nrow=0)
  pen <- matrix(ncol=1,nrow=0)
  fit <- matrix(ncol=1,nrow=0)
  
  for (i in 1:param$N){
    
    Yi <- Y[i,]
    
    jX <- c(1:param$N)
    Xi <- X[,jX]
    
    ZYi <- t(Z)%*%Yi ## (NxT)x(Tx1) = (Nx1)
    ZXi <- t(Z)%*%Xi ## (NxT)x(TxN*) = (NxN*_i)
    
    ZY <- rbind(ZY,ZYi) ## (N^2x1)
    ZX <- Matrix::bdiag(ZX,ZXi) ## (N^{2}XN*)
    
    fiti <- W[i,jX]
    fit <- c(fit,fiti)
    
  }
  
  # compute error and prepare return
  
  err = ZY-ZX%*%fit
  
  f_obj = sum(err^2)
  
  return(f_obj)
  
}

# data management ==============================================================

organiseData <- function(data,param){
  
  # sort data by id and time 
  
  iX <- base::order(x<-data$time,y<-data$id)
  data <- data[iX,]
  
  uId <- base::sort(unique(data$id))
  uTime <- base::sort(unique(data$time))
  
  id <- NA * data$id
  time <- NA * data$time
  for (i in 1:length(uId)) {id[which(data$id==uId[i])] <- i}
  for (t in 1:length(uTime)) {time[which(data$time==uTime[t])] <- t}
  
  data$id <- id
  data$time <- time
  
  # remove fixed effects 
  
  iX <- grep("[xyz]",colnames(data))
  
  data_dm <- fixest:::demean(X=data[,iX],f=data[,c("id","time")])
  data[,iX] <- data_dm
  
  iXt <- which(data$time==min(data$time))
  data <- data[-iXt,]
  data$time <- data$time-1
  
  # standardise
  
  if (param$standardise==1){
    
    sdx <- matrix(NA,nrow=1,ncol=dim(data)[2])
    colnames(sdx) <- colnames(data)
    
    for (i in iX){
      
      x <- data[,i]
      sdx[,i] <- sqrt(var(data[,i]))
      data[,i] <- x/as.vector(sdx[,i])
      
      sdx_y <- sdx[which(colnames(sdx)=="y")]
      sdx_x <- sdx[grep("[x]",colnames(sdx))]
      sdxnorm <- sdx_y * (sdx_x^(-1))
      
    }
    
  } else {
    
    sdx <- 1
    sdxnorm <- 1
    
  }
  
  # add weights
  
  if (length(param$timeweights)>1) {
    
    for (i in iX){
      
      data[,i] <- data[,i]*sqrt(rep(param$timeweights,each=param$N))
      
    }
    
  }
  
  # square the data
  
  datalist <- list()
  
  datalist$datasq <- squareData(data,param)
  datalist$data <- data
  
  datalist$sdx <- sdx
  datalist$sdxnorm <- sdxnorm
  
  return(datalist)
  
}

squareData <- function(data,param){
  
  cnames <- colnames(data)[grep("[xyz]",colnames(data))]
  
  datasq <- list()
  
  for (i in 1:length(cnames)) {
    
    txteval <- paste0("datasq$sq",cnames[i], " = base::as.matrix(tidyr::pivot_wider(data,id_cols=id,names_from=time,values_from=",cnames[i],")[,-1])")
    eval(parse(text=txteval))
    
    
  }
  
  return(datasq)
  
}

popIV <- function(What,data,datasq,param) {
  
  # creates peer-of-peers instrument
  
  WhatSq <- What %*% What
  
  Wy_pop <- matrix(0,nrow=param$N*(param$T-1))
  #Wx_pop <- matrix(0,nrow=length(data$x))
  sqnx <- names(data)[grep("[x]",names(data))]
  
  for (i in 1:length(sqnx)){
    texteval <- paste0("W",sqnx[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
    eval(parse(text=texteval))
    texteval <- paste0("WSq",sqnx[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
    eval(parse(text=texteval))
  }
  
  #Wz_pop <- matrix(0,nrow=length(data$x))
  sqnz <- names(data)[grep("[z]",names(data))]
  
  if (length(sqnz)>0) {
    
    for (i in 1:length(sqnz)){
      texteval <- paste0("W",sqnz[i],"_pop <- matrix(0,nrow=param$N*(param$T-1))")
      eval(parse(text=texteval))
    }
    
  }
  
  for (t in 1:(param$T-1)) {
    
    # select the observations for time period t
    Wy_pop[which(data$time==t),] <- What %*% datasq$sqy[,t]
    
    for (i in 1:length(sqnx)) {
      texteval <- paste0("W",sqnx[i],"_pop[which(data$time==t),] <- What %*% datasq$sq",sqnx[i],"[,t]")
      eval(parse(text=texteval))
      texteval <- paste0("WSq",sqnx[i],"_pop[which(data$time==t),] <- WhatSq %*% datasq$sq",sqnx[i],"[,t]")
      eval(parse(text=texteval))
    }
    
    
    if (length(sqnz)>0){
      for (i in 1:length(sqnz)) {
        texteval <- paste0("W",sqnz[i],"_pop[which(data$time==t),] <- What %*% datasq$sq",sqnz[i],"[,t]")
        eval(parse(text=texteval))
      }
    }
    
  }
  
  data$Wy <- Wy_pop
  for (i in 1:length(sqnx)) {
    texteval <- paste0("data$W",sqnx[i]," <- W",sqnx[i],"_pop")
    eval(parse(text=texteval))
    texteval <- paste0("data$WSq",sqnx[i]," <- WSq",sqnx[i],"_pop")
    eval(parse(text=texteval))
  }
  if (length(sqnz)>0){
    for (i in 1:length(sqnz)) {
      texteval <- paste0("data$W",sqnz[i]," <- W",sqnz[i],"_pop")
      eval(parse(text=texteval))
    }
  }
  
  return(data)
  
}

# moment cond derivative =======================================================

computeDerivatives <- function(data,datasq,param){
  
  # compute beta hat without considering the spatial effects
  formula <- paste0("y~",paste(colnames(data)[grep("[x]",colnames(data))],collapse = "+"),"|id+time")
  panelest <- fixest::feols(as.formula(formula),data)
  betahat <- base::as.vector(panelest$coefficients)
  
  # compute the derivative starting from the empty network, estimatate beta as above, rho = .5 and gamma = 0
  Wstart <- base::matrix(0,param$N,param$N)
  
  #if (param$ncores==0){
  
  # IF SINGLE-CORE THREAD
  g <- momentCond2(rho=.5,beta=betahat,gamma=rep(0,length(betahat)),W=Wstart,datasq,param)
  WDer <- base::lapply(1:(param$N^2),momentDer_i,rho=.5,
                       beta=betahat[1],gamma=0,W=Wstart,datasq=datasq,param=param,g=g)
  WDer <- base::unlist(WDer)
  dim(WDer) <- c(param$N,param$N)
  
  #}
  
  r <- unlist(WDer)
  
  # store
  param$WDer <- WDer
  return(WDer)
  
}

momentCond2 <- function(rho,beta,gamma,W,datasq,param) {
  
  tt <- 1:(dim(datasq$sqy)[2])
  g <- base::sapply(tt,function(t) momentCond2_aux(rho=rho,beta=beta,gamma=gamma,W=W,datasq=datasq,param=param,t=t) )
  g <- Matrix::rowSums(g)/(dim(datasq$sqy)[2])
  
  return(g)
  
}


momentCond2_aux <- function(rho,beta,gamma,W,datasq,param,t){
  
  # select the X
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  xt <- sapply(sqn, function(i) eval(parse(text=paste0("datasq$",i,"[,t]"))))
  sqn <- names(datasq)[grep("[z]",names(datasq))]
  zt <- xt
  dim(zt) <- c(dim(zt)[1]*dim(zt)[2],1)
  
  # calculate the residuals
  residT <- (diag(dim(W)[1]) - rho*W) %*% datasq$sqy[,t] - xt %*% beta - W %*% xt %*% gamma
  
  # compute the moment condition
  g <- base::kronecker(zt,residT)
  
  return(g)
  
}

momentDer_i <- function(rho,beta,gamma,W,datasq,param,g,i) {
  
  if (is.na(param$Wfixed[i])) {
    
    # create matrices used in the computations
    I <- diag(param$N)
    F1 <- solve(I-rho*W)
    F2 <- beta*I+gamma*W
    F3 <- F1*W
    
    dDer <- matrix(0,param$N^2,1)
    A <- matrix(0,param$N,param$N)
    A[i] <- 1
    
    xT <- as.matrix(datasq$sqx1)
    Gk <- (rho * F1 %*% A %*% F1 %*% F2 + gamma * F1 %*% A) %*% xT
    
    for (t in 1:(dim(datasq$sqy)[2])) {
      
      gt <- as.matrix(kronecker(xT[,t],Gk[,t]))
      
      dDer <- dDer + gt
      
    }
    
    dDer <- dDer / (dim(datasq$sqy)[2])
    
    g <- g[1:(length(g)/param$K)]
    WDer <- -2 * t(g) %*% dDer
    
  } else {
    
    WDer <- 0
    
  }
  
  return(WDer)
  
}


# asymptotic distribution ==================================================

momentDer_asyvar <- function(rho,beta,gamma,W,datasq,param,s,i) {
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  dDer <- matrix(0,param$N^2*param$K,1)
  
  for (t in 1:(dim(datasq$sqy)[2])) {
    
    yT <- as.matrix(datasq$sqy)
    yT <- yT[,t]
    
    if (s=="rho"){
      
      Gk <- -1 * W %*% yT
      
    }
    
    if (s=="beta"){
      
      evaltext <- paste0("xT <- as.matrix(datasq$sqx",i,")")
      eval(parse(text=evaltext))
      xT <- xT[,t]
      
      Gk <- -xT
      
    }
    
    if (s=="gamma"){
      
      evaltext <- paste0("xT <- as.matrix(datasq$sqx",i,")")
      eval(parse(text=evaltext))
      xT <- xT[,t]
      
      Gk <- -1 * W %*% xT
      
    }
    
    if (s=="W"){
      
      A <- matrix(0,param$N,param$N)
      A[i] <- 1
      
      sqn <- names(datasq)[grep("[x]",names(datasq))]
      
      texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
      eval(parse(text=texteval))
      
      Gk <- - rho * A %*% yT - A %*% xT %*% gamma
      
    }
    
    texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
    eval(parse(text=texteval))
    
    dim(xT) <- c(dim(xT)[1]*dim(xT)[2],1)
    
    gt <- as.matrix(kronecker(xT,Gk))
    
    dDer <- dDer + gt
    
  }
  
  dDer <- dDer / dim(datasq$sqy)[2]
  
  return(dDer)
  
}

moment_asyvar <- function(rho,beta,gamma,W,datasq,param,i) {
  
  sqn <- names(datasq)[grep("[x]",names(datasq))]
  
  d <- matrix(0,param$N^2*param$K,param$N^2*param$K)
  
  for (t in 1:(dim(datasq$sqy)[2])) {
    
    texteval <- paste0("xT<-cbind(",paste0("datasq$",sqn,"[,t]",collapse=","),")")
    eval(parse(text=texteval))
    
    yT <- as.matrix(datasq$sqy)
    yT <- yT[,t]
    
    Gk <- yT - rho*W%*%yT - xT%*%beta - W%*%xT%*%gamma
    
    dim(xT) <- c(dim(xT)[1]*dim(xT)[2],1)
    
    gt <- as.matrix(kronecker(xT,Gk))
    
    d <- d + gt %*% t(gt)
    
  }
  
  d <- d / dim(datasq$sqy)[2]
  
  return(d)
  
}


# cross-validation =============================================================

cv <- function(data,param){
  
  # compute to different folds
  
  l <- sapply(1:4,function(f) {cvinner(data,param,fold=f)})
  return(l)
  
}


cvinner <- function(data,param,fold) {
  
  # split points
  
  splitpoints <- seq(from=min(data$time),to=(max(data$time)+1),length.out=6)
  splitpoints <- splitpoints[-1]
  
  st <- splitpoints[fold]
  en <- splitpoints[fold+1]
  
  # slice the data
  
  data_estimation <- data %>% filter(time<st)
  
  data_evaluation <- data %>% filter(time>=st & time<en)
  data_estimation$time <- data_estimation$time-min(data_estimation$time)+1 
  
  # run procedure
  
  lambda <- c(param$lambdaL1,param$lambdaL2,param$lambdaL1Star)
  res <- recoverNetwork(data_estimation,lambda,docv=0)
  
  # compute msqpe
  
  W <- res$adaen$W
  rho <- res$adaen$rho
  beta <- res$adaen$beta
  gamma <- res$adaen$gamma
  
  # compute msqpe 
  
  msqpe <- 0
  
  for (t in unique(data_evaluation$time)) {
    
    yt <- data_evaluation %>% filter(time==t) %>% dplyr::select(y)
    xt <- data_evaluation %>% filter(time==t) %>% dplyr::select(starts_with("x"))
    
    sqn <- names(data)[grep("[x]",names(data))]
    
    fitt <- rho * W %*% as.matrix(yt) + as.matrix(xt) %*% as.matrix(beta) + W %*% as.matrix(xt) %*% as.matrix(gamma)
    errort <- sum((as.matrix(yt)-fitt)^2)
    
    msqpe <- msqpe + errort
    
  }
  
  return(msqpe)
  
}

# generates simulated data =====================================================

gendata <- function(setting,seed){
  
  # SET SEED
  base::set.seed(seed)
  
  # load all possible settings
  allsettings <- loadsettings()
  param <- allsettings[allsettings$setting==setting,]
  param$setting <- setting
  
  if (setting<=999){
    
    data <- gendatai(param) 
    
  } else if (setting>1000 & setting<=1999) {
    
    # map the param vectors
    param_start <- param
    param_start$W0 <- strsplit(param$W0,"-")[[1]][1]
    
    param_end <- param
    param_end$W0 <- strsplit(param$W0,"-")[[1]][2]
    
    # generate the data
    data_start <- gendatai(param_start)
    data_end <- gendatai(param_end)
    data_end$time <- data_end$time + max(data_start$time)
    
    # combine
    data <- rbind(data_start,data_end)
    
  } else if (setting>2000) {
    
    data <- gendatai(param) 
    
  }
  
  return(data)
  
}

# loads the simulation parameters ==============================================

gendatai <- function(param){
  
  setting <- param$setting
  
  # create fixed effects
  if (param$FEsetting == "standard") {
    FE <- stats::rnorm(param$N)+1
    TE <- stats::rnorm(param$N)+1
  }
  
  
  # create fixed effects
  if (param$FEsetting == "zero") {
    FE <- 0*stats::rnorm(param$N)+1
    TE <- 0*stats::rnorm(param$N)+1
  }
  
  # define variance-covariance matrix
  if (param$varcov == "orthogonal") {
    Sigma <- base::diag(param$N)
  }
  
  if (param$varcov == "q3") {
    q <- 0.3
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q5") {
    q <- 0.5
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q8") {
    q <- 0.8
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (param$varcov == "q10") {
    q <- 1.0
    Sigma <- (1-q)*base::diag(param$N)+q
  }
  
  if (setting <= 2000){
    
    # create the reduced-form matrix
    base::eval(base::parse(text=paste('W0 <- W0$W0_', param$W0, '_N',param$N,sep='')))
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  } else {
    
    base::eval(base::parse(text=paste('W0erdos <- W0$W0_', 'erdosrenyi', '_N',param$N,sep='')))
    
    W0bilateral <- 0*W0erdos
    for (i in seq(1,29,by=2)){
      W0bilateral[i,i+1]<-1
      W0bilateral[i+1,i]<-1
    }
    
    mix <- (setting-2001)/5
    
    W0 <- (1-mix)*W0erdos + mix*(W0bilateral)
    
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  }
  
  
  # generate data for each time period
  data = dplyr::tibble(y=double(),x=double(),z=double(),id=integer(),time=integer())
  
  for (t in 1:param$T) {
    
    idt <- 1:param$N
    timet <- base::rep(t,param$N)
    
    et <- MASS::mvrnorm(n=1,mu=rep(0,param$N),Sigma)
    xt <- stats::rnorm(param$N*param$K)
    dim(xt) <- c(param$N,param$K)
    
    yt <- Pi0inv %*% xt %*% rep(param$beta,param$K) + Pi0inv %*% W0 %*% xt %*% rep(param$gamma,param$K) + Pi0inv %*% FE + Pi0inv %*% TE + Pi0inv %*% et
    xt <- as_tibble(as.data.frame(xt))
    names(xt) <- paste0('x',1:ncol(xt))
    
    datat <- dplyr::bind_cols(dplyr::tibble(id=idt,time=timet),dplyr::tibble(y=as.vector(yt)),xt)
    data <- base::rbind(data,datat)
    
  }
  
  # return data
  return(data)
  
}

# loads the simulation parameters ==============================================

loadsettings <- function(setting){
  
  K <- 1
  
  W0 <- 'erdosrenyi'
  N <- 15
  param <- data.frame(setting=1,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform',stringsAsFactors=FALSE)
  param <- base::rbind(param,data.frame(setting=0002,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0003,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0004,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0005,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0006,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0007,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0008,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0009,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0501,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0502,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0503,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0801,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0802,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0803,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0011,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0012,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0013,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0014,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0015,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0016,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0017,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0018,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0019,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0511,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0512,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0513,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0811,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0812,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0813,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0601,W0=W0,N=N,T=150,K=K,rho=0.1,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0602,W0=W0,N=N,T=150,K=K,rho=0.5,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0603,W0=W0,N=N,T=150,K=K,rho=0.7,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0604,W0=W0,N=N,T=150,K=K,rho=0.9,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0611,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.2,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0612,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.6,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0621,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.3,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0622,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.7,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0631,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q3',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0632,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q5',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0633,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q8',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0634,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q10',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=2001,W0='100erdosrenyi-000bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2002,W0='080erdosrenyi-020bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2003,W0='060erdosrenyi-040bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2004,W0='040erdosrenyi-060bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2005,W0='020erdosrenyi-080bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=2006,W0='000erdosrenyi-100bilateral',N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  
  N <- 50
  param <- base::rbind(param,data.frame(setting=0021,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0022,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0023,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0024,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0025,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0026,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0027,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0028,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0029,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0521,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0522,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0523,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0821,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0822,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0823,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0031,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0032,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0033,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0034,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0035,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0036,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0037,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0038,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0039,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'politicalparty'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0041,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0042,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0043,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0044,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0045,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0046,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0047,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0048,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0049,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0541,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0542,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0543,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0841,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0842,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0843,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0051,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0052,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0053,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0054,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0055,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0056,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0057,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0058,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0059,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0551,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0552,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0553,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0851,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0852,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0853,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0651,W0=W0,N=N,T=150,K=K,rho=0.1,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0652,W0=W0,N=N,T=150,K=K,rho=0.5,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0653,W0=W0,N=N,T=150,K=K,rho=0.7,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0654,W0=W0,N=N,T=150,K=K,rho=0.9,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0661,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.2,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0662,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.6,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0671,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.3,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0672,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.7,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0681,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q3',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0682,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q5',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0683,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q8',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0684,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='q10',Wknown='no',timeweights='uniform'))
  
  N <- 50
  param <- base::rbind(param,data.frame(setting=0061,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0062,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0063,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0064,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0065,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0066,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0067,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0068,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0069,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0561,W0=W0,N=N,T=500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0562,W0=W0,N=N,T=1000,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0563,W0=W0,N=N,T=1500,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=0861,W0=W0,N=N,T=50,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0862,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0863,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='zero',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0071,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0072,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0073,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0074,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0075,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0076,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0077,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0078,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0079,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'highschool'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0081,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0082,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0083,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0084,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0085,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0086,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0087,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0088,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0089,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0091,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0092,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0093,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0094,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0095,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0096,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0097,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0098,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0099,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 50
  param <- base::rbind(param,data.frame(setting=0101,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0102,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0103,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0104,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0105,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0106,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0107,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0108,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0109,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0111,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0112,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0113,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0114,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0115,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0116,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0117,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0118,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0119,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'duflo'
  N <- 15
  param <- base::rbind(param,data.frame(setting=0121,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0122,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0123,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0124,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0125,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0126,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0127,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0128,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0129,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 30
  param <- base::rbind(param,data.frame(setting=0131,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0132,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0133,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0134,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0135,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0136,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0137,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0138,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0139,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 50
  param <- base::rbind(param,data.frame(setting=0141,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0142,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0143,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0144,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0145,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0146,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0147,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0148,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0149,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  N <- 70
  param <- base::rbind(param,data.frame(setting=0151,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0152,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0153,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0154,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0155,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0156,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0157,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0158,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=0159,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  ##
  
  W0 <- 'erdosrenyi-politicalparty'
  N <- 30
  
  param <- base::rbind(param,data.frame(setting=1011,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1012,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1013,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1014,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1015,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1016,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1017,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1018,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1019,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=1011+10,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1012+10,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1013+10,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1014+10,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1015+10,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1016+10,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1017+10,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1018+10,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1019+10,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  
  param <- base::rbind(param,data.frame(setting=1011+20,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1012+20,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1013+20,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1014+20,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1015+20,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1016+20,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1017+20,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1018+20,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1019+20,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  
  W0 <- 'highschool-duflo'
  N <- 70
  
  param <- base::rbind(param,data.frame(setting=1111,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1112,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1113,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1114,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1115,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1116,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1117,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1118,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  param <- base::rbind(param,data.frame(setting=1119,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  param <- base::rbind(param,data.frame(setting=1111+10,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1112+10,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1113+10,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1114+10,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1115+10,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1116+10,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1117+10,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1118+10,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  param <- base::rbind(param,data.frame(setting=1119+10,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='start'))
  
  param <- base::rbind(param,data.frame(setting=1111+20,W0=W0,N=N,T=005,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1112+20,W0=W0,N=N,T=010,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1113+20,W0=W0,N=N,T=015,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1114+20,W0=W0,N=N,T=025,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1115+20,W0=W0,N=N,T=050,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1116+20,W0=W0,N=N,T=075,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1117+20,W0=W0,N=N,T=100,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1118+20,W0=W0,N=N,T=125,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  param <- base::rbind(param,data.frame(setting=1119+20,W0=W0,N=N,T=150,K=K,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='end'))
  
  W0 <- 'geoneighuniform'
  param <- base::rbind(param,data.frame(setting=901,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration'
  param <- base::rbind(param,data.frame(setting=902,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration_trim1'
  param <- base::rbind(param,data.frame(setting=903,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighmigration_trim2'
  param <- base::rbind(param,data.frame(setting=904,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighlength'
  param <- base::rbind(param,data.frame(setting=905,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  W0 <- 'geoneighdecay'
  param <- base::rbind(param,data.frame(setting=906,W0=W0,N=48,T=53,K=1,rho=0.3,beta=0.4,gamma=0.5,FEsetting='standard',varcov='orthogonal',Wknown='no',timeweights='uniform'))
  
  return(param)
  
}

# true coefficient =============================================================

truecoeffs <- function(param){
  
  # # load all possible settings
  # allsettings <- loadsettings()
  # param <- allsettings[allsettings$setting==setting,]
  
  setting <- param$setting
  
  # create the reduced-form matrix
  if (setting <= 2000){
    
    # create the reduced-form matrix
    base::eval(base::parse(text=paste('W0 <- W0$W0_', param$W0, '_N',param$N,sep='')))
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  } else {
    
    base::eval(base::parse(text=paste('W0erdos <- W0$W0_', 'erdosrenyi', '_N',param$N,sep='')))
    
    W0bilateral <- 0*W0erdos
    for (i in seq(1,29,by=2)){
      W0bilateral[i,i+1]<-1
      W0bilateral[i+1,i]<-1
    }
    
    mix <- (setting-2001)/5
    
    W0 <- (1-mix)*W0erdos + mix*(W0bilateral)
    
    Pi0inv <- base::solve(diag(param$N)-param$rho*W0)
    
  }
  
  # build list to return
  true <- list()
  true$N <- param$N
  true$T <- param$T
  true$rho <- param$rho
  true$beta <- param$beta
  true$gamma <- param$gamma
  true$W0 <- W0
  true$Pi0inv <- Pi0inv
  
  return(true)
  
}



