

### adjacent categories category-specific effects

library("fastmatrix")


####################  programs ##############
###############################
fitadjcats<-function(resp,k,pred,deriv,hessian,maxit,start, lambda,weightlasso,single){
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  
  #### if lambda >0 lasso fusion
  #### weightlasso means on lasso smoothing lambda*weightlasso (typical 0.3)
  #### single (set) specifies variables that are considered global (set lambda<-1000, weightlasso=0) 
     # if single==0 all are category specific
     # if single == "smooth" smoothing of effects across categories
  
  p<-dim(pred)[2]
  
  ## single
  singleset<-single
  if (single[1] == 0) singleset<-rep(1:p)
  
  betamat<-matrix(.01,p,k-1)  ## matrix betas p x k-1
  #betamat<-matrix(1:((k-1)*p),p,k-1)
  beta<-c(betamat)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  #matrix(beta,p,k-1)
  
  
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv !='der')fitopt <- optim(parst, loglikadjcatcats, gr = NULL,resp=resp,k=k,pred=pred,lambda=lambda,weightlasso=weightlasso, 
                                      singleset=singleset,method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt<- optim(parst, loglikadjcatcats, gr = derloglikadjcatcats,resp=resp,k=k,pred=pred,lambda=lambda,
                                     weightlasso=weightlasso,singleset=singleset, method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  
   loglikunpen<- -loglikadjcatcats(fitopt$par,resp,k,pred, lambda=0,weightlasso,singleset=0)
    
  AIC<- -2*(loglikunpen-length(parst))
  
  #derloglikadjcatcats(fitopt$par,resp,k,pred)
  #parst<-fitopt$par
  
  stderr<-0
  stddisp<-0
  zval<-0
  stdmat<-0
  
  #location<-fitopt$par[1:p]
  beta<-fitopt$par[1:((k-1)*p)]
  location<-matrix(beta,p,k-1)
  thresh<-fitopt$par[((k-1)*p+1):length(fitopt$par)]
  
    if(hessian ==TRUE){
      hessinv<-solve(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  
  stdpar<-std[1:((k-1)*p)]  ### first are covariate weights
  stdmat<-matrix(stdpar,p,k-1)
  #zval<-fitopt$par[1:p]/stderr
  zval<-location/stdmat
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  #location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= loglikunpen,"Pen Loglik"= -fitopt$value,"AIC"=AIC, "parametermatrix"=location,  
                      "parameter"=fitopt$par, "Intercepts"=thresh,  "Stderr"=stdmat, "zval"=zval, "convergence"= fitopt$convergence)  
  return(newList)}


##################################
loglikadjcatcats<-function(par,resp,k,pred, lambda,weightlasso,singleset){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:((k-1)*p)]  ### first are covariate weights
  betamat<-matrix(beta,p,k-1)
  beta0<-c(0,par[((k-1)*p+1):((k-1)*p+k-1)])   ### thresholds
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%betamat[,r-1])}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum}
    sum(prob[i,])
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  #### penalty fusion L_1
  pen<-0
  c<-.01
  if(lambda>0 ){
      for (v in singleset){
      for (r1 in 1:(k-2)){
        #for (r2 in (r1+1):(k-1)) 
       pen<- pen + NormShifto((betamat[v,r1]-betamat[v,r1+1]),c)
         }
      }
  loglik<-loglik- lambda*pen
   } ### end lambd
  
  #### penalty lasso
  pen<-0
  c<-.01
  if(lambda>0){
    for (v in 1:p){
              pen<- pen + (betamat[v, ]%*%as.matrix(betamat[v, ])+c)^(1/2)
          }
    loglik<-loglik- weightlasso*lambda*pen   ### factor!
  } ### end lambd
  
  
  
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  #return(newList)
  return(loglikneg)
  }

#######################################################

##################################
derloglikadjcatcats<-function(par,resp,k,pred, lambda,weightlasso,singleset){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  #p<-dim(pred)[2]
  #beta<-par[1:p]  ### first are covariate weights
  #beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  
  p<-dim(pred)[2]
  beta<-par[1:((k-1)*p)]  ### first are covariate weights
  betamat<-matrix(beta,p,k-1)
  betamatext<-cbind(rep(0,p),betamat)   ### extended
  beta0<-c(0,par[((k-1)*p+1):((k-1)*p+k-1)])   ### thresholds
  
  
  
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,k-1)  ## matrix
  der2<-matrix(0,(k-1),1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%betamat[,r-1])}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    #sum(prob[i,])
    
    cat<-resp[i]
    
    ### covariates
    #der1<-der1+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der1[,j-1]<-der1[,j-1]+(ind1-sum(prob[i,j:k]))*as.matrix(pred[i,])
    }
    
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
  }  ## end i
  
  #### penalty fusion
  #pen<-0
  c<-.01
  if(lambda>0){
    
    for (s in 2:k){
        for (j in singleset){
          if(s==2)  der1[j,1]<-der1[j,1]-lambda*NormShiftoder((betamat[j,1]-betamat[j,2]),c)
          if(s==k)  der1[j,k-1]<-der1[j,k-1]+-lambda*NormShiftoder((betamat[j,k-2]-betamat[j,k-1]),c)
          
          if((s >2) & (s <k))  {der1[j,s-1]<-der1[j,s-1]-lambda*NormShiftoder((betamat[j,s-1]-betamat[j,s]),c)
                                der1[j,s]<-der1[j,s]+lambda*NormShiftoder((betamat[j,s-1]-betamat[j,s]),c)}
         
              }} ### end k
    } # enf if
    
  #### penalty lasso
  #pen<-0
  c<-.01
  if(lambda>0){
    
    for (s in 2:k){
      for (j in 1:p){
        der1[j,s-1]<-der1[j,s-1]-weightlasso*lambda*(betamat[j, ]%*% as.matrix(betamat[j, ])+c)^(-1/2)*betamat[j,s-1 ]
      }} ### end k
  } # enf if
  
    der<-c(c(der1),der2)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

#####

###############################
fitadjcatsmooth<-function(resp,k,pred,deriv,hessian,maxit,start, lambda,weightlasso){
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  
  #### if lambda >0 lasso fusion
  #### weightlasso means on lasso smoothing lambda*weightlasso (typical 0.3)
  #### single (set) specifies variables that are considered global (set lambda<-1000, weightlasso=0) 
  # if single==0 all are category specific
  # if single == "smooth" smoothing of effects across categories
  
  p<-dim(pred)[2]
  
  ## single
  singleset<-rep(1:p)  ### all smooth
  
  betamat<-matrix(.01,p,k-1)  ## matrix betas p x k-1
  #betamat<-matrix(1:((k-1)*p),p,k-1)
  beta<-c(betamat)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  #matrix(beta,p,k-1)
  
  
  parst<-c(beta,beta0)
  if(sum(start!=0))parst<-start
  if (deriv !='der')fitopt <- optim(parst, loglikadjcatsmooth, gr = NULL,resp=resp,k=k,pred=pred,lambda=lambda,weightlasso=weightlasso, 
                                    singleset=singleset,method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv =='der')fitopt<- optim(parst, loglikadjcatsmooth, gr = derloglikadjcatsmooth,resp=resp,k=k,pred=pred,lambda=lambda,
                                   weightlasso=weightlasso,singleset=singleset, method = "BFGS",
                                   lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  
  loglikunpen<- -loglikadjcatcats(fitopt$par,resp,k,pred, lambda=0,weightlasso,singleset=0)
  
  AIC<- -2*(loglikunpen-length(parst))
  
  #derloglikadjcatcats(fitopt$par,resp,k,pred)
  #parst<-fitopt$par
  
  stderr<-0
  stddisp<-0
  zval<-0
  stdmat<-0
  
  #location<-fitopt$par[1:p]
  beta<-fitopt$par[1:((k-1)*p)]
  location<-matrix(beta,p,k-1)
  thresh<-fitopt$par[((k-1)*p+1):length(fitopt$par)]
  
  if(hessian ==TRUE){
    hessinv<-solve(fitopt$hessian)
    std<-sqrt(diag(hessinv))  
    stderr<-std[1:p]
    
    stdpar<-std[1:((k-1)*p)]  ### first are covariate weights
    stdmat<-matrix(stdpar,p,k-1)
    #zval<-fitopt$par[1:p]/stderr
    zval<-location/stdmat
    
    #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
    #location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= loglikunpen,"Pen Loglik"= -fitopt$value,"AIC"=AIC, "parametermatrix"=location,  
                  "parameter"=fitopt$par, "Intercepts"=thresh,  "Stderr"=stdmat, "zval"=zval, "convergence"= fitopt$convergence)  
  return(newList)}


##################################
loglikadjcatsmooth<-function(par,resp,k,pred, lambda,weightlasso,singleset){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  beta<-par[1:((k-1)*p)]  ### first are covariate weights
  betamat<-matrix(beta,p,k-1)
  beta0<-c(0,par[((k-1)*p+1):((k-1)*p+k-1)])   ### thresholds
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%betamat[,r-1])}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum}
    sum(prob[i,])
    loglik<-loglik+log(prob[i,])[resp[i]]
  }
  
  #### penalty fusion L_2
  pen<-0
  #c<-.01
  if(lambda>0 ){
    d1<- -diag(k-1)
    d2<-cbind(as.matrix(rep(0,k-1)),diag(k-1))
    dif<-d1[1:(k-2),1:(k-1)]+d2[1:(k-2),1:(k-1)]  ### difference generating matrix 
    M<-t(dif)%*%dif
    for (v in singleset){
      #for (r1 in 1:(k-2)){
        #for (r2 in (r1+1):(k-1)) 
        #pen<- pen + (betamat[v,r1]-betamat[v,r1+1])^2
        pen<- pen + t(dif%*%betamat[v,])%*%dif%*%betamat[v,]
                #}
    }
    loglik<-loglik- lambda*pen
  } ### end lambd
  
  
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  #return(newList)
  return(loglikneg)
}

#######################################################



##################################
derloglikadjcatsmooth<-function(par,resp,k,pred, lambda,weightlasso,singleset){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  #p<-dim(pred)[2]
  #beta<-par[1:p]  ### first are covariate weights
  #beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  
  p<-dim(pred)[2]
  beta<-par[1:((k-1)*p)]  ### first are covariate weights
  betamat<-matrix(beta,p,k-1)
  betamatext<-cbind(rep(0,p),betamat)   ### extended
  beta0<-c(0,par[((k-1)*p+1):((k-1)*p+k-1)])   ### thresholds
  
  
  
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,k-1)  ## matrix
  der2<-matrix(0,(k-1),1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%betamat[,r-1])}
    
    ### etasumterm
    etasum<- matrix(0,1,k)
    for (r in 2:k){etasum[1,r]<-sum(eta[1,2:r])}  
    
    ### expterm  
    expterm<- matrix(0,1,k)
    for (r in 1:k){expterm[1,r]<-exp(etasum[1,r])}
    
    #  sum denominator
    sum <-1
    for (r in 2:k){sum<-sum+expterm[1,r]}
    
    sum1<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) sum1<-sum1 +prob[i,r]*(r-1) }
    #sum(prob[i,])
    
    cat<-resp[i]
    
    ### covariates
    #der1<-der1+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der1[,j-1]<-der1[,j-1]+(ind1-sum(prob[i,j:k]))*as.matrix(pred[i,])
    }
    
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
  }  ## end i
  
  #### penalty smooth
  #pen<-0
  #c<-.01
  if(lambda>0){
    
    #for (s in 2:k){
    #  for (j in singleset){
    #    if(s==2)  der1[j,1]<-der1[j,1]-lambda*NormShiftoder((betamat[j,1]-betamat[j,2]),c)
    #    if(s==k)  der1[j,k-1]<-der1[j,k-1]+-lambda*NormShiftoder((betamat[j,k-2]-betamat[j,k-1]),c)
        
     #   if((s >2) & (s <k))  {der1[j,s-1]<-der1[j,s-1]-lambda*NormShiftoder((betamat[j,s-1]-betamat[j,s]),c)
     #   der1[j,s]<-der1[j,s]+lambda*NormShiftoder((betamat[j,s-1]-betamat[j,s]),c)}
        
     # }} ### end k
  #} # enf if
  
  ####
  d1<- -diag(k-1)
  d2<-cbind(as.matrix(rep(0,k-1)),diag(k-1))
  dif<-d1[1:(k-2),1:(k-1)]+d2[1:(k-2),1:(k-1)]  ### difference generating matrix 
  M<-t(dif)%*%dif
  
  for (v in singleset){
  #  for (r1 in 1:(k-2)){
    #    pen<- pen + t(dif%*%betamat[v,])%*%dif%*%betamat[v,]
  der1[v, ]<-der1[v, ] - lambda*2*M%*%as.matrix(betamat[v,])
  }
  } # enf if
  ####
  
  der<-c(c(der1),der2)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

#####
