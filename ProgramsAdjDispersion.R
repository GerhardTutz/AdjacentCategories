


##### programs for ADJACENT dispersion model (location-shift) with lasso


###############################
fitadjlocshift<-function(resp,k,pred,disp,deriv,hessian,maxit,start,lambdl){
  
  ### fits location-shift model
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  ### pred. matrix of predictors
  ### lambdl:  if >0 lasso
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  
  if(sum(disp)==0){
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv !='der')fitopt <- optim(parst, loglikadjcat, gr = NULL,resp=resp,k=k,pred=pred,lambdl=lambdl, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt <- optim(parst, loglikadjcat, gr = derloglikadjcat,resp=resp,k=k,pred=pred, lambdl=lambdl,method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  if(sum(disp!=0)){
    betadisp<-rep(.01,dim(disp)[2])
    parst<-c(beta,beta0,betadisp)
    if(sum(start!=0))parst<-start
    if (deriv !='der')fitopt <- optim(parst, loglikadjcatlocshift, gr = NULL,resp=resp,k=k,pred=pred,disp=disp,lambdl=lambdl, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt <- optim(parst, loglikadjcatlocshift, gr = derloglikadjlocshift,resp=resp,k=k,pred=pred,disp=disp,lambdl=lambdl, method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  pardisp<-0
  if(sum(disp!=0)) pardisp<-fitopt$par[(p+k):(p+k+dim(disp)[2]-1)]
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){hessinv<-solve(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  zval<-fitopt$par[1:p]/stderr
  if(sum(disp!=0)){stderrdisp<-std[(p+k):(p+k+dim(disp)[2]-1)]
  zvaldisp<-pardisp/stderrdisp
  pardisp<-cbind(pardisp,stderrdisp,zvaldisp)}
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "pardisp"= pardisp,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}


##################################################

loglikadjcatlocshift<-function(par,resp,k,pred,disp,lambdl){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # disp is matrix with dispersion covariates
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)+(k/2-r+1)*(disp[i,]%*%betadisp)}
    
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
  
  pen<-0
  c<-.001
  if(lambdl>0){
    for (v in 1:p){
      pen<- pen + (beta[v]^2+c)^(1/2) }
    for (v in 1:pdisp){
      pen<- pen + (betadisp[v]^2+c)^(1/2)  }
    loglik<-loglik- lambdl*pen   ### factor!
  } ###
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  #return(newList)}
  return(loglikneg)}
##################################################


##################################
derloglikadjlocshift<-function(par,resp,k,pred,disp,lambdl){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  
  ###### probabilities
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)   #### !
  der2<-matrix(0,(k-1),1)
  der3<-matrix(0,pdisp,1)   ### disp
  
  predcatspec<-matrix(0,k,dim(disp)[2])
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    
    #for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)}
    for (r in 2:k){predcatspec[r,]<-(k/2-r+1)*(disp[i,])   #### !!1
      eta[1,r]<-beta0[r]+(pred[i,]%*%beta)+(k/2-r+1)*(disp[i,]%*%betadisp)}
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
    sum2<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) {sum1<-sum1 +prob[i,r]*(r-1) 
              sum2<-sum2 +prob[i,r]*as.matrix(colSums(predcatspec[1:r,])) }
    }
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]] 
    cat<-resp[i]
    
    ### covariates
    der1<-der1+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der1<-der1+(cat-1)*(as.matrix(c(pred[i,],predcatspec[cat,])))-sum1*(as.matrix(c(pred[i,],predcatspec[cat,])))  ###!
    if(cat >1)der3<-der3+as.matrix(colSums(predcatspec[1:cat,]))-sum2  ###!
    if(cat ==1)der3<-der3-sum2 
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
  }  ## end i
  
  c<-.001
  if(lambdl>0){
    for (j in 1:p) der1[j,1]<-der1[j,1]-lambdl*(beta[j]^2+c)^(-1/2)*beta[j]
      for (j in 1:pdisp) der3[j,1]<-der3[j,1]-lambdl*(betadisp[j]^2+c)^(-1/2)*betadisp[j]
     } ###
  
  
  
  ##### log-lik
  der<-c(der1,der2,der3)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}
  
  
  #############


###adaptive lasso version
###############################
fitadjlocshiftadapt<-function(resp,k,pred,disp,deriv,hessian,maxit,start,lambdl,adaptind){
  
  
  ### fits location-shift model
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  ### pred. matrix of predictors
  ### lambdl:  if >0 lasso
  ### adaptind:  if 1 adaptive lasso used
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  
  
  if(sum(disp)==0){
    parst<-c(beta,beta0)
    if(sum(start!=0))parst<-start
    if (deriv !='der')fitopt <- optim(parst, loglikadjcat, gr = NULL,resp=resp,k=k,pred=pred,lambdl=lambdl, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt <- optim(parst, loglikadjcat, gr = derloglikadjcat,resp=resp,k=k,pred=pred, lambdl=lambdl,method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  if(sum(disp!=0)){
    betadisp<-rep(.01,dim(disp)[2])
    parst<-c(beta,beta0,betadisp)
    if(sum(start!=0))parst<-start
    
    #### estimates for adaptive
    adaptpar<-rep(1,length(parst))
    if(adaptind==1){
    if (deriv !='der')fitoptad <- optim(parst, loglikadjcatlocshift, gr = NULL,resp=resp,k=k,pred=pred,disp=disp,lambdl=0, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitoptad <- optim(parst, loglikadjcatlocshift, gr = derloglikadjlocshift,resp=resp,k=k,pred=pred,disp=disp,lambdl=0, method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    adaptpar<-fitoptad$par
    } # end adapt
    #####
    
    if (deriv !='der')fitopt <- optim(parst, loglikadjcatlocshiftadapt, gr = NULL,resp=resp,k=k,pred=pred,disp=disp,lambdl=lambdl,adaptpar=adaptpar, method = "Nelder-Mead",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
    if (deriv =='der')fitopt <- optim(parst, loglikadjcatlocshiftadapt, gr = derloglikadjlocshiftadapt,resp=resp,k=k,pred=pred,disp=disp,lambdl=lambdl,adaptpar=adaptpar, method = "BFGS",
                                      lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
    
  }
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  pardisp<-0
  if(sum(disp!=0)) pardisp<-fitopt$par[(p+k):(p+k+dim(disp)[2]-1)]
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  
  if(hessian ==TRUE){hessinv<-solve(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  zval<-fitopt$par[1:p]/stderr
  if(sum(disp!=0)){stderrdisp<-std[(p+k):(p+k+dim(disp)[2]-1)]
  zvaldisp<-pardisp/stderrdisp
  pardisp<-cbind(pardisp,stderrdisp,zvaldisp)}
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "pardisp"= pardisp,    "parameter"=fitopt$par, "convergence"= fitopt$convergence)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}

#############################################################################
loglikadjcatlocshiftadapt<-function(par,resp,k,pred,disp,lambdl,adaptpar){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # disp is matrix with dispersion covariates
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  
  betaadapt<-adaptpar[1:p]  ### first are covariate weights
  betadispadapt<-adaptpar[(p+k):(p+k+pdisp-1)]
  
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)+(k/2-r+1)*(disp[i,]%*%betadisp)}
    
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
  
  pen<-0
  c<-.001
  if(lambdl>0){
    for (v in 1:p){
      pen<- pen + (beta[v]^2+c)^(1/2) /abs(betaadapt[v])}
    for (v in 1:pdisp){
      pen<- pen + (betadisp[v]^2+c)^(1/2) /abs(betadispadapt[v]) }
    loglik<-loglik- lambdl*pen   ### factor!
  } ###
  
  ##### log-lik
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  #return(newList)}
  return(loglikneg)}
##################################################


##################################
derloglikadjlocshiftadapt<-function(par,resp,k,pred,disp,lambdl,adaptpar){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #zgamma<- 1  # not yet needed'  
  
  p<-dim(pred)[2]
  pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  betadisp<-par[(p+k):(p+k+pdisp-1)]
  
  betaadapt<-adaptpar[1:p]  ### first are covariate weights
  betadispadapt<-adaptpar[(p+k):(p+k+pdisp-1)]
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)   #### !
  der2<-matrix(0,(k-1),1)
  der3<-matrix(0,pdisp,1)   ### disp
  
  predcatspec<-matrix(0,k,dim(disp)[2])
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    
    #for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)}
    for (r in 2:k){predcatspec[r,]<-(k/2-r+1)*(disp[i,])   #### !!1
    eta[1,r]<-beta0[r]+(pred[i,]%*%beta)+(k/2-r+1)*(disp[i,]%*%betadisp)}
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
    sum2<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    if(r >1) {sum1<-sum1 +prob[i,r]*(r-1) 
    sum2<-sum2 +prob[i,r]*as.matrix(colSums(predcatspec[1:r,])) }
    }
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]] 
    cat<-resp[i]
    
    ### covariates
    der1<-der1+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    #der1<-der1+(cat-1)*(as.matrix(c(pred[i,],predcatspec[cat,])))-sum1*(as.matrix(c(pred[i,],predcatspec[cat,])))  ###!
    if(cat >1)der3<-der3+as.matrix(colSums(predcatspec[1:cat,]))-sum2  ###!
    if(cat ==1)der3<-der3-sum2 
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
  }  ## end i
  
  c<-.001
  if(lambdl>0){
    for (j in 1:p) der1[j,1]<-der1[j,1]-lambdl*(beta[j]^2+c)^(-1/2)*beta[j]/abs(betaadapt[j])
    for (j in 1:pdisp) der3[j,1]<-der3[j,1]-lambdl*(betadisp[j]^2+c)^(-1/2)*betadisp[j]/abs(betadispadapt[j])
  } ###
  
  
  
  ##### log-lik
  der<-c(der1,der2,der3)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}


#############


