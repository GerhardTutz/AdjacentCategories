

#### stereotype programs

###########fitstereo unordered
###########fitstereomon ordered

###############################
fitstereo<-function(resp,k,pred,deriv,hessian,maxit,start,lasso,c){
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  ##### lasso: penalty for lasso selection, used if greater 0
  ##### c: small value (.01) for lasso estimates
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  scale<-rep(1,k-2)
  
  
  parst<-c(beta,beta0,scale)
  if(sum(start!=0))parst<-start
  if (deriv !='der')fitopt <- optim(parst, loglikstereo, gr = NULL,resp=resp,k=k,pred=pred,lasso=lasso, c=c, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv =='der')fitopt <- optim(parst, loglikstereo, gr = derlogstereo,resp=resp,k=k,pred=pred,lasso=lasso, c=c, method = "BFGS",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  
  
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  scale<-fitopt$par[(p+k):(p+k+k-3)]
  scalefull<-c(0,scale,1)
  
  scalesum<-0
  for (r in 2:k)scalesum<- c(scalesum,sum(scalefull[1:r]))
  scalesum<-scalesum/scalesum[k] 
  
  locm<- matrix(0,length(beta))
  for (r in 2:k)locm<-cbind(locm,as.matrix(scalefull[r]*location))
  
  if(hessian ==TRUE){hessinv<-solve(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  zval<-fitopt$par[1:p]/stderr
  
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "scalefull"=scalefull,  "parameter"=fitopt$par, "convergence"= fitopt$convergence, "scalesum"=scalesum,"locationmatrix"=locm)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}


##################################
loglikstereo<-function(par,resp,k,pred,lasso,c){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # pred is matrix with predictors
  
  #c<-.01  ## lasso  
  
  p<-dim(pred)[2]
  #pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  scale<-c(0,par[(p+k):(p+k+k-3)],1)
  #betadisp<-par[(p+k):(p+k+pdisp-1)]
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)*scale[r]}
    
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
    loglik<-loglik+log(prob[i,resp[i]]) 
  }
  
  ##### log-lik
  
  ### penalty scales monotone
  #pen1<-0
  #for (r in 2:(k-1)) pen1<-pen1+100*Indf(scale[r])
  ##oglik<-loglik-pen1
  
  ##############
  
  ### penalty lasso
  if (lasso >0){
    pen2<-0
    for (p in 1:dim(pred)[2]) {for (r in 2:(k)){
      pen2<-pen2+NormShifto(beta[p]*scale[r],c)}  ### weighted effect
    }
    loglik<-loglik-lasso*pen2}
  
  ##############
  
  
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}
###############


Indf<-function(x){
  ret<-0
  if(x <0)ret<--x
  return(ret)}
################




##################################

derlogstereo<-function(par,resp,k,pred,lasso,c){  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #c<-.01  ## lasso 
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  scale<-c(0,par[(p+k):(p+k+k-3)],1)
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  #der1d<-matrix(0,p,1)
  der2<-matrix(0,(k-1),1)
  der3<-matrix(0,(k-2),1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)*scale[r]}
    
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
    #sum1d<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    #if(r >1) sum1d<-sum1d +prob[i,r]*(r-1) 
    #if(r >1) sum1<-sum1 +prob[i,r]*sum(scale[2:r])
    sum1<-sum1 +prob[i,r]*sum(scale[1:r])
    }
    
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]]
    
    cat<-resp[i]
    
    ### covariates
    
    #der1d<-der1d+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    der1<-der1+sum(scale[1:cat])*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
    
    ### scores
    for (j in 2:(k-1)){
      ind1<-0
      if(cat>=j)ind1<-1
      der3[j-1]<-der3[j-1]+(ind1-sum(prob[i,j:k]))*pred[i,]%*%beta
    } 
    
    
  }  ## end i
  
  ### covariates penalty lasso
  if (lasso >0){
    pen2<-0 
    for (p in 1:dim(pred)[2]) {for (r in 2:(k)){
      der1[p]<-der1[p]-lasso*NormShiftoder(beta[p]*scale[r],c)*scale[r]}  ### weighted effect
    }}
  ### scores penalty lasso
  if (lasso >0){
    for (p in 1:dim(pred)[2]) {for (j in 2:(k-1)){
      der3[j-1]<-der3[j-1]-lasso*NormShiftoder(beta[p]*scale[j],c)*beta[p]}  ### weighted effect
    }}
  
  ##### log-lik
  der<-c(der1,der2,der3)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

###### 



#####################Oelker variante!
NormShifto <-function(x,c){
  ## used argument x^Tx
  y  <-sqrt(x^2+c)-sqrt(c)
  return(y)}
##############

c<-.01 
xpl<-seq(-1,1,.01)
ypl<-NormShifto(xpl[1],c)
for(i in 2:length(xpl))ypl<-c(ypl,NormShifto(xpl[i],c))
plot(xpl,ypl)


NormShiftoder <-function(x,c){
  ## used argument x^Tx
  y  <-(x^2+c)^(-1/2)*x
  return(y)}
##############

c<-.01 
xpl<-seq(-1,1,.01)
ypl<-NormShiftoder(xpl[1],c)
for(i in 2:length(xpl))ypl<-c(ypl,NormShiftoder(xpl[i],c))
plot(xpl,ypl)



######  monotone versions

###############################
fitstereomon<-function(resp,k,pred,deriv,hessian,maxit,start,lasso,c){
  
  ### resp:responses as matrix
  ### pred: predictors as matrix
  ### k: number of categories, 1,2,...,k
  ### der: if 'der' derivatives used
  ### Hessian: TRUE or FALSE
  ### maxit: number of iterations
  ##### lasso: penalty for lasso selection, used if greater 0
  ##### c: small value (.01) for lasso estimates
  
  
  p<-dim(pred)[2]
  beta<-rep(.01,p)  ### first are covariate weights
  beta0<-rep(.01,(k-1)) 
  scale<-rep(1,k-2)   ### mod!
  
  
  parst<-c(beta,beta0,scale)
  if(sum(start!=0))parst<-start
  if (deriv !='der')fitopt <- optim(parst, loglikstereomon, gr = NULL,resp=resp,k=k,pred=pred,lasso=lasso, c=c, method = "Nelder-Mead",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  if (deriv =='der')fitopt <- optim(parst, loglikstereomon, gr = derlogstereomon,resp=resp,k=k,pred=pred,lasso=lasso, c=c, method = "BFGS",
                                    lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
  
  
  
  
  AIC<- -2*(-fitopt$value-length(parst))
  
  
  
  stderr<-0
  stddisp<-0
  location<-fitopt$par[1:p]
  scale<-exp(fitopt$par[(p+k):(p+k+k-3)])  ### exp!
  scalefull<-c(0,scale,1)
  
  scalesum<-0
  for (r in 2:k)scalesum<- c(scalesum,sum(scalefull[1:r]))
  scalesum<-scalesum/scalesum[k] 
  
  locm<- matrix(0,length(beta))
  for (r in 2:k)locm<-cbind(locm,as.matrix(scalefull[r]*location))
  
  if(hessian ==TRUE){hessinv<-solve(fitopt$hessian)
  std<-sqrt(diag(hessinv))  
  stderr<-std[1:p]
  zval<-fitopt$par[1:p]/stderr
  
  
  #stddisp <- std[(p+k):(p+k+dim(disp)[2]-1)]
  location<-cbind( fitopt$par[1:p], stderr,zval)
  }
  
  
  
  #derloglikadjcat(fitopt$par,resp,k,pred)   ## compute derivatives
  
  newList <- list("Loglik"= -fitopt$value,"AIC"=AIC, "location"=location,  
                  "scalefull"=scalefull,  "parameter"=fitopt$par, "convergence"= fitopt$convergence, "scalesum"=scalesum,"locationmatrix"=locm)  # "pare"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}





###################

loglikstereomon<-function(par,resp,k,pred,lasso,c){
  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  # pred is matrix with predictors
  
  #c<-.01  ## lasso  
  
  p<-dim(pred)[2]
  #pdisp<-dim(disp)[2]
  
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  scale<-c(0,exp(par[(p+k):(p+k+k-3)]),1)    ### exp!
  #betadisp<-par[(p+k):(p+k+pdisp-1)]
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  
  
  loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)*scale[r]}
    
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
    loglik<-loglik+log(prob[i,resp[i]]) 
  }
  
  ##### log-lik
  
  ### penalty scales monotone
  #pen1<-0
  #for (r in 2:(k-1)) pen1<-pen1+100*Indf(scale[r])
  ##oglik<-loglik-pen1
  
  ##############
  
  ### penalty lasso
  if (lasso >0){
    pen2<-0
    for (p in 1:dim(pred)[2]) {for (r in 2:(k)){
      pen2<-pen2+NormShifto(beta[p]*scale[r],c)}  ### weighted effect
    }
    loglik<-loglik-lasso*pen2}
  
  ##############
  
  
  
  loglikneg<- -loglik
  
  newList <- list("Loglik"= loglikneg)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(newList)}
###############

##################################

derlogstereomon<-function(par,resp,k,pred,lasso,c){  
  ### negative loglikelihood  
  # par is betavector,beta_02,... 
  
  #c<-.01  ## lasso 
  
  p<-dim(pred)[2]
  beta<-par[1:p]  ### first are covariate weights
  beta0<-c(0,par[(p+1):(p+k-1)])   ### thresholds
  scale<-c(0,exp(par[(p+k):(p+k+k-3)]),1)  ### exp
  ###### probabilities
  
  n<-dim(resp)[1]
  
  prob<- matrix(0,n,k)
  
  der1<-matrix(0,p,1)
  #der1d<-matrix(0,p,1)
  der2<-matrix(0,(k-1),1)
  der3<-matrix(0,(k-2),1)
  
  #loglik<-0
  for(i in 1:n){
    
    ### etaterm
    eta<- matrix(0,1,k) 
    for (r in 2:k){eta[1,r]<-beta0[r]+(pred[i,]%*%beta)*scale[r]}
    
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
    #sum1d<-0
    for (r in 1:k){prob[i,r]<-expterm[1,r]/sum
    #if(r >1) sum1d<-sum1d +prob[i,r]*(r-1) 
    #if(r >1) sum1<-sum1 +prob[i,r]*sum(scale[2:r])
    sum1<-sum1 +prob[i,r]*sum(scale[1:r])
    }
    
    sum(prob[i,])
    #loglik<-loglik+log(prob[i,])[resp[i]]
    
    cat<-resp[i]
    
    ### covariates
    
    #der1d<-der1d+(cat-1)*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    der1<-der1+sum(scale[1:cat])*(as.matrix(pred[i,]))-sum1*(as.matrix(pred[i,]))
    
    
    ### thresholds
    
    for (j in 2:k){
      ind1<-0
      if(cat>=j)ind1<-1
      der2[j-1]<-der2[j-1]+ind1-sum(prob[i,j:k])
    }
    
    ### scores
    for (j in 2:(k-1)){
      ind1<-0
      if(cat>=j)ind1<-1
      der3[j-1]<-der3[j-1]+(ind1-sum(prob[i,j:k]))*pred[i,]%*%beta*scale[j]  ### corr exp!
    } 
    
    
  }  ## end i
  
  ### covariates penalty lasso
  if (lasso >0){
    pen2<-0 
    for (p in 1:dim(pred)[2]) {for (r in 2:(k)){
      der1[p]<-der1[p]-lasso*NormShiftoder(beta[p]*scale[r],c)*scale[r]}  ### weighted effect
    }}
  ### scores penalty lasso
  if (lasso >0){
    for (p in 1:dim(pred)[2]) {for (j in 2:(k-1)){
      der3[j-1]<-der3[j-1]-lasso*NormShiftoder(beta[p]*scale[j],c)*beta[p]}  ### weighted effect
    }}
  
  ##### log-lik
  der<-c(der1,der2,der3)
  
  der<--der
  
  
  #newList <- list("derivative"= der)  # "par"=fits$par,"conv"=fits$convergence,"hessian"=fits$hessian,"sum"=parm)
  return(der)}

###### 


