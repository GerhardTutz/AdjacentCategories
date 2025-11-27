

#### adjacent dispersion model for heart data with regularization

library("ordinalForest")
library("effectsize")

data(hearth)
names(hearth)[names(hearth)=="Class"] <- "Cat"
#### with heart
dat<- hearth

summary(dat)
dat[1:10,]

for (i in 1:11)dat[,i]<-as.numeric(dat[,i])
dat$sex<-dat$sex-1
dat$fbs<-dat$fbs-1
dat$exang<-dat$exang-1

dat$age<-standardize(dat$age)
dat$chol<-standardize(dat$chol)
dat$thalach<-standardize(dat$thalach)
dat$chest_pain<-standardize(dat$chest_pain)


pred <- as.matrix((dat[,c(1,2,5,6,9,3)])) #age1,gender2,chol5, fbs6, exang9 , chestpain3
resp<-as.matrix(as.numeric(dat$Cat))

k<-5

######################

source("ProgramsAdjacent.R")

fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=1500,start=0)
fit

###  
source("ProgramsAdjDispersion.R")

disp<-pred

# fits with derivative 

lambdl<-10
fitlocshder<-fitadjlocshift(resp,k,pred,disp=disp,der='der',hessian = TRUE,maxit=2500,start=0,lambdl=lambdl)
fitlocshder
#derloglikadjlocshift(fitlocshder$parameter,resp,k,pred,disp,lambdl=lambdl)

### adapt version
lambdl<-10
fitlocshderadn<-fitadjlocshiftadapt(resp,k,pred,disp=disp,der='der',hessian = TRUE,maxit=2500,start=0,lambdl=lambdl,adaptind = 1)
fitlocshderadn


round(fitlocshder$location,3)
round(fitlocshder$pardisp,3)



##########

###loops 
################################
### works for adaptive and simple lasso
maxit<-500
hessian<-FALSE

#lambdav<- c(0,.5,1,4,6,8,12) ##ok with 5
lambdav<- c(0,.5,1,1.5,2,4,6,8,12,20,50,80,200)

adaptvers<-TRUE  ### uses adapt version if true
adaptind<-0
adaptpar<-rep(1,dim(pred)[2]+dim(disp)[2]+k-1)

if(adaptvers){adaptind<-1
fitoptad <- optim(rep(.01,dim(pred)[2]+dim(disp)[2]+k-1), loglikadjcatlocshift, gr = derloglikadjlocshift,resp=resp,k=k,pred=pred,disp=disp,lambdl=0, method = "BFGS",
                  lower = -Inf, upper = Inf,control = list(maxit=maxit), hessian = hessian)
adaptpar<-fitoptad$par}

p<-dim(pred)[2]
pdisp<-dim(disp)[2]
cov<-matrix(0,length(lambdav),p)
covdisp<-matrix(0,length(lambdav),pdisp)

cross<-10

shortversion<-0  ## if 1 means only one split but with number of observations computed with cross
loglikv<-matrix(0,length(lambdav),1)
loglikcrossm<-matrix(0,length(lambdav),cross)



if(cross>0){
  psect<-floor(dim(pred)[1]/cross) 
  seqsec<-seq(1,dim(pred)[1],psect)  
  seqsec<-c(seqsec,dim(pred)[1])  ##sequence cutpoints (last is n)
  set.seed(1)  ####set seed
  seqsample<-sample(c(1:dim(pred)[1]), size=dim(pred)[1], replace = FALSE, prob = NULL)
  respcross<-resp[seqsample,]
  predcross<-pred[seqsample,]
  dispcross<-disp[seqsample,]
}

for (l in 1:length(lambdav)){
  if(l==1){
    fitl<-fitadjlocshift(resp,k,pred,disp=disp,der='der',hessian = FALSE,maxit=2500,start=0,lambdl=lambdav[l])
    cov[l,]<-fitl$location
    covdisp[l,]<-fitl$pardisp
    print(fitl$convergence)
  }
  
  if(l>1){maxit <- 3000
  #fitl<-fitadjlocshift(resp,k,pred,disp=disp,der='der',hessian = FALSE,maxit=2500,start=fitl$parameter,lambdl=lambdav[l])
  fitl<-fitadjlocshiftadapt(resp,k,pred,disp=disp,der='der',hessian = FALSE,maxit=2500,start=fitl$parameter,lambdl=lambdav[l],adaptind = adaptind)
  
  
  cov[l,]<-fitl$location
  covdisp[l,]<-fitl$pardisp
  print(fitl$convergence)
  }
  
  #### cross validation
  
  if(cross>0){
    loglikcross<-0
    crossnum<-cross
    if(shortversion==1)crossnum<-1
    
    for (ll in 1:crossnum){
      indsec<-seqsec[ll]:seqsec[ll+1]
      respcrossnow<-as.matrix(respcross[indsec])
      predcrossnow<-predcross[indsec,]
      dispcrossnow<-dispcross[indsec,]
      respfit<-as.matrix(respcross[-(indsec)])
      predfit<-predcross[-indsec,]
      dispfit<-dispcross[-indsec,]
      #maxnow<-floor(maxiterations)/4
      dim(respfit)
      dim(predfit)
      # fit
      if(l==1){
        #fitc<-fitadjcats(respfit,k,predfit,deriv='der',hessian=FALSE,maxit=1000,start=fitl$parameter,lambda=lambdav[1],weightlasso,single=0)
        fitc<-fitadjlocshift(respfit,k,predfit,disp=dispfit,der='der',hessian = FALSE,maxit=2500,start=fitl$parameter,lambdl=lambdav[l])
      }
      if(l>1){
        #fitc<-fitadjlocshift(respfit,k,predfit,disp=dispfit,der='der',hessian = FALSE,maxit=2500,start=fitc$parameter,lambdl=lambdav[l])
        fitc<-fitadjlocshiftadapt(respfit,k,predfit,disp=dispfit,der='der',hessian = FALSE,maxit=2500,start=fitc$parameter,lambdl=lambdav[l],adaptind = adaptind)
        
      }
      
      ## loss 
      
      ###compute loglikhood
      #if(l==1)
        secfit<- -loglikadjcatlocshift(fitc$parameter,respcrossnow,k,predcrossnow,dispcrossnow,lambdl=0)  ### 1
      
      #if(l>1)
      #{secfit<- -loglikadjcatlocshiftadapt(fitc$parameter,respcrossnow,k,predcrossnow,dispcrossnow,lambdl=lambdav[l],adaptind = adaptind)}
      
      print(secfit)
      loglikcrossm[l,ll]<-secfit
      loglikcross<-loglikcross+secfit/crossnum
    } ## end ll
    
    loglikv[l,1]<-loglikcross
    
    
  }### end cross
  
} # end l







maxpl<-max(cov)
minpl<-min(cov)
plot(log(1+lambdav),cov[,1], type="b", ylim=c(minpl,maxpl),ylab="",xlab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
#title(ylab=(paste("variable",v)), line=2.32, cex.lab=1.7)
for( r in 2:p)lines(log(1+lambdav),cov[,r], type="b", cex.lab=1.7,lwd=2, cex=1.2)
title(xlab=(paste("log(1+lambda)")), cex.lab=1.7,lwd=2, cex=1.2)
title(ylab=(paste("Location")), cex.lab=1.7,lwd=2, cex=1.2)

maxpl<-max(covdisp)
minpl<-min(covdisp)
plot(log(1+lambdav),covdisp[,1], type="b", ylim=c(minpl,maxpl),ylab="",xlab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
#title(ylab=(paste("variable",v)), line=2.32, cex.lab=1.7)
for( r in 2:p)lines(log(1+lambdav),covdisp[,r], type="b", cex.lab=1.7,lwd=2, cex=1.2)
title(xlab=(paste("log(1+lambda)")), cex.lab=1.7,lwd=2, cex=1.2)
title(ylab=(paste("Dispersion")), cex.lab=1.7,lwd=2, cex=1.2)

#### plot cross
plot(log(1+lambdav),-loglikv,type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.5)
title(xlab=(paste("log(1+lambda)")), cex.lab=1.7,lwd=2, cex=1.2)
title(ylab=(paste("Cross validation")), cex.lab=1.7,lwd=2, cex=1.2)

maxpl<-max(-loglikcrossm)
minpl<-min(-loglikcrossm)
plot(log(1+lambdav),-loglikcrossm[,1],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.5, ylim=c(minpl,maxpl))
for( r in 2:cross)lines(log(1+lambdav),-loglikcrossm[,r], type="b", cex.lab=1.7,lwd=2, cex=1.2)
title(xlab=(paste("log(1+lambda)")), cex.lab=1.7,lwd=2, cex=1.2)
title(ylab=(paste("Cross validation")), cex.lab=1.7,lwd=2, cex=1.2)



### plot separate empty

maxpl<-max(cov)
minpl<-min(cov)
plot(log(1+lambdav),cov[,dim(pred)[2]-empty+1], type="b", ylim=c(minpl,maxpl),ylab="",xlab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
#title(ylab=(paste("variable",v)), line=2.32, cex.lab=1.7)
for( r in (dim(pred)[2]-empty+2):dim(pred)[2])lines(log(1+lambdav),cov[,r], type="b", cex.lab=1.7,lwd=2, cex=1.2)
maxpl<-max(covdisp)
minpl<-min(covdisp)
plot(log(1+lambdav),covdisp[,dim(pred)[2]-empty+1], type="b", ylim=c(minpl,maxpl),ylab="",xlab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
#title(ylab=(paste("variable",v)), line=2.32, cex.lab=1.7)
for( r in (dim(pred)[2]-empty+2):dim(pred)[2])lines(log(1+lambdav),covdisp[,r], type="b", cex.lab=1.7,lwd=2, cex=1.2)


