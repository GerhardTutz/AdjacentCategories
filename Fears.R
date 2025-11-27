

###  adjacent categories, stereotype for fears data

library("VGAM")
library("effectsize")

#load("C:/Users/tutz/LRZ Sync+Share/ABookOrdinal/R/Families/GLES17angst.rda")

load("./GLES17angst.rda")
summary(GLES)

GLES$Age<-standardize(GLES$Age)   ### age standardized  !!!
#GLES$Gender<-standardize(GLES$Gender)
#GLES$Unemployment<-standardize(GLES$Unemployment)
#GLES$EastWest<-standardize(GLES$EastWest)
#GLES$Abitur<-standardize(GLES$Abitur)
dim(GLES)

### subsample, not used
#set.seed(1)
#sam<-sample(dim(GLES)[1], 300, replace = FALSE, prob = NULL)

### all observations:
sam<-1:2036

pred <- as.matrix(GLES[sam,c(7,8,9,10,11)]) #age,gender,eastwest, abitur unemployment
resp<-as.matrix(GLES[sam,2]) ### climate change
#resp<-as.matrix(GLES[sam,6]) ### nuclear energy
#resp<-as.matrix(GLES[sam,4]) ### globalization


summary(GLES[sam,])





#### fit models

#### vglm fit adjacent Climate change!

k<-7
GLESred<-GLES[sam,]
fitvglm <- vglm(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,family=acat(parallel=TRUE),data=GLESred)
fitvglm
summary(fitvglm)

#non parallel
GLESred<-GLES[sam,]
fitvglmext <- vglm(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,family=acat(parallel=FALSE),data=GLESred)
fitvglmext
summary(fitvglmext)


### adjacent from own programs, also category specific
source("ProgramsAdjacent.R")
source("ProgramsAdjacentCatSpec.R")
source("ProgramsStereotype.R")


fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit
fit$Loglik
round(fit$location,digits=3) 

fitcats <- fitadjcats(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=0,lambda=0,single=0)
fitcats 
round(fitcats$parametermatrix,digits=3)

lr<--2*(fit$Loglik-fitcats$Loglik)
1-pchisq(lr, df=25)



### visualization
i<-4 ## variable
ylims<-c(min(fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,]),max(fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[1,]))
plot(seq(2,k),fitcats$parametermatrix[i,],  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2,ylim=ylims)

lines(seq(2,k),fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,],lwd=.5)
lines(seq(2,k),fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,],lwd=.5)
lines(c(2,k),c(0,0),lwd=1.5,type="l",lty=2)
title(ylab=(paste("variable",i)), line=2.32, cex.lab=1.7,lwd=2, cex=1.7)




## cumulative models
library("ordinal")
library("VGAM")
library("xtable")

dat<-GLESred
dat$ClimateChange<-as.ordered(dat$ClimateChange)
fit1 <- vglm(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,
             family=cumulative(parallel=TRUE,reverse=TRUE),
             data=dat)
summary(fit1)
xtable(coef(summary(fit1)), digits=3)
#round(exp(coef(fit1)),3)

fit2 <- vglm(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,
             family=cumulative(parallel=FALSE,reverse=TRUE),
             data=dat)
summary(fit2)
xtable(coef(summary(fit2)), digits=3)
#round(exp(coef(fit2)),3)


# alternatives 
dat$ClimateChange <- as.factor(dat$ClimateChange)
fit3 <- clm(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,  data=dat, link = "logit")
summary(fit3)



########################################
#### stereotype ordinal (monotone)
#############################################

fitstereom<-fitstereomon(resp,k,pred,deriv='der',hessian=FALSE,maxit=2000,start=c(fit$parameter,rep(0,k-2)),lasso=0,c=0.01)
fitstereom

lr<--2*(fit$Loglik-fitstereom$Loglik)
1-pchisq(lr, df=5)

plot(seq(1,k),fitstereom$scalefull[1:k], ylab=" ",xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
title(ylab=expression(paste(gamma[ ]," - scales")), line=2.32, cex.lab=1.7)


plot(seq(1,k),fitstereom$scalesum,  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
lines(c(1,k),c(0,1),lwd=1.0)
title(ylab=expression(paste(phi[ ]," - scales")), line=2.32, cex.lab=1.7)



fitstd<-fitstereom$scalefull/(sum(fitstereom$scalefull))
plot(seq(1,k),fitstd, ylab=" ",xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
title(ylab=expression(paste(gamma[ ]," - scales")), line=2.32, cex.lab=1.7)


### stereotype non-ordinal

fitstereod<-fitstereo(resp,k,pred,deriv='der',hessian=FALSE,maxit=2000,start=c(fit$parameter,rep(1,k-2)),lasso=000,c=0.01)
fitstereod

plot(seq(1,k),fitstereod$scalefull[1:k], ylab=" ",xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
title(ylab=expression(paste(gamma[ ]," - scales")), line=2.32, cex.lab=1.7)


plot(seq(1,k),fitstereod$scalesum,  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
lines(c(1,k),c(0,1),lwd=1.0)
title(ylab=expression(paste(phi[ ]," - scales")), line=2.32, cex.lab=1.7)



#############################
### category specific model
#################################

###loops lamda

#### cross validation  
cross<-5
shortversion<-1  ## if 1 means only one split but with number of observations computed with cross

if(cross>0){
  psect<-floor(dim(pred)[1]/cross) 
  seqsec<-seq(1,dim(pred)[1],psect)  
  seqsec<-c(seqsec,dim(pred)[1])  ##sequence cutpoints (last is n)
  set.seed(10)  ####set seed
  seqsample<-sample(c(1:dim(pred)[1]), size=dim(pred)[1], replace = FALSE, prob = NULL)
  respcross<-resp[seqsample,]
  predcross<-pred[seqsample,]
}

#### lambda values

#lambdav<- c(0,0.1,0.5,20,100,1000)
#lambdav<- c(0.01,.3,.5,1,3,5,10,30,50,100,300,1000,30000) ## used!
lambdav<- c(0.01,.5,1,3,5,10,30,50,100,300,1000,30000)

#lambdav<- c(0,0.1,.5)

p<-dim(pred)[2]
pararr<-array(0,dim =c(length(lambdav),p,k-1))

weightlasso<-.3

loglikv<-matrix(0,length(lambdav),1)
loglikcrossm<-matrix(0,length(lambdav),cross)

for (l in 1:length(lambdav)){
  if(l==1){
    fitl<-fitadjcats(resp,k,pred,deriv='der',hessian=FALSE,maxit=1000,start=0,lambda=lambdav[1],weightlasso,single=0)
    pararr[l,,]<-fitl$parametermatrix
    print(fitl$convergence)
  }
  
  if(l>1){maxit <- 3000
  fitl<-fitadjcats(resp,k,pred,deriv='der',hessian=FALSE,maxit=maxit,start=fitl$parameter,lambda=lambdav[l],weightlasso,single=0)
  pararr[l,,]<-fitl$parametermatrix
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
      respfit<-as.matrix(respcross[-(indsec)])
      predfit<-predcross[-indsec,]
      
      #maxnow<-floor(maxiterations)/4
      dim(respfit)
      dim(predfit)
      # fit
      if(l==1){
        fitc<-fitadjcats(respfit,k,predfit,deriv='der',hessian=FALSE,maxit=1000,start=fitl$parameter,lambda=lambdav[1],weightlasso,single=0)}
      if(l>1){
        fitc<-fitadjcats(respfit,k,predfit,deriv='der',hessian=FALSE,maxit=1000,start=fitc$parameter,lambda=lambdav[l],weightlasso,single=0)}
      
      ## loss 
      
      ###compute loglikhood
      #secfit<-fitadjcats(respcrossnow,k,predcrossnow,deriv='non-der',hessian=FALSE,maxit=1,start=fitc$parameter,lambda=lambdav[l],weightlasso,single=0)
      secfit<--loglikadjcatcats(fitc$parameter,respcrossnow,k,predcrossnow, lambda,weightlasso,single=rep(1:p))
        
      print(secfit)
      loglikcrossm[l,ll]<-secfit
      loglikcross<-loglikcross+secfit/crossnum
    } ## end ll
    
    loglikv[l,1]<-loglikcross
    
    
  }### end cross
  
  } # end loop



maxpl<-max(pararr)
minpl<-min(pararr)

for(v in 1:p){
  var<-v
  maxpl<-max(pararr[,var,])
  minpl<-min(pararr[,var,])
  plot(log(1+lambdav),pararr[,var,1], type="b", ylim=c(minpl,maxpl),ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
  title(ylab=(paste("variable",v)), line=2.32, cex.lab=1.7,lwd=2, cex=1.7)
  for( r in 2:(k-1))lines(log(1+lambdav),pararr[,var,r], type="b",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
}


#### plot cross
plot(log(1+lambdav),-loglikv,type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)
title(xlab=(paste("log(1+lambda)")), cex.lab=1.7,lwd=2, cex=1.7)

ylimn<-c(min(loglikcrossm),max(loglikcrossm))
plot(log(1+lambdav),loglikcrossm[,1],type='p',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2,ylim=ylimn)
for( l in 2:cross){
  lines(log(1+lambdav),loglikcrossm[,l],type='p',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)  
}
lines(log(1+lambdav),loglikv,type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)



### single parameters


coll<-matrix(0,p,4)
for(v in 1:p){
fitt<-fitadjcats(resp,k,pred,deriv='der',hessian=FALSE,maxit=maxit,start=fitcats$parameter,lambda=800,weightlasso=0,single=c(v) )
fitt
lr<--2*(fitt$Loglik-fitcats$Loglik)
print (v)
pv<-1-pchisq(lr, df=5)
print(pv)
coll[v,1]<-v
coll[v,2]<-fitt$Loglik
coll[v,3]<-lr
coll[v,4]<-pv
}

fit$Loglik
fitcats$Loglik
round(coll, 3)


fitt<-fitadjcats(resp,k,pred,deriv='der',hessian=FALSE,maxit=maxit,start=fitcats$parameter,lambda=100,weightlasso=0,single=5 )
fitt
fitt$Loglik


###############################################
# location-shift model
######################################


library("ordDisp")
 

dat<-GLES
summary(dat)
dat$ClimateChange <- as.factor(dat$ClimateChange)
dat$ClimateChange <- as.ordered(dat$ClimateChange)
 
cumshift <- ordDisp(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment|Age+Gender+EastWest+Abitur+Unemployment,
                     data=dat,family="cumulative",reverse=TRUE)
summary(cumshift)
xtable(coef(summary(cumshift)), digits=3)

#pdf("./Dispersion/starplotnuclear.pdf", width=8, height=8)
plot(cumshift,names=c("Age","Gender","EastWest", "Abitur","Unemployment"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)

plot(cumshift,names=c("Age","Gender","EastWest", "Abitur"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)


### acat
adj <- ordDisp(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment,
                    data=dat,family="acat",reverse=FALSE)

summary(adj)

adjshift <- ordDisp(ClimateChange~Age+Gender+EastWest+Abitur+Unemployment|Age+Gender+EastWest+Abitur+Unemployment,
                    data=dat,family="acat",reverse=FALSE)

summary(adjshift)
AIC<- -2*(-3585.3 -16)

plot(adjshift,names=c("Age","Gender","EastWest", "Abitur"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)
plot(adjshift,names=c("Age","Gender","EastWest" ),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,     cex=1.4)

#plot(adjshift,names=c("Abitur","Unemployment"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,    cex=1.4)




###loops

maxit<-5000
lambda<- c(0,1,5,10,20,100)
#lambda<- c(0,10,1000)
lambda<- c(200,100,80,70,50,25,20,15,5,3,1,0.01)

scales<-matrix(0,length(lambda),k)
betas<-matrix(0,length(lambda),dim(pred)[2])


#derivn<-0
derivn<-'der'

c<-0.001
parst<-c(fit$parameter,rep(.5,k-2))
#parst[1:dim(pred)[2]]<-rep(0,dim(pred)[2])
for (l in 1:length(lambda)){
  if(l==1){
  fitstereo3<-fitstereomon(resp,k,pred,deriv=derivn,hessian=FALSE,maxit=maxit,start=parst,lasso=lambda[l],c)
  scales[l,]<-fitstereo3$scalefull
  betas[l,]<-fitstereo3$location
  print(fitstereo3$convergence)
  }
  
  if(l>1){maxit <- 5000
  fitstereo3<-fitstereomon(resp,k,pred,deriv=derivn,hessian=FALSE,maxit=maxit,start=fitstereo3$parameter,lasso=lambda[l],c)
  #fitstereo3<-fitstereo(resp,k,pred,deriv=0,hessian=FALSE,maxit=maxit,start=(parst+fitstereo3$parameter)/2,lasso=lambda[l])
  #fitstereo3<-fitstereo(resp,k,pred,deriv=derivn,hessian=FALSE,maxit=maxit,start=parst,lasso=lambda[l])
  scales[l,]<-fitstereo3$scalefull 
  betas[l,]<-fitstereo3$location
  print(fitstereo3$convergence)
  }
} # end l

fitstereo3

maxpl<-max(scales)
minpl<-min(scales)
plot(log(1+lambda),scales[,2], type="b", ylim=c(minpl,maxpl))
for( r in 3:k)lines(log(1+lambda),scales[,r], type="b")

maxpl<-max(betas)
minpl<-min(betas)
plot(log(1+lambda),betas[,1], type="b", ylim=c(minpl,maxpl))
for( r in 2:dim(pred)[2])lines(log(1+lambda),betas[,r], type="b")


