

######## heart data

library("ordinalForest")
library("effectsize")

data(hearth)
names(hearth)[names(hearth)=="Class"] <- "Cat"


table(hearth$Cat)
dim(hearth)
head(hearth)
summary(hearth)
dat<- hearth


###############################################
###### exploratory: cumulative models

############### fits  package ordinal

library("ordinal")

formula <- Cat ~   age+sex+trestbps+chol+fbs+restecg+thalach+exang+oldpeak +chest_pain

#formula <- Cat ~   sex+chol+exang+oldpeak +chest_pain

dat$Cat<- as.ordered(dat$Cat)
fitcumordinal <- clm(formula,  data=dat, doFit = TRUE, model = TRUE,  link = "logit")
summary(fitcumordinal)

############ vgam: vglm
library("VGAM")

fitp <- vglm(formula,family=cumulative(parallel=TRUE,reverse=TRUE),data=dat)
summary(fitp)


fitseq <- vglm(formula, family=sratio(parallel=TRUE),data=dat)
summary(fitseq)


fitacat <- vglm(formula,family=acat(parallel=TRUE),data=dat)
summary(fitacat)


####additive
#formula <- Cat ~   sex+chol+exang+oldpeak +chest_pain
fitsm <- vgam(Cat ~   sex+s(chol)+fbs+restecg+exang+s(oldpeak) +chest_pain,
              family=cumulative(parallel=TRUE,reverse=TRUE),
              data=dat)
summary(fitsm)
plot(fitsm,se=TRUE)



#######################################
##### adjacent categories models,  covariates standardized
################################

### prepare data
dat<- hearth
summary(dat)


for (i in 1:11)dat[,i]<-as.numeric(dat[,i])
dat$sex<-dat$sex-1
dat$fbs<-dat$fbs-1
dat$exang<-dat$exang-1

dat$age<-standardize(dat$age)
dat$chol<-standardize(dat$chol)
dat$thalach<-standardize(dat$thalach)
dat$chest_pain<-standardize(dat$chest_pain)

#### select response (resp) and predictors (pred)

#pred <- as.matrix((dat[,c(1,2,5,6,8,9,3)])) #age1,gender2,chol5, fbs6,thalach8, exang9 , chestpain3
###without thalach:

pred <- as.matrix((dat[,c(1,2,5,6,9,3)])) #age1,gender2,chol5, fbs6, exang9 , chestpain3
resp<-as.matrix(as.numeric(dat$Cat))

#### fits
source("ProgramsAdjacent.R")
source("ProgramsAdjacentCatSpec.R")

k<-5
fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit


fitcats <- fitadjcats(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=0,lambda=0,single=0)
fitcats 
round(fitcats$parametermatrix, 3)

lr<--2*(fit$Loglik-fitcats$Loglik)
1-pchisq(lr, df=18)

### visualization

i<-6 ## variable
ylims<-c(min(fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,]),max(fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,]))
plot(seq(2,k),fitcats$parametermatrix[i,],  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2,ylim=ylims)

lines(seq(2,k),fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,],lwd=.5)
lines(seq(2,k),fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,],lwd=.5)
lines(c(2,k),c(0,0),lwd=1.5,type="l",lty=2)
title(ylab=(paste("variable",i)), line=2.32, cex.lab=1.7,lwd=2, cex=1.7)


### single parameters

p<-dim(pred)[2]
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




### check: fit with vglm

fitacat <- vglm(Cat~age+sex+chol+fbs+exang+chest_pain,
             family=acat(parallel=TRUE,reverse=FALSE),data=dat)
summary(fitacat)
xtable(coef(summary(fitacat)), digits=3)

### log-lik not available for category-specific

fit2 <- vglm(Cat~age+sex+chol+fbs+exang+chest_pain,
             family=cumulative(parallel=FALSE,reverse=TRUE),data=dat)
summary(fit2)




############################################
#### stereotype monotone
###################################

source("ProgramsStereotype.R")

fitstereom<-fitstereomon(resp,k,pred,deriv='der',hessian=FALSE,maxit=2000,start=c(fit$parameter,rep(0,k-2)),lasso=0,c=0.01)
fitstereom

lr<--2*(fit$Loglik-fitstereom$Loglik)
1-pchisq(lr, df=5)

plot(seq(1,k),fitstereom$scalefull[1:k], ylab=" ",xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
title(ylab=expression(paste(gamma[ ]," - scales")), line=2.32, cex.lab=1.7)


plot(seq(1,k),fitstereom$scalesum,  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
lines(c(1,k),c(0,1),lwd=1.0)
title(ylab=expression(paste(phi[ ]," - scales")), line=2.32, cex.lab=1.7)


### barplot
low<-2
up<-8
plot(fitstereom$scalesum, rep(up,k), ylab=" ",xlab="Scaling",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2,ylim=c(0,10), cex=1.7)
#lines(seq(0,1,1/(k-1)),rep(.2,k),lwd=2, type="b", cex=1.7 )
pcdum <-c ('1','2','3','4','5','6','7','8','9','k')
lines(seq(0,1,1/(k-1)),rep(low,k),lwd=2, type="b",pch =pcdum, cex=1.7)

for (i in 1:k) {
  lines(c((i-1)/(k-1),fitstereom$scalesum[i]),c(low,up),lwd=1,  cex=1.0 )
}

text(x=.5, y= 9.5,label=expression('Fitted scale'   ),lty=2, cex=1.3)
text(x=.5, y= 0.5,label=expression('Equidistance scale'   ),lty=2, cex=1.3)

### for reduced data 
#pred <- as.matrix((dat[,c(1,2,5)])) #age1,gender2,chol5, fbs6, exang9 , chestpain3
#resp<-as.matrix(as.numeric(dat$Cat))
fitstereom<-fitstereomon(resp,k,pred,deriv='der',hessian=FALSE,maxit=2000,start=0,lasso=0,c=0.01)
fitstereom

##############################

# dispersion model : location-shift
library("ordDisp")


summary(dat)
dat$Cat <- as.ordered(dat$Cat)
#age1,gender2,chol5, fbs6,thalach8, exang9 , chestpain3

### acat

adj <- ordDisp(Cat~age+sex+chol+fbs+exang+chest_pain,
               data=dat,family="acat",reverse=FALSE)
summary(adj)

dat$pain<-dat$chest_pain
adjshift <- ordDisp(Cat~age+sex+chol+fbs+exang+chest_pain|age+sex+chol+fbs+exang+chest_pain,
                    data=dat,family="acat",reverse=FALSE)

summary(adjshift)
xtable(coef(summary(adjshift)), digits=3)
coefficients(adjshift)
 
#AIC<- -2*(-3585.3 -16)
#AIC<- -2*(-253.8934 -16)
#lr<--2*(-253.8934-fitcats$Loglik)
#lr<--2*(fit$Loglik-fitcats$Loglik)
#1-pchisq(lr, df=12)

plot(adjshift,names=c("age","sex","chol", "fbs","chest_pain"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)


#### cumulative dispersion

cumshift <- ordDisp(Cat~age+sex+chol+fbs+exang+chest_pain|age+sex+chol+fbs+thalach+exang+chest_pain,
                    data=dat,family="cumulative",reverse=TRUE)
summary(cumshift)
xtable(coef(summary(cumshift)), digits=3)

#pdf("./Dispersion/starplotnuclear.pdf", width=8, height=8)
plot(cumshift,names=c("age","sex","chol", "fbs","thalach"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)
## interessing
plot(cumshift,names=c("age","chol","thalach"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)

plot(cumshift,names=c("fbs"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)
plot(cumshift,names=c("sex"),colorvec=rep(1,5),KI = TRUE, KIfactor = 11/12,
     cex=1.4)



####################################
### category specific with regularization: selection

###loops lamda

#### cross validation  
cross<-5
shortversion<-0  ## if 1 means only one split but with number of observations computed with cross

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

#lambdav<- c(0,1,5)
#lambdav<- c(0.01,.3,.5,1,3,5,10,30,50,100,300,1000,30000) ## used!
lambdav<- c(0.0,.5,1,3,5,8,10,15,20,30,50,100,300,1000)


p<-dim(pred)[2]
pararr<-array(0,dim =c(length(lambdav),p,k-1))

weightlasso<-1

loglikv<-matrix(0,length(lambdav),1)
loglikcrossm<-matrix(0,length(lambdav),cross)



###  loop
set.seed(1)
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
 
    ## loop
    
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
      secfit<--loglikadjcatcats(fitc$parameter,respcrossnow,k,predcrossnow, lambda=lambdav[l],weightlasso,single=rep(1:p))
      
      print(secfit)
      loglikcrossm[l,ll]<-secfit
      loglikcross<-loglikcross+secfit/crossnum
    } ## end ll
    
    loglikv[l,1]<-loglikcross
    
    
  }### end cross
  
} # end loop



#### plots

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

ylimn<-c(min(-loglikcrossm),max(-loglikcrossm))
plot(log(1+lambdav),-loglikcrossm[,1],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2,ylim=ylimn)
for( l in 2:cross){
  lines(log(1+lambdav),-loglikcrossm[,l],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)  
}
lines(log(1+lambdav),loglikv,type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)




####exploratory: additive models

#smoothing variables must have at least 7 unique values!

fitadd <- vgam(Cat ~   age+sex+s(chol)+fbs+exang+chest_pain,
              family=acat(parallel=TRUE,reverse=FALSE),
              data=dat)
summary(fitadd)
plot(fitadd,se=TRUE)
 

fitaddc <- vgam(Cat ~   age+sex+s(chol)+fbs+exang+chest_pain,
               family=acat(parallel=FALSE,reverse=FALSE),
               data=dat)
summary(fitaddc)
plot(fitaddc,se=TRUE)


#### additive alternative
numknotsstart <-4
ord<-4
dep=c("Cat")
lin =c("sex","fbs","restecg","exang","chest_pain")
basis=c("chol", "oldpeak")
k<-5
data=dat

fitad<-FitcumAdd(dep,k,lin, basis,data,numknotsstart,ord,lambd=0.100)

PlotAdditivedata(fitad,lin,basis,data)

PlotAdditiveBoost(fitad,dep,k,lin,basis,data,numknotsstart,ord,lambd=0.01,numsim=2,nquant <-10,curves=1,
                  ylims <- c(-2.8,2.8),alpha=0.10)

##### smooth

k<-5

fitcats <- fitadjcats(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=0,lambda=0,single=0)
fitcats$parametermatrix 


lambda<-10 
fitcatsm <- fitadjcatsmooth(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=0,lambda=lambda)
fitcatsm$parametermatrix
#derloglikadjcatsmooth(fitcatsm$parameter,resp,k,pred, lambda=lambda,weightlasso,singleset)

v<-1
plot(seq(1,(k-1)),fitcats$parametermatrix[v,],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2,lty=2)
lines(seq(1,(k-1)),fitcatsm$parametermatrix[v,],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2)

