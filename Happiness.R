
### happiness: adjacent categories and stereotype model

 
#getwd() # check 
#setwd("C:\\Users\\Gerhard Tutz\\LRZ Sync+Share\\ABookRasch\\R\\BinaryFits")


happy<-readRDS("relgood2M")
summary(happy)

#### happiness   

##fits basic

source("ProgramsAdjacent.R")

pred <- as.matrix(happy[,c(9,1,4,5,6)]) #age,gender,education,walk, neighbors
resp<-(as.matrix(as.numeric(happy[,10]))) ### happiness


k<-10
### fit adjacent category model

fit<-fitadj(resp,k,pred,disp=0,deriv='der',hessian = TRUE,maxit=500,start=0)
fit
fit$location
round(fit$location,digits=3)

## exploratory: uncertainty model 
disp<-pred
fitdispder<-fitadj(resp,k,pred,disp=disp,deriv='der',hessian = TRUE,maxit=500,start=c(fit$parameter,rep(.1,dim(disp)[2])))  ## with der
fitdispder

### with R packages (check)

library("VGAM")
fitvglmh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=acat(parallel=TRUE),data=relgood2)
summary(fitvglmh)

fitcumh <- vglm(happicat~age+Gender+Education+Walk+Neighbours,family=cumulative(parallel=TRUE),data=relgood2)
fitcumh
summary(fitcumh)

### stereo monoton
source("ProgramsStereotype.R")

fitstereom<-fitstereomon(resp,k,pred,deriv='der',hessian=FALSE,maxit=2000,start=c(fit$parameter,rep(0.2,k-2)),lasso=0,c=0.01)
fitstereom


plot(seq(1,k),fitstereom$scalefull[1:k], ylab=" ",xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
title(ylab=expression(paste(gamma[ ]," - scales")), line=2.32, cex.lab=1.7)


plot(seq(1,k),fitstereom$scalesum,  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
lines(c(1,k),c(0,1),lwd=1.0)
title(ylab=expression(paste(phi[ ]," - scales")), line=2.32, cex.lab=1.7)

lr<--2*(fit$Loglik-fitstereom$Loglik)
1-pchisq(lr, df=8)

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

#### category-specific version

source("ProgramsAdjacentCatSpec.R")

k<-10
fit<-fitadj(resp,k,pred,disp=0,der='der',hessian = TRUE,maxit=500,start=0)
fit


fitcats <- fitadjcats(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=0,lambda=0,single=0)
fitcats 
round(fitcats$parametermatrix,digits=3)

lr<--2*(fit$Loglik-fitcats$Loglik)
1-pchisq(lr, df=45)


### visualization without smoothing
i<-1 ## variable
for (i in 1:dim(pred)[2]){
ylims<-c(min(fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,]),max(fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,]))
plot(seq(2,k),fitcats$parametermatrix[i,],  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2,ylim=ylims)
#lines(seq(2,k),fitcatsm$parametermatrix[i,],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2,lty=2)

lines(seq(2,k),fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,],lwd=.5)
lines(seq(2,k),fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,],lwd=.5)
lines(c(2,k),c(0,0),lwd=1.5,type="l",lty=2)
title(ylab=(paste("variable",i)), line=2.32, cex.lab=1.7,lwd=2, cex=1.7)
}

### visualization with smoothing

lambda<-30 
fitcatsm <- fitadjcatsmooth(resp,k,pred,deriv='der',hessian=TRUE,maxit=500,start=fitcats$parameter,lambda=lambda)




i<-1 ## variable
for (i in 1:dim(pred)[2]){
  ylims<-c(min(fitcats$parametermatrix),max(fitcats$parametermatrix))
  plot(seq(2,k),fitcats$parametermatrix[i,],  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2,ylim=ylims)
  #plot(seq(2,k),fitcats$parametermatrix[i,],  ylab="" ,xlab="Categories",cex.lab=1.6,cex.axis=1.4, type="b",lwd=2)
  lines(seq(2,k),fitcatsm$parametermatrix[i,],type='b',xlab="",ylab="",cex.lab=1.6,cex.axis=1.4,lwd=1.5, cex=1.2,lty=2)
  
  #lines(seq(2,k),fitcats$parametermatrix[i,]-1.96*fitcats$Stderr[i,],lwd=.5)
  #lines(seq(2,k),fitcats$parametermatrix[i,]+1.96*fitcats$Stderr[i,],lwd=.5)
  lines(c(2,k),c(0,0),lwd=1.0,type="l",lty=2)
  title(ylab=(paste("variable",i)), line=2.32, cex.lab=1.7,lwd=2, cex=1.7)
}



