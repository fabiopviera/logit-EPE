
## Test
require(numDeriv)
require(gamlss)

source("logit_EPE_gamlss.R")


par(mfrow=c(2,2))
curve(dEPE01(x,mu=0.3,sigma=2,nu=0.65,tau=1.2),add=F,xlim=c(0,1),
      lwd=4,ylab="f(y)",xlab="y")
curve(pEPE01(x,mu=0.3,sigma=2,nu=0.65,tau=1.2),add=F,xlim=c(0,1),
      lwd=4,ylab="F(y)",xlab="y")
curve(qEPE01(x,mu=0.3,sigma=2,nu=0.65,tau=1.2),add=F,xlim=c(0,1),
      lwd=4,ylab="Q(y)",xlab="y")
hist(rEPE01(n=5000,mu=0.3,sigma=2,nu=0.65,tau=1.2),40,freq=F,xlim=c(0,1),col="gray65",ylim=c(0,3.5),main="",xlab="y",cex.lab=1.35,cex.axis=1.35);box()



par(mfrow=c(2,2))
curve(dEPE01(x,mu=0.45,sigma=0.65,nu=4,tau=0.8),add=F,xlim=c(0,1),
      lwd=4,ylab="f(y)",xlab="y")
curve(pEPE01(x,mu=0.45,sigma=0.65,nu=4,tau=0.8),add=F,xlim=c(0,1),
      lwd=4,ylab="F(y)",xlab="y")
curve(qEPE01(x,mu=0.45,sigma=0.65,nu=4,tau=0.8),add=F,xlim=c(0,1),
      lwd=4,ylab="Q(y)",xlab="y")
hist(rEPE01(n=5000,mu=0.45,sigma=0.65,nu=4,tau=0.8),40,freq=F,xlim=c(0,1),col="gray65",ylim=c(0,3.5),main="",xlab="y",cex.lab=1.35,cex.axis=1.35);box()
