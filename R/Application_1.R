require(gamlss)
require(gamlss.dist)
require(e1071)

source("R/logit_EPE_gamlss.R")

data1 <- read.table("data_applications_1_2.txt",header=T,dec=",")
str(data1)


y <- data1$rh/100
sd(log(y/(1-y)))
#Descriptive
cbind(round(c(mean=mean(y),median=median(y),sd=sd(y),min=min(y),max=max(y),skew=skewness(y),kurtosi=kurtosis(y)),3))


##logit-PE
fit.logit.PE <- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity"),method=CG(1000),c.crit=0.01, 
                       tau.start = 1,tau.fix = T,mu.start =0.5,sigma.start = 0.7,nu.start = 4)
summary(fit.logit.PE)

##logit-EL
fit.logit.EL <- gamlss(y~1,family = EPE01(mu.link = "identity",sigma.link = "identity",tau.link = "identity"),
                       method=CG(1000),c.crit=0.01,nu.start = 1,nu.fix = T,mu.start =0.544)
summary(fit.logit.EL)


##logit-ENO
fit.logit.ENO <- gamlss(y~1,family = EPE01(mu.link = "identity",sigma.link = "identity",tau.link = "identity"),
                        method=CG(1000),c.crit=0.01,nu.fix = T,mu.start =0.5,sigma.start = 0.7)
summary(fit.logit.ENO)


##logit-L
fit.logit.L <- gamlss(y~1,family = EPE01(mu.link = "identity",sigma.link = "identity"),
                      method=CG(1000),c.crit=0.01,nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T,mu.start =0.5,sigma.start = 0.7)
summary(fit.logit.L)


##logit-ENO
fit.logit.NO <- gamlss(y~1,family = EPE01(mu.link = "identity",sigma.link = "identity"),
                       method=CG(1000),c.crit=0.01,nu.fix = T,tau.start = 1,tau.fix = T,mu.start =0.54,sigma.start = 0.712)
summary(fit.logit.NO)


##logit - EPE
fit.logit.EPE <- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"),
                        method=CG(1000),c.crit=0.01,mu.start =0.5,sigma.start = 0.7)
summary(fit.logit.EPE)


#LR test

#LR.test logit-EPE x logit-PE
LR.test(fit.logit.PE,fit.logit.EPE)
(LR1 <- (2*(logLik(fit.logit.EPE)-logLik(fit.logit.PE)))[1])
1-pchisq(LR1,df=1)

#LR.test logit-EPE x logit-EL
LR.test(fit.logit.EL,fit.logit.EPE)
(LR2 <- (2*(logLik(fit.logit.EPE)-logLik(fit.logit.EL)))[1])
1-pchisq(LR2,df=1)


#LR.test logit-EPE x logit-ENO
LR.test(fit.logit.ENO,fit.logit.EPE)
(LR3 <- (2*(logLik(fit.logit.EPE)-logLik(fit.logit.ENO)))[1])
1-pchisq(LR3,df=1)


#LR.test logit-EPE x logit-L
LR.test(fit.logit.L,fit.logit.EPE)
(LR4 <- (2*(logLik(fit.logit.EPE)-logLik(fit.logit.L)))[1])
1-pchisq(LR4,df=2)

#LR.test logit-EPE x logit-NO
LR.test(fit.logit.NO,fit.logit.EPE)
(LR5 <- (2*(logLik(fit.logit.EPE)-logLik(fit.logit.NO)))[1])
1-pchisq(LR5,df=2)




## GBE
fit.GBE <- gamlss(y~1,family = GB1(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),n.cyc=220,c.crit=0.01)
summary(fit.GBE)

## BE
fit.BE <- gamlss(y~1,family = BE(mu.link = "identity",sigma.link = "identity"),method=CG(1000),c.crit=0.01)
summary(fit.BE)

# Simplex
fit.simplex <- gamlss(y~1,family = SIMPLEX(mu.link = "identity",sigma.link = "identity"),method=CG(1000),c.crit=0.01,mu.start = 0.08)
summary(fit.simplex)




x11()
hist(y,freq=F,xlab = "Relative humidity",ylab = "f(y)",cex.lab=1.35,cex.axis=1.35,main="");box()
curve(dEPE01(x, mu=fit.logit.EPE$mu.coefficients,sigma=fit.logit.EPE$sigma.coefficients,nu=fit.logit.EPE$nu.coefficients,tau=fit.logit.EPE$tau.coefficients), add=T,lwd=4,col="gray15")
curve(dBE(x, mu=fit.BE$mu.coefficients,sigma=fit.BE$sigma.coefficients), add=T,lwd=4,col="gray30",lty=2)
curve(dSIMPLEX(x, mu=fit.simplex$mu.coefficients,sigma=fit.simplex$sigma.coefficients), add=T,lwd=4,col="gray45",lty=4)
curve(dGB1(x, mu=fit.GBE$mu.coefficients,sigma=fit.GBE$sigma.coefficients,nu=fit.GBE$nu.coefficients,tau=fit.GBE$tau.coefficients), add=T,lwd=3,col="gray60")
legend("topleft",c("logit-EPE","BE","Simplex","GBE"),lwd=c(3,3,3,3),lty=c(1,2,4,1),col=c("gray15","gray30","gray45","gray60"),bty='n',cex=1.2)

x11()
plot(ecdf(y),col="gray40",xlab = "Relative humidity",ylab = "F(y)",cex.lab=1.35,cex.axis=1.35,main="")
curve(pEPE01(x, mu=fit.logit.EPE$mu.coefficients,sigma=fit.logit.EPE$sigma.coefficients,nu=fit.logit.EPE$nu.coefficients,tau=fit.logit.EPE$tau.coefficients), add=T,lwd=4,col="gray15")
curve(pBE(x, mu=fit.BE$mu.coefficients,sigma=fit.BE$sigma.coefficients), add=T,lwd=4,col="gray30",lty=2)
curve(pSIMPLEX(x, mu=fit.simplex$mu.coefficients,sigma=fit.simplex$sigma.coefficients), add=T,lwd=4,col="gray45",lty=4)
curve(pGB1(x, mu=fit.GBE$mu.coefficients,sigma=fit.GBE$sigma.coefficients,nu=fit.GBE$nu.coefficients,tau=fit.GBE$tau.coefficients), add=T,lwd=3,col="gray60")
legend("topleft",c("Empirical distribution","logit-EPE","BE","Simplex","GBE"),lwd=c(3,3,3,3,3),lty=c(1,1,2,4,1),col=c("gray40","gray15","gray30","gray45","gray60"),bty='n',cex=1.2)



Res.q1 <- fit.logit.EPE$residuals
Res.q2 <- fit.GBE$residuals
Res.q3 <- fit.BE$residuals
Res.q4 <- fit.simplex$residuals
x11()
par(mfrow=c(2,2))
qqnorm(Res.q1,pch=19,col="gray15",ylim = c(-3.5,3),xlim=c(-3,3),main="logit-EPE",
       ylab="Sample Quantiles",xlab="Theorical Quantiles",cex.lab=1.2)
qqline(Res.q1,col="gray55",lwd=2)
qqnorm(Res.q2,pch=19,col="gray15",ylim = c(-3.5,3),xlim=c(-3,3),
       ylab="Sample Quantiles",xlab="Theorical Quantiles",main = "GBE",cex.lab=1.2)
qqline(Res.q2,col="gray55",lwd=2)
qqnorm(Res.q3,pch=19,col="gray15",ylim = c(-3.5,3),xlim=c(-3,3),
       ylab="Sample Quantiles",xlab="Theorical Quantiles",main = "BE",cex.lab=1.2)
qqline(Res.q3,col="gray55",lwd=2)
qqnorm(Res.q4,pch=19,col="gray15",ylim = c(-3.5,3),xlim=c(-3,3),
       ylab="Sample Quantiles",xlab="Theorical Quantiles",main = "Simplex",cex.lab=1.2)
qqline(Res.q4,col="gray55",lwd=2)



#Cramer Von Mises:
n<-length(y)
x<-1:n # para calcular o CVM


## EPE01
cvm1<-1/(12*n)+sum(((2*x-1)/(2*n)-pEPE01(sort(y),mu=fit.logit.EPE$mu.coefficients,sigma=fit.logit.EPE$sigma.coefficients,nu=fit.logit.EPE$nu.coefficients,tau=fit.logit.EPE$tau.coefficients))^2)
cvm1



## Beta
cvm2<-1/(12*n)+sum(((2*x-1)/(2*n)-pBE(sort(y), mu=fit.BE$mu.coefficients,sigma=fit.BE$sigma.coefficients))^2)
cvm2

## G Beta
cvm22<-1/(12*n)+sum(((2*x-1)/(2*n)-pGB1(sort(y), mu=fit.GBE$mu.coefficients,sigma=fit.GBE$sigma.coefficients,nu=fit.GBE$nu.coefficients,tau=fit.GBE$tau.coefficients))^2)
cvm22


## Simplex
cvm5<-1/(12*n)+sum(((2*x-1)/(2*n)-pSIMPLEX(sort(y), mu=fit.simplex$mu.coefficients,sigma=fit.simplex$sigma.coefficients))^2)
cvm5


## logit-EPE
ks.test(y,"pEPE01",mu=fit.logit.EPE$mu.coefficients,sigma=fit.logit.EPE$sigma.coefficients,nu=fit.logit.EPE$nu.coefficients,tau=fit.logit.EPE$tau.coefficients)

## G Beta
ks.test(y,"pGB1",mu=fit.GBE$mu.coefficients,sigma=fit.GBE$sigma.coefficients,nu=fit.GBE$nu.coefficients,tau=fit.GBE$tau.coefficients)

## Beta
ks.test(y,"pBE",mu=fit.BE$mu.coefficients,sigma=fit.BE$sigma.coefficients)

## Simplex
ks.test(y,"pSIMPLEX",mu=fit.simplex$mu.coefficients,sigma=fit.simplex$sigma.coefficients)

