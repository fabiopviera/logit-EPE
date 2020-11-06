require(gamlss)
require(gamlss.dist)
require(e1071)

source("logit_EPE_gamlss.R")

data1 <- read.table("data_paper_1_2.txt",header=T,dec=",")
str(data1)


y <- data1$sw/100

#Descriptive
cbind(round(c(mean=mean(y),median=median(y),sd=sd(y),min=min(y),max=max(y),skew=skewness(y),kurtosi=kurtosis(y)),3))



##logit-PE
#For mu = logit link, sigma and nu log link and tau = 1
#fit.logit.PE <- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01, tau.start = 1,tau.fix = T)
#For identity link
fit.logit.PE <- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),method=CG(1000),c.crit=0.01,
tau.start = 1,tau.fix = T,control=gamlss.control(tau.step=1))

summary(fit.logit.PE)
LR.test(fit02,fit0)



##logit - EPE

#For mu = logit link, sigma, nu and tau log link
#fit.logit.EPE <- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01) 
#For identity link
fit.logit.EPE <- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"),method=CG(1000),c.crit=0.01,mu.start = 0.08,sigma.start = 0.02)
summary(fit.logit.EPE)

##logit-EL
#For mu = logit link, sigma and tau log link and nu = 1
#fit.logit.EL <- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01,nu.start = 1,nu.fix = T)
fit.logit.EL <- gamlss(y~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"),method=CG(1000),c.crit=0.01,nu.start = 1,nu.fix = T,mu.start =0.0544,sigma.start = 0.2348)
summary(fit.logit.EL)

#LR.test
LR.test(fit.logit.EL,fit.logit.EPE)




## logit-ENO
fit03<- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01,
               nu.start = 2,nu.fix = T)
fit03<- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),method=CG(1000),c.crit=0.01,
nu.start = 2,nu.fix = T)

LR.test(fit03,fit0)


## logit-L
fit04<- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01,
               nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

fit04<- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),method=CG(1000),c.crit=0.01,
nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

LR.test(fit04,fit0)

LR.test(fit04,fit01)


## logit-NO
fit05<- gamlss(y ~1,family = EPE01(),method=CG(1000),c.crit=0.01,
               nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)
fit05<- gamlss(y ~1,family = EPE01(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),method=CG(1000),c.crit=0.01,
nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)

LR.test(fit05,fit0)

LR.test(fit05,fit03)

fit3<- gamlss(y~1,family = GB1(),n.cyc=220,c.crit=0.01)

fit3<- gamlss(y~1,family = GB1(mu.link = "identity",sigma.link = "identity",nu.link = "identity",tau.link = "identity"
),n.cyc=220,c.crit=0.01)
plot(fit3)
summary(fit3)

fit4<- gamlss(y~1,family = BE(),n.cyc=220,c.crit=0.01)

fit4<- gamlss(y~1,family = BE(mu.link = "identity",sigma.link = "identity"),n.cyc=220,c.crit=0.01)
plot(fit4)
summary(fit4)

fit5<- gamlss(y~1,family = SIMPLEX(),n.cyc=220,c.crit=0.01)


fit5<- gamlss(y~1,family = SIMPLEX(mu.link = "identity",sigma.link = "identity"),n.cyc=220,c.crit=0.01,data=dados)
plot(fit5)
summary(fit5)

AIC(fit0,fit3,fit4,fit5)
BIC(fit0,fit3,fit4,fit5)



x11()
hist(y,freq=F,xlab = "Soil water content in the 0-1 m soil layer",ylab = "f(y)",cex.lab=1.35,cex.axis=1.35,main="")
#curve(dEPE01(x, mu=exp(fit0$mu.coefficients)/(1+exp(fit0$mu.coefficients)),sigma=exp(fit0$sigma.coefficients),nu=exp(fit0$nu.coefficients),tau=exp(fit0$tau.coefficients)), add=T,lwd=4,col="gray15")
curve(dEPE01(x, mu=exp(fit01$mu.coefficients)/(1+exp(fit01$mu.coefficients)),sigma=exp(fit01$sigma.coefficients),nu=1,tau=exp(fit01$tau.coefficients)), add=T,lwd=4,col="gray15")
#curve(dEPE01(x, mu=fit01$mu.coefficients,sigma=fit01$sigma.coefficients,nu=1,tau=fit01$tau.coefficients), add=T,lwd=4,col="gray15")

curve(dBE(x, mu=exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients))), add=T,lwd=4,col="gray30",lty=2)
curve(dSIMPLEX(x, mu=exp(fit5$mu.coefficients)/(1+exp(fit5$mu.coefficients)),sigma=exp(fit5$sigma.coefficients)), add=T,lwd=4,col="gray45",lty=4)
curve(dGB1(x, mu=exp(fit3$mu.coefficients)/(1+exp(fit3$mu.coefficients)),sigma=exp(fit3$sigma.coefficients)/(1+exp(fit3$sigma.coefficients)),nu=exp(fit3$nu.coefficients),tau=exp(fit3$tau.coefficients)), add=T,lwd=3,col="gray60")
legend("topright",c("logit-EL","BE","Simplex","GBE"),lwd=c(3,3,3,3),lty=c(1,2,4,1),col=c("gray15","gray30","gray45","gray60"),bty='n',cex=1.2)
box()
curve(dEPE01(x, mu=exp(fit01$mu.coefficients)/(1+exp(fit01$mu.coefficients)),sigma=exp(fit01$sigma.coefficients),nu=1,tau=exp(fit01$tau.coefficients)), add=T,lwd=4,col="gray15")

x11()
plot(ecdf(dados$sw/100),col="gray40",xlab = "Soil water content in the 0-1 m soil layer",ylab = "F(y)",cex.lab=1.35,cex.axis=1.35,main="")
#curve(pEPE01(x, mu=exp(fit0$mu.coefficients)/(1+exp(fit0$mu.coefficients)),sigma=exp(fit0$sigma.coefficients),nu=exp(fit0$nu.coefficients),tau=exp(fit0$tau.coefficients)), add=T,lwd=4,col="gray15")
curve(pEPE01(x, mu=exp(fit01$mu.coefficients)/(1+exp(fit01$mu.coefficients)),sigma=exp(fit01$sigma.coefficients),nu=1,tau=exp(fit01$tau.coefficients)), add=T,lwd=4,col="gray15")

curve(pBE(x, mu=exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients))), add=T,lwd=4,col="gray30",lty=2)
curve(pSIMPLEX(x, mu=exp(fit5$mu.coefficients)/(1+exp(fit5$mu.coefficients)),sigma=exp(fit5$sigma.coefficients)), add=T,lwd=4,col="gray45",lty=4)
curve(pGB1(x, mu=exp(fit3$mu.coefficients)/(1+exp(fit3$mu.coefficients)),sigma=exp(fit3$sigma.coefficients)/(1+exp(fit3$sigma.coefficients)),nu=exp(fit3$nu.coefficients),tau=exp(fit3$tau.coefficients)), add=T,lwd=3,col="gray60")
legend("topleft",c("Empirical distribution","logit-EL","BE","Simplex","GBE"),lwd=c(3,3,3,3,3),lty=c(1,1,2,4,1),col=c("gray40","gray15","gray30","gray45","gray60"),bty='n',cex=1.2)


Res.q1 <- fit.logit.EL$residuals
Res.q2 <- fit3$residuals
Res.q3 <- fit4$residuals
Res.q4 <- fit5$residuals
x11()
par(mfrow=c(2,2))
qqnorm(Res.q1,pch=19,col="gray15",ylim = c(-3.5,3),xlim=c(-3,3),main="logit-EL",
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
n<-length(dados$sw/100)
x<-1:n # para calcular o CVM


## EPE01
cvm1<-1/(12*n)+sum(((2*x-1)/(2*n)-pEPE01(sort(dados$sw/100), mu=exp(fit01$mu.coefficients)/(1+exp(fit01$mu.coefficients)),sigma=exp(fit01$sigma.coefficients),nu=1,tau=exp(fit01$tau.coefficients)))^2)
cvm1

cvm1<-1/(12*n)+sum(((2*x-1)/(2*n)-pEPE01(sort(dados$sw/100), mu=fit0$mu.coefficients,sigma=fit0$sigma.coefficients,nu=fit0$nu.coefficients,tau=fit0$tau.coefficients))^2)
cvm1

## Beta
cvm2<-1/(12*n)+sum(((2*x-1)/(2*n)-pBE(sort(dados$sw/100), mu=exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients))))^2)
cvm2

## G Beta
cvm22<-1/(12*n)+sum(((2*x-1)/(2*n)-pGB1(sort(dados$sw/100), mu=exp(fit3$mu.coefficients)/(1+exp(fit3$mu.coefficients)),sigma=exp(fit3$sigma.coefficients)/(1+exp(fit3$sigma.coefficients)),nu=exp(fit3$nu.coefficients),tau=exp(fit3$tau.coefficients)))^2)
cvm22


## Simplex
cvm5<-1/(12*n)+sum(((2*x-1)/(2*n)-pSIMPLEX(sort(dados$sw/100), mu=exp(fit5$mu.coefficients)/(1+exp(fit5$mu.coefficients)),sigma=exp(fit5$sigma.coefficients)))^2)
cvm5


## logit-EPE
ks.test(dados$sw/100,"pEPE01",mu=exp(fit01$mu.coefficients)/(1+exp(fit01$mu.coefficients)),sigma=exp(fit01$sigma.coefficients),nu=1,tau=exp(fit01$tau.coefficients))

## G Beta
ks.test(dados$sw/100,"pGB1",mu=exp(fit3$mu.coefficients)/(1+exp(fit3$mu.coefficients)),sigma=exp(fit3$sigma.coefficients)/(1+exp(fit3$sigma.coefficients)),nu=exp(fit3$nu.coefficients),tau=exp(fit3$tau.coefficients))

## Beta
ks.test(dados$sw/100,"pBE",exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients)))

## Simplex
ks.test(dados$sw/100,"pSIMPLEX",mu=exp(fit5$mu.coefficients)/(1+exp(fit5$mu.coefficients)),sigma=exp(fit5$sigma.coefficients))


