require(betareg)
require(gamlss)
require(gamlss.dist)

dados <- read.table("umidade_alininha.txt",header=T,dec=",")
str(dados)
dados$trat <- as.factor(dados$trat)
dados$prof  <- as.factor(dados$prof )
dim(dados)
a=dados$umid[dados$trat=="t1"]
b=dados$umid[dados$trat=="t2"]
c=dados$umid[dados$trat=="t3"]
d=dados$umid[dados$trat=="t4"]
e=dados$umid[dados$trat=="t5"]
f=dados$umid[dados$trat=="t6"]
g=dados$umid[dados$trat=="Cerrado"]
H=list(g,a,b,c,d,e,f)
# Change the names of the elements of the list :
names(H)=c(paste("Control"),paste("T1"), paste("T2"), paste("T3"), paste("T4"), paste("T5"), paste("T6"))



x11()
boxplot(H,ylab="Soil water content (cm³/cm³)",xlab="Treatments",
        col=c("gray15","gray25","gray40","gray50","gray65","gray75","gray95"),
        cex.axis=1.2,cex.lab=1.2)
points(c(1,2,3,4,5,6,7),c(tapply(dados$umid, dados$trat, mean)),pch=19,lwd=4,col="black")
legend("topright",c("Mean"),pch=19,lty=0,col=c("black"),bty='n',cex=1.2)
abline(h=0.284,lwd=3,lty=3)

a=dados$umid[dados$prof =="0-10"]
b=dados$umid[dados$prof =="10-20"]
c=dados$umid[dados$prof =="20-30"]
d=dados$umid[dados$prof =="30-40"]
e=dados$umid[dados$prof =="40-50"]

H=list(a,b,c,d,e)
# Change the names of the elements of the list :
names(H)=c(paste("0-10"),paste("10-20"), paste("20-30"), paste("30-40"), paste("40-60"))



x11()
boxplot(H,ylab="Soil water content (cm³/cm³)",xlab="Depth (cm)",
        col=c("gray15","gray25","gray40","gray50","gray65"),
        cex.axis=1.2,cex.lab=1.2)
points(c(1,2,3,4,5),c(tapply(dados$umid, dados$prof, mean)),pch=19,lwd=4,col="black")
#abline(h=0.2837619,lwd=3,lty=3)
legend("topright",c("Mean"),pch=19,lty=0,col=c("black"),bty='n',cex=1.2)




fit0<- gamlss(umid ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit0)
summary(fit0)



fit3<- gamlss(umid~1,family = LOGITNO(),n.cyc=300,c.crit=0.01,
              data=dados)
plot(fit3)
summary(fit3)


fit4<- gamlss(umid~1,family = BE,n.cyc=220,c.crit=0.01,data=dados)
plot(fit4)
summary(fit4)


fit5<- gamlss(umid~1,family = SIMPLEX(),n.cyc=220,c.crit=0.01,data=dados)
plot(fit5)
summary(fit5)

fit5<- gamlss(umid~1,family = SIMPLEX(),n.cyc=220,c.crit=0.01,data=dados)
plot(fit5)
summary(fit5)

AIC(fit0,fit3,fit4,fit5)
BIC(fit0,fit2,fit3,fit4,fit5)

x11()
hist(dados$umid,freq=F,xlab = "Soil water content (cm³/cm³)",ylab = "f(y)",cex.lab=1.35,cex.axis=1.35,main="")
curve(dEPE01(x, mu=exp(fit0$mu.coefficients)/(1+exp(fit0$mu.coefficients)),sigma=exp(fit0$sigma.coefficients),nu=exp(fit0$nu.coefficients),tau=exp(fit0$tau.coefficients)), add=T,lwd=4,col="gray15")
curve(dBE(x, mu=exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients))), add=T,lwd=4,col="gray30",lty=2)
curve(dSIMPLEX(x, mu=exp(fit5$mu.coefficients)/(1+exp(fit5$mu.coefficients)),sigma=exp(fit5$sigma.coefficients)), add=T,lwd=4,col="gray45",lty=4)
curve(dGB1(x, mu=exp(fit.GBE$mu.coefficients)/(1+exp(fit.GBE$mu.coefficients)),sigma=exp(fit.GBE$sigma.coefficients)/(1+exp(fit.GBE$sigma.coefficients)),nu=exp(fit.GBE$nu.coefficients),tau=exp(fit.GBE$tau.coefficients)), add=T,lwd=3,col="gray60")
legend("topright",c("logit-EPE","BE","Simplex","GBE"),lwd=c(3,3,3,3),lty=c(1,2,4,1),col=c("gray15","gray30","gray45","gray60"),bty='n',cex=1.2)
box()


plot(ecdf(dados$umid))
curve(pEPE01(x, mu=exp(fit0$mu.coefficients)/(1+exp(fit0$mu.coefficients)),sigma=exp(fit0$sigma.coefficients),nu=exp(fit0$nu.coefficients),tau=exp(fit0$tau.coefficients)), add=T,xlim=c(0,1),lwd=3,col=2)
curve(pBE(x, mu=exp(fit4$mu.coefficients)/(1+exp(fit4$mu.coefficients)),sigma=exp(fit4$sigma.coefficients)/(1+exp(fit4$sigma.coefficients))), add=T,lwd=3,col=4)

########################################################################################
#fit.EPE01.mu<- gamlss(umid ~as.factor(trat),family = EPE01(mu.link="probit"),n.cyc=8,c.crit=0.01,data=dados)
#fit0<- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(10),c.crit=0.01,data=dados)
fit.EPE01 <- gamlss(umid ~1,
                    sigma.formula = ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.EPE01)
summary(fit.EPE01)

fit.EPE01.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                       sigma.formula = ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.EPE01.mu)
summary(fit.EPE01.mu)

fit.EPE01.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(10),c.crit=0.01,data=dados)

plot(fit.EPE01.sig)
summary(fit.EPE01.sig)


fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                           sigma.formula = ~as.factor(trat)+as.factor(prof),
                           family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados)
fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                           sigma.formula = ~as.factor(trat)+as.factor(prof),
                           family = EPE01(),method = CG(1000),c.crit=0.01,data=dados)

plot(fit.EPE01.mu.sig)
summary(fit.EPE01.mu.sig)

plot(dados$umid,fitted(fit.EPE01.mu.sig,what="mu"),xlim=c(min(dados$umid),max(dados$umid)))
lines(c(min(dados$umid),1),c(0,1))
lines(dados$umid,fitted(fit.EPE01.mu.sig,what="mu"),type = "p",pch=20,col=2)
lines(dados$umid,fitted(fit.BE.mu.sig,what="mu"),type = "p",pch=20)
lines(c(min(dados$umid),max(dados$umid)),c(min(dados$umid),max(dados$umid)))



###########################################################################

########################################################################################
#PE
gen.Family(family = PE(mu.link="logit"),type="logit")
fit.PE01 <- gamlss(umid ~1,
                   sigma.formula = ~1,family = EPE01(),method = CG(100),c.crit=0.01,data=dados)

plot(fit.PE01)
summary(fit.PE01)

fit.PE01.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                      sigma.formula = ~1,family = EPE01(),method=CG(100),c.crit=0.01,data=dados)

plot(fit.PE01.mu)
summary(fit.PE01.mu,type="qr")

fit.PE01.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = logitPE,c.crit=0.01,data=dados,method = CG(100))

plot(fit.PE01.sig)
fit.PE01.sig


fit.PE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                          sigma.formula = ~as.factor(trat)+as.factor(prof),
                          family = EPE01(nu.link="identity",tau.link="identity"),method = CG(100),c.crit=0.01,data=dados)

#fit.PE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
#                          sigma.formula = ~as.factor(trat)+as.factor(prof),
#                         family = logitPE(),n.cyc=1000,c.crit=0.01,data=dados)
plot(fit.PE01.mu.sig)
summary(fit.PE01.mu.sig)



###########################################################################
##logit-ENO
fit.ENO01 <- gamlss(umid ~1,
                    sigma.formula = ~1,family = EPE01(nu.link="identity",tau.link="identity"),n.cyc=1000,c.crit=0.01,data=dados,
                    nu.start = 2,nu.fix = T)

plot(fit.ENO01)
summary(fit.ENO01)

fit.ENO.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                     sigma.formula = ~1,family = EPE01(),n.cyc=150,c.crit=0.01,data=dados,
                     nu.start = 2,nu.fix = T)

plot(fit.ENO.mu)
summary(fit.ENO.mu)

fit.ENO.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(1000),c.crit=0.01,data=dados,
                      nu.start = 2,nu.fix = T)

plot(fit.ENO.sig)
summary(fit.ENO.sig)


fit.ENO.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                         sigma.formula = ~as.factor(trat)+as.factor(prof),
                         family = EPE01(tau.link="identity"),method = CG(30),c.crit=0.01,data=dados,
                         nu.start = 2,nu.fix = T,tau.start =4.55)

plot(fit.ENO.mu.sig)
summary(fit.ENO.mu.sig)

LR.test(fit.ENO.mu.sig,fit.PE01.mu.sig)
###########################################################################
###########################################################################
##logit-NO
fit.NO01 <- gamlss(umid ~1,
                   sigma.formula = ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados,
                   nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)

fit.NO01 <- gamlss(umid ~1,
                   sigma.formula = ~1,family = LOGITNO(),n.cyc=1000,c.crit=0.01,data=dados)
plot(fit.NO01)
summary(fit.NO01)

fit.NO.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                    sigma.formula = ~1,family = EPE01(),n.cyc=150,c.crit=0.01,data=dados,
                    nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.NO.mu)
fit.NO.mu

fit.NO.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(1000),c.crit=0.01,data=dados,
                     nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.NO.sig)
fit.NO.sig


fit.NO.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                        sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(nu.link="identity",tau.link="identity"),nmethod=CG(100),c.crit=0.01,data=dados,
                        nu.start = 2,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.NO.mu.sig)
fit.NO.mu.sig

LR.test(fit.NO.mu.sig,fit.PE01.mu.sig)
LR.test(fit.NO.mu.sig,fit.ENO.mu.sig)
###########################################################################
###########################################################################
##logit-EL
fit.EL01 <- gamlss(umid ~1,
                   sigma.formula = ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados,
                   nu.start = 1,nu.fix = T)

plot(fit.EL01)
summary(fit.EL01)

fit.EL.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                    sigma.formula = ~1,family = EPE01(),n.cyc=150,c.crit=0.01,data=dados,
                    nu.start = 1,nu.fix = T)

plot(fit.EL.mu)
fit.EL.mu

fit.EL.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(1000),c.crit=0.01,data=dados,
                     nu.start = 1,nu.fix = T)

plot(fit.EL.sig)
fit.EL.sig


fit.EL.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                        sigma.formula = ~as.factor(trat)+as.factor(prof),
                        family = EPE01(nu.link="identity",tau.link="identity"),
                        method=CG(100),c.crit=0.01,data=dados,
                        nu.start = 1,nu.fix = T)

plot(fit.EL.mu.sig)
fit.EL.mu.sig

LR.test(fit.EL.mu.sig,fit.PE01.mu.sig)

###########################################################################

###########################################################################
##logit-L
fit.L01 <- gamlss(umid ~1,
                  sigma.formula = ~1,family = EPE01(),n.cyc=1000,c.crit=0.01,data=dados,
                  nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.L01)
summary(fit.L01)

fit.L.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                   sigma.formula = ~1,family = EPE01(),n.cyc=150,c.crit=0.01,data=dados,
                   nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.L.mu)
fit.L.mu

fit.L.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = EPE01(),method=CG(1000),c.crit=0.01,data=dados,
                    nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.L.sig)
fit.L.sig


fit.L.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                       sigma.formula = ~as.factor(trat)+as.factor(prof),
                       family = EPE01(nu.link="identity",tau.link="identity"),method = CG(1000),c.crit=0.01,data=dados,
                       nu.start = 1,nu.fix = T,tau.start = 1,tau.fix = T)

plot(fit.L.mu.sig)
fit.L.mu.sig


LR.test(fit.L.mu.sig,fit.PE01.mu.sig)
LR.test(fit.L.mu.sig,fit.EL.mu.sig)
###########################################################################
##GBE
fit.GBE <- gamlss(umid ~1,
                  sigma.formula = ~1,family = GB1,n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.GBE)
summary(fit.GBE)

fit.GBE.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                     sigma.formula = ~1,family = GB1(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.GBE.mu)
summary(fit.GBE.mu)

fit.GBE.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = GB1(),c.crit=0.01,data=dados)

plot(fit.GBE.sig)
summary(fit.GBE.sig)


fit.GBE.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                         sigma.formula = ~as.factor(trat)+as.factor(prof),family = GB1(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.GBE.mu.sig)
summary(fit.GBE.mu.sig)


###########################################################################
###########################################################################
##BE
fit.BE <- gamlss(umid ~1,
                 sigma.formula = ~1,family = BE,n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.BE)
summary(fit.BE)

fit.BE.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                    sigma.formula = ~1,family = BE(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.BE.mu)
summary(fit.BE.mu)

fit.BE.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = BE(),c.crit=0.01,data=dados)

plot(fit.BE.sig)
summary(fit.BE.sig)


fit.BE.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                        sigma.formula = ~as.factor(trat)+as.factor(prof),family = BE(),n.cyc=1000,method = CG(100),data=dados)

plot(fit.BE.mu.sig)
summary(fit.BE.mu.sig)




fit.Simplex <- gamlss(umid ~1,
                      sigma.formula = ~1,family = SIMPLEX(),method = CG(100),c.crit=0.01,data=dados)

plot(fit.Simplex)
summary(fit.Simplex)

fit.Simplex.mu <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                         sigma.formula = ~1,family = SIMPLEX(),n.cyc=1000,c.crit=0.01,data=dados)

plot(fit.Simplex.mu)
summary(fit.Simplex.mu)

fit.Simplex.sig <- gamlss(umid ~1,  sigma.formula = ~as.factor(trat)+as.factor(prof),family = SIMPLEX(),c.crit=0.01,data=dados)

plot(fit.Simplex.sig)
summary(fit.Simplex.sig)


fit.Simplex.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                             sigma.formula = ~as.factor(trat)+as.factor(prof),family = SIMPLEX(),n.cyc=1000,method = CG(100),data=dados)

plot(fit.Simplex.mu.sig)
summary(fit.Simplex.mu.sig)


###########################################################################

wp(fit.EPE01.mu.sig,ylim.all = 1.5)
wp(fit.Simplex.mu.sig,ylim.all = 1.5)
wp(fit.BE.mu.sig,ylim.all = 1.5)
wp(fit.GBE.mu.sig,ylim.all = 1.5)


x11()
wp1(fit.Simplex.mu.sig,ylim.all = 1.5,cor1 = "gray65",lwd=4)

x11()
wp1(fit.GBE.mu.sig,ylim.all = 1.5,cor1 = "gray65",lwd=4)

x11()
wp1(fit.BE.mu.sig,ylim.all = 1.5,cor1 = "gray65",lwd=4)



fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                           sigma.formula = ~as.factor(trat)+as.factor(prof),
                           family = EPE01(),method = CG(100),c.crit=0.01,data=dados)

summary(fit.EPE01.mu.sig)

n <- length(dados$umid)
index <- 1:n

Res.q1 <- fit.EPE01.mu.sig$residuals
x11()
plot(index,Res.q1,col="gray11",pch=20,
     ylab="Quantile residuals",xlab="Index",main = "",cex.lab=1.2,
     ylim=c(-4,4))
abline(h=-3,lwd=4,lty=2,col="gray28")
abline(h=0,lwd=3,lty=1,col="gray50")
abline(h=3,lwd=4,lty=2,col="gray28")
identify(index,Res.q1)


Fitted.1 <- fitted(fit.EPE01.mu.sig)
x11()
plot(Fitted.1,Res.q1,col="gray11",pch=20,
     ylab="Quantile residuals",xlab="Fitted values for mu",main = "",cex.lab=1.2,
     ylim=c(-4,4))
abline(h=-3,lwd=4,lty=2,col="gray28")
abline(h=0,lwd=3,lty=1,col="gray50")
abline(h=3,lwd=4,lty=2,col="gray28")

Fitted.2 <- fitted(fit.EPE01.mu.sig,what="sigma")
x11()
plot(Fitted.2,Res.q1,col="gray11",pch=20,
     ylab="Quantile residuals",xlab="Fitted values for sigma",main = "",cex.lab=1.2,
     ylim=c(-4,4))
abline(h=-3,lwd=4,lty=2,col="gray28")
abline(h=0,lwd=3,lty=1,col="gray50")
abline(h=3,lwd=4,lty=2,col="gray28")


x11()
wp1(fit.EPE01.mu.sig,ylim.all = 1.5,cor1 = "gray65",lwd=4)



Res.q1 <- fit.EPE01.mu.sig$residuals
Res.qo <- sort(Res.q1)
index <- 1:length(dados$umid)
plot(index,Res.q1, ylim=c(-4,4))

iter<-0  
j<-1
B<-1000

set.seed(78)
mrq <- matrix(0,  ncol = B, nrow = n)
while(j<B+1){
  simula <- rEPE01(n,mu= fitted(fit.EPE01.mu.sig,what="mu"),sigma=fitted(fit.EPE01.mu.sig,what="sigma"),nu=fitted(fit.EPE01.mu.sig,what="nu"),tau=fitted(fit.EPE01.mu.sig,what="tau"))
  m1s <- try(gamlss(simula ~trat  + prof,sigma.formula =~trat  + prof,family=EPE01(),c.crit=0.01,n.cyc=25,data=dados
                    ,mu.start =fit.EPE01.mu.sig$mu.fv,sigma.start =fit.EPE01.mu.sig$sigmat.fv ))
  if((class(m1s) != "try-error")==T){
    
    Res.qs <- m1s$residuals
    mrq[,j] <- Res.qs
    j=j+1
  }
  cat("iteration = ", iter <- iter + 1, j,"\n")
}


Res.q.Ords <- apply(mrq[,1:340],2,sort)

i <- 1:n
Z <- qnorm((i - 3/8) / (n + 1/4))
#Z <- qnorm((i - 1/8) / (2*n + 1/2))
#Z <- qnorm((i +n - 1/8) / (2*n + 1/2))

rqi.m  <-  apply(Res.q.Ords,1, mean)
rqi.min <- apply(Res.q.Ords,1, min)
rqi.max <- apply(Res.q.Ords,1, max)


(res.out <- sum(c(sum(Res.qo>rqi.max),sum(Res.qo<rqi.min))))
(per.out <- round(res.out/n*100,2))

x11()
dd <- qqnorm(Res.qo,pch=20,col="gray1",ylim = c(-3.5,3.5),xlim=c(-3.5,3.5),
             ylab="Quantile Residuals",xlab="N(0,1) quantiles",main = "",cex.lab=1.35,cex.axis=1.35)
lines(Z,rqi.max,col="gray40",lwd=3)
lines(Z,rqi.m, lty = 2,col="gray60",lwd=3)
lines(Z,rqi.min,col="gray40",lwd=3)
identify(dd$x,dd$y)
legend("topleft", c(paste("Total points:",n), paste("Points out of envelope:",res.out,"(",per.out,"%)")), bty="n", cex=1.5)

###########################################

n=length(y)
#####Global influence
vet <- numeric()
estima1 <- matrix(c(rep(0,(n)*6)), ncol = 6, nrow = n)

for(i in 1:n){
  m1 <-gamlss(y[-i]~x1[-i]+x2[-i],
              family = GOLLBE(),n.cyc=300,
              c.crit=0.01)
  v <- as.numeric(logLik(m1)); vet <- c(vet,v)
  estima1[i,] <- as.matrix(c(m1$mu.coefficients,m1$sigma.coefficients,
                             m1$nu.coefficients,m1$tau.coefficients))
}

v1<- as.numeric(logLik(fit11))
vcomp <- c(rep(v1,length(vet)))
LDp1 <- (2*(vcomp-vet))

x11()
par(pch=19, col="black")
plot(1:length(LDp1),abs(LDp1),ylim=c(0,6.5),ylab="Likelihood distance",xlab="Index",col="1",cex.lab=1.5,cex.axis=1.5)
points(1:n,abs(LDp1[1:n]),col="gray25")
lines(1:n,abs(LDp1[1:n]),type="h",col="gray25")
identify(1:length(vet),abs(LDp1))


hh <- vcov(fit11)

dim(hh)
Dif <- matrix(c(rep(0,(n)*6)), ncol = 6, nrow = n)
for(i in 1:n){
  Dif[i,] <- estima1[i,]-c(fit11$mu.coefficients,fit11$sigma.coefficients,
                           fit11$nu.coefficients,fit11$tau.coefficients)
}

dim(Dif)
GDi <- as.numeric()
for(i in 1:n){
  GDi[i] <- t(Dif[i,])%*%(hh)%*%Dif[i,]
}

index = 1:n
x11()
par(pch=19, col=1)
plot(index,GDi,ylab = "Generalized Cook distance",xlab="Index",
     cex.lab=1.5,cex.axis=1.2,ylim = c(0,0.25))
points(1:n,GDi[1:n],col="gray25")
lines(1:n,GDi[1:n],type="h",col="gray25")
identify(index,GDi)



(trat.1<- relevel(dados$trat, ref = "Cerrado"))
(prof.1<- relevel(dados$prof, ref = "0-10"))
fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat.1)+as.factor(prof.1),
                           sigma.formula = ~as.factor(trat.1)+as.factor(prof.1),
                           family = EPE01(),method=CG(200),c.crit=0.01,data=dados)
plot(fit.EPE01.mu.sig)
#fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
#                          sigma.formula = ~as.factor(trat)+as.factor(prof),
#                          family = EPE01(),method = CG(1000),c.crit=0.01,data=dados)


summary(fit.EPE01.mu.sig)

cbind(round(fit.EPE01.mu.sig$mu.coefficients,3),round(confint(fit.EPE01.mu.sig,what = "mu"),3))

cbind(round(fit.EPE01.mu.sig$sigma.coefficients,3),round(confint(fit.EPE01.mu.sig,what = "sigma"),3))

tapply(dados$umid, dados$trat, mean)
tapply(dados$umid, dados$trat, sd)

0.0228-0.009

plot(c(0.2,0.42), c(0,1),type="n", xlab="Time (months)",
     ylab="Suvival function",cex.lab=1.2,cex.axis=1.2)
lines(ecdf(dados$umid[dados$trat=="t1"]),col=1)
lines(ecdf(dados$umid[dados$trat=="t2"]),col=2)
lines(ecdf(dados$umid[dados$trat=="t3"]),col=3)
lines(ecdf(dados$umid[dados$trat=="t4"]),col=4)
lines(ecdf(dados$umid[dados$trat=="t5"]),col=5)
lines(ecdf(dados$umid[dados$trat=="t6"]),col=6)
lines(ecdf(dados$umid[dados$trat=="Cerrado"]),col=7)


exp(-1.169 -0.925)

0.298-0.284

exp(0.238-1.466)/(1+exp(0.238-1.466))


require(ggplot2)

dim(dados)

fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                           sigma.formula = ~as.factor(trat)+as.factor(prof),
                           family = EPE01(),n.cyc=1000,c.crit=0.001,data=dados)

fit.EPE01.mu.sig <- gamlss(umid ~as.factor(trat)+as.factor(prof),
                           sigma.formula = ~as.factor(trat)+as.factor(prof),
                           family = EPE01,method=CG(1000),c.crit=0.0001,data=dados)
summary(fit.EPE01.mu.sig)

cbind(round(fit.EPE01.mu.sig$mu.coefficients,3),round(confint(fit.EPE01.mu.sig,what = "mu"),3))

d <- data.frame(SWC=c(dados$umid,fitted(fit.EPE01.mu.sig)),Treatments = as.factor(c(dados$trat,dados$trat)),Observations = as.factor(c(rep("True",105),rep("Fitted",105))))

x11()
ggplot(d, aes(Treatments, SWC)) +
  geom_boxplot(aes(fill = Observations)) +
  labs(x = "Treatments", y = "MDC (NTU/gl-1)",cex=1.5) +
  theme_bw()+   scale_fill_manual(values = c("gray95", "gray45"))
