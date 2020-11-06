EPE01 <- function (mu.link = "logit", sigma.link="log", nu.link = "log", tau.link = "log")
{
  mstats <- checklink(   "mu.link", "Exponentiated Power Exponential (0,1)", substitute(mu.link), 
                         c("logit", "log", "identity", "own"))
  dstats <- checklink("sigma.link", "Exponentiated Power Exponential (0,1)", substitute(sigma.link), 
                      c("inverse", "log", "identity", "own"))
  vstats <- checklink(   "nu.link", "Exponentiated Power Exponential (0,1)", substitute(nu.link),    
                         c("1/nu^2", "log", "identity", "own"))
  tstats <- checklink(  "tau.link", "Exponentiated Power Exponential (0,1)", substitute(tau.link),   
                        c("1/tau^2", "log", "identity", "own")) 
  structure(
    list(family = c("EPE01", "Exponentiated Power Exponential (0,1)"),
         parameters = list(mu=TRUE, sigma=TRUE, nu=TRUE, tau=TRUE), 
         nopar = 4, 
         type = "Continuous",
         mu.link = as.character(substitute(mu.link)),  
         sigma.link = as.character(substitute(sigma.link)), 
         nu.link = as.character(substitute(nu.link)), 
         tau.link = as.character(substitute(tau.link)), 
         mu.linkfun = mstats$linkfun, 
         sigma.linkfun = dstats$linkfun, 
         nu.linkfun = vstats$linkfun,
         tau.linkfun = tstats$linkfun,  
         mu.linkinv = mstats$linkinv, 
         sigma.linkinv = dstats$linkinv,
         nu.linkinv = vstats$linkinv,
         tau.linkinv = tstats$linkinv, 
         mu.dr = mstats$mu.eta, 
         sigma.dr = dstats$mu.eta, 
         nu.dr = vstats$mu.eta,
         tau.dr = tstats$mu.eta, 
         
         dldm = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           dldm
         },
         d2ldm2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = T),"mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient")) 
           d2ldm2 = -dldm * dldm
         },     
         dldd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok  
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           dldd
         } ,
         d2ldd2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))
           d2ldd2 = -dldd*dldd
           d2ldd2 
         },   
         dldv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           dldv 
         },
         d2ldv2 = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok 
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"nu", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldv2 = -dldv * dldv
         },
         dldt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"tau", delta = 1e-04)
           dldt = as.vector(attr(nd1, "gradient"))           
           dldt
         } ,
         d2ldt2 = function(y,mu,sigma,nu,tau){ #----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log =T),"tau", delta = 1e-04)
           dldt = as.vector(attr(nd1, "gradient"))  
           d2ldt2 = -dldt * dldt
           d2ldt2
         },
         d2ldmdd = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu, tau,log= TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu, tau,log=TRUE), "sigma", delta = 1e-04)
           dldd = as.vector(attr(nd1, "gradient"))           
           d2ldmdd = -dldm * dldd
           d2ldmdd               
         },
         d2ldmdv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv			
         },
         
         d2ldmdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "mu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv         },
         
         d2ldddv = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv	
         },
         d2ldddt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "sigma", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv
         },
         d2ldvdt = function(y,mu,sigma,nu,tau){#----------------------------------------------------- ok
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "nu", delta = 1e-04)
           dldm = as.vector(attr(nd1, "gradient"))
           nd1 = gamlss:::numeric.deriv(dEPE01(y, mu, sigma, nu,tau,log = TRUE), "tau", delta = 1e-04)
           dldv = as.vector(attr(nd1, "gradient"))
           d2ldmdv = -dldm * dldv
           d2ldmdv 
         },
         #----------------------------------------------------- ok
         G.dev.incr  = function(y,mu,sigma,nu,tau,...) 
         { 
           -2*dEPE01(y,mu,sigma,nu,tau,log=TRUE)
         } ,                     
         rqres = expression(   
           rqres(pfun="pEPE01", type="Continuous", y=y, mu=mu, sigma=sigma, nu=nu, tau=tau)) ,
         mu.initial = expression(mu <- rep(mean(y), length(y))),
         #mu.initial = expression( mu <- rep(mean(y), length(y))), 
         #sigma.initial = expression(sigma <- rep(sd(log(y/(1-y))),length(y))), # 
         sigma.initial = expression( sigma <- rep(sd(y),length(y)) ), #
         nu.initial = expression( nu <- rep(2, length(y))), 
         tau.initial = expression(tau <-rep(2, length(y))), 
         mu.valid = function(mu) all(mu > 0 & mu < 1)  , 
         sigma.valid = function(sigma)  all(sigma > 0),
         nu.valid = function(nu) all(nu > 0), 
         tau.valid = function(tau) all(tau > 0), 
         y.valid = function(y)  all(y > 0 & y < 1) 
    ),
    class = c("gamlss.family","family"))
}
#-----------------------------------------------------------------
dEPE01 <- function(x, mu = 0.5, sigma = 1, nu = 1, tau = 1, log = FALSE){
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("tau must be positive", "\n", ""))  
  
  z <- log(x)-log(1-x)
  pdfep <- dPE(z, mu = log(mu/(1-mu)), sigma = sigma, nu = nu)
  cdfep <- pPE(z, mu = log(mu/(1-mu)), sigma = sigma, nu = nu)
  fy1 <- (tau*pdfep*cdfep^(tau-1))*(1/(x*(1-x)))
  if(log==FALSE) fy<-fy1 else fy<-log(fy1)
  fy
}    
#-----------------------------------------------------------------  
pEPE01 <- function(q, mu = 0.5, sigma = 1, nu = 1, tau = 1, lower.tail = TRUE, log.p = FALSE){  
  if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(nu <= 0))  stop(paste("tau must be positive", "\n", ""))  
  
  qz <- log(q)-log(1-q)
  cdfep <-  pPE(qz, mu = log(mu/(1-mu)), sigma = sigma, nu = nu)
  cdf1 <- cdfep^tau
  
  if(lower.tail==TRUE) cdf<-cdf1 else cdf<- 1-cdf1
  if(log.p==FALSE) cdf<- cdf else cdf<- log(cdf)
  cdf
}
#-----------------------------------------------------------------  

qEPE01 <-  function(p, mu=0.5, sigma=1, nu=1, tau=1, lower.tail = TRUE, log.p = FALSE){   
  if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(p < 0)|any(p > 1))  stop(paste("p must be between 0 and 1", "\n", ""))    
  if (log.p==TRUE) p <- exp(p) else p <- p
  if (lower.tail==TRUE) p <- p else p <- 1-p
  #p <- runif(1000)
  u <- p^(1/tau)
  qz <- qPE(p=u, mu=log(mu/(1-mu)), sigma=sigma, nu=nu)
  qy <- 1/(1+exp(-qz))
  qy
}
#-----------------------------------------------------------------  
#-----------------------------------------------------------------  


rEPE01 <- function(n, mu=0.5, sigma=1, nu=1, tau=1){
  if (any(mu < 0) | any(mu > 1))  stop(paste("mu must be between 0 and 1", "\n", "")) 
  if (any(sigma <= 0))  stop(paste("sigma must be positive", "\n", "")) 
  if (any(nu <= 0))  stop(paste("tau must be positive", "\n", ""))  
  if (any(tau <= 0))  stop(paste("tau must be positive", "\n", ""))
  u <- runif(n ,0,1)
  uni <- u^(1/tau)
  rz <- qPE(p=uni, mu=log(mu/(1-mu)), sigma=sigma, nu=nu)
  rz
  #mu=0.5; sigma=0.3; nu=1; tau=0.15
  #rz <- rEPE(n, mu=log(mu/(1-mu)), sigma=sigma, nu=nu,tau=tau)
  ry <- 1/(1+exp(-rz))
  
  ry1 <- NULL
  for(i in 1:length(ry)){
    if(ry[i]==0){ry1[i]=0.00001}
    if(ry[i]==1){ry1[i]=0.99999}
    if(ry[i]>0&ry[i]<1){ry1[i]=ry[i]}
  }
  ry1
  
}


