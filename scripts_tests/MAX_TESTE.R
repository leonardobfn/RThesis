rm(list=ls())
library(lattice)
beta = as.matrix(c(1,0.15,-0.1))
lambda = as.matrix(c(1,15))
alpha=.75
ar=0.5
ma=0.5
nn=1
j=1:(120*nn)
x1=rep(1,120*nn)
x2=sin(2*pi*j/12)
x3=cos(2*pi*j/12)
w1=rep(1,120*nn)
w2=rep(c(1,1,1,1,1,1,0,0,0,0,0,0),10)
x=cbind(x1,x2,x3)
w=cbind(w1,w2)
xbeta=x%*%beta
wlambda = w%*%lambda
m = exp(-exp(-xbeta[1]))
delta=exp(wlambda)
a=pweibull(-log(0.5), shape=1/alpha, scale = delta[1]^(alpha), lower.tail = TRUE, log.p = TRUE)/log(m)

R=1
emv=matrix(0,R,8)
set.seed(30)
for(k in 1:R){

  u=runif(1)
  y1 = (pweibull(-log(1-u), shape=1/alpha, scale = delta[1]^(alpha), lower.tail = TRUE, log.p = FALSE))^(1/a)
  y=NULL
  mi=NULL
  ai=NULL
  y[1]=y1
  mi[1] = m
  ai[1]=pweibull(-log(0.5), shape=1/alpha, scale = delta[1]^(alpha), lower.tail = TRUE, log.p = TRUE)/log(mi[1])
  for(t in 2:(120*nn)){
    mi[t] = exp(-exp(-( xbeta[t] + ar*(-log(-log(y[t-1])) - xbeta[t-1]) + ma*(-log(-log(y[t-1])) + log(-log(mi[t-1])))  )))
    ai[t]=pweibull(-log(0.5), shape=1/alpha, scale = delta[t]^(alpha), lower.tail = TRUE, log.p = TRUE)/log(mi[t])
    repeat{
      u=runif(1)
      y[t]=(pweibull(-log(1-u), shape=1/alpha, scale = delta[t]^(alpha), lower.tail = TRUE, log.p = FALSE))^(1/ai[t])
      if(format(y[t], scientific = FALSE)!=1){
        break
      }
    }
  }

  Q=function(theta){
    beta=c(theta[1],theta[2],theta[3])
    phi =theta[4]
    gama=theta[5]
    lambda=c(theta[6],theta[7])
    alpha=theta[8]
    #alpha=alpha
    xbeta = x%*%beta
    wlambda = w%*%lambda
    delta = exp(wlambda)
    for(t in 2:(120*nn)){
      mi[t] = exp(-exp(-(xbeta[t]  + phi*(-log(-log(y[t-1])) - xbeta[t-1]) + gama*(-log(-log(y[t-1])) + log(-log(mi[t-1])))  )))
      ai[t]=pweibull(-log(0.5), shape=1/alpha, scale = delta[t]^(alpha), lower.tail = TRUE, log.p = TRUE)/log(mi[t])
    }

    l=sum( log(ai) + log(alpha) + alpha*log(delta) + (ai-1)*log(y) - log(1-y^(ai)) + (alpha-1)*log(-log(1-y^(ai))) - (-delta*log(1-y^(ai)))^(alpha)   )
    return(l)
  }
  start=c(1,0.15,-0.1,0.5,0.5,1,1,0.5)
  EMVteta=try(optim(start,Q, method = "BFGS", control= list(fnscale=-1),
                   ),TRUE)
  # EMVteta=try(optim(start,Q, method = "L-BFGS-B", control= list(fnscale=-1),
  #                   lower=c(-.5,-.5,-.5,-.5,-.5,.-5,.-5,0.45),
  #                   upper=c(1.5,1.5,1.5,1.5,1.5,1.5,50,0.90)),TRUE)
  EMVteta
  if (class(EMVteta) == "try-error") {
    EMVteta <- list()
    EMVteta$par <- rep(NA,length(start))
    #try_error <- 1
  }

  emv[k,]=EMVteta$par
  #rm(EMVteta)
  cat(k,emv[k,],"\r")
} ###fim k
EMVest=colMeans(emv[1:8,],na.rm = T)
apply(emv,2,quantile,c(0.025,0.5,0.975),na.rm=T) %>% t()
EMVest
