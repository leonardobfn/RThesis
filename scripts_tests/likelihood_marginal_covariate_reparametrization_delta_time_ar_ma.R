

rm(list=ls())
Q=function(theta){
  #theta = start
  Beta = theta[1:3]
  ar = theta[4]
  ma = theta[5]
  lambda = theta[c(6:7)]
  alpha = theta[8]
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  delta = exp(wlambda)
  b = 1/delta
  a=NULL

  #a[1] = a.real[1]
  a[1] = exp(xbeta[1])
  for(t in 2:N[k]){

    a[t] = exp(
      xbeta[t] + ar*(-log(-log(y[t-1])) - xbeta[t-1]) + ma*(-log(-log(y[t-1])) - log(a[t-1]) )
    )

  }

  # a = a[-1]
  # m = m[-1]
  # y = y[-1]

  l = sum( log(a) + log(alpha) +
             alpha*log(b) + (a-1)*log(y) -
             log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) -
             (-b*log(1-y^(a)))^(alpha))
  return(l)
}

N=c(1,12,24,36,48)*120
#N=240
emv = NULL
B=100

betas = as.matrix(c(1,-1,-.5))
lambdas = as.matrix(c(-1,-3))
alpha.real = .55
ar.real = .5
ma.real = .5
#start=c(betas,ar.real,ma.real,lambdas,alpha.real)
start=c(1,0.15,-0.1,0.5,0.5,1,1,0.5)
emv.b = matrix(0,B,length(start))
for(k in 1:length(N)){
  #k=5
  j=1:N[k]
  x1=rep(1,N[k])
  x2=sin(2*pi*j/12)
  x3=cos(2*pi*j/12)
  w1=rep(1,N[k])
  w2=rep(c(1,1,1,1,1,1,0,0,0,0,0,0),N[k]/12)
  x=cbind(x1,x2,x3)
  w=cbind(w1,w2)
  xbetas=x%*%betas
  wlambdas = w%*%lambdas
  delta.real = exp(wlambdas)
  b.real = 1/delta.real
  a.real = NULL
  y = NULL

  a.real[1] = exp(xbetas)[1]

  for(b in 1:B){

    y[1] = (pweibull(-log(1-runif(1)), shape=1/alpha.real, scale = b.real[1]^(alpha.real), lower.tail = TRUE, log.p = FALSE))^(1/a.real[1])
    for(t in 2:N[k]){
      a.real[t] = exp(
        xbetas[t] + ar.real*(-log(-log(y[t-1])) - xbetas[t-1]) + ma.real*(-log(-log(y[t-1])) - log(a.real[t-1]))
      )
      repeat{
        u=runif(1)
        y[t]=(pweibull(-log(1-u), shape=1/alpha.real, scale = b.real[t]^(alpha.real), lower.tail = TRUE, log.p = FALSE))^(1/a.real[t])
        if(format(y[t], scientific = FALSE)!=1){
          break
        }
      }
    }
    hist(y)
    plot.ts(y)
    results = try(
      optim(start,Q,method = "BFGS",control=list(fnscale=-1)),silent = T)
    results
    # results = try(
    #   optim(start,Q,method = "L-BFGS-B",control=list(fnscale=-1)),silent = T)
    # results
    if(class(results)=="try-error"){
      results <- list()
      results$par <- rep("NA",length(start))
    }
    results
    emv.b[b,] = results$par
  }
  emvs = data.frame(N=N[k],emv.b = emv.b)
  emv = rbind(emv,emvs)
}
emv$N.f = factor(emv$N,levels = N)
write.table(emv,"likelihood_marginal_covariate_reparametrization_delta_time_ar_ma.txt",col.names = T)
emv.results = read.table("likelihood_marginal_covariate_reparametrization_delta_time_ar_ma.txt")
alpha.real = .55
boxplot(emv.results$emv.b.8~emv.results$N)
abline(h= .55)
H=results$hessian
s.d = sqrt(diag( -solve(H)))
EMV = results$par
EMV -1.96*s.d
EMV +1.96*s.d

require(tidyr)
emv.results[,c(1,9)] %>% tapply(N,mean,na.rm=T)
