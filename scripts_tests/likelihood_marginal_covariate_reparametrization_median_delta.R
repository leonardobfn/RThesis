

rm(list=ls())
Q=function(theta){
  #theta = start
  Beta = theta[1:3]
  lambda = theta[c(4:5)]
  alpha = theta[6]
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  m = exp(-exp(-( xbeta )))
  delta = exp(wlambda)
  b = 1/delta
  a = pweibull(-log(0.5), shape=1/alpha, scale = b^(alpha), lower.tail = TRUE, log.p = TRUE)/log(m)

  l = sum( log(a) + log(alpha) +
             alpha*log(b) + (a-1)*log(y) -
             log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) -
             (-b*log(1-y^(a)))^(alpha))
  return(l)
}

N=c(1,12,24,36,48)*120
emv = NULL
B=60

betas = as.matrix(c(1,-1,-.5))
lambdas = as.matrix(c(-1,-3))
alpha.real = .25
#start=c(betas,lambdas,alpha.real)
start=c(.1,.1,.1,.1,-.1,.5)
emv.b = matrix(0,B,length(start))
for(k in 1:length(N)){
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
  m.real = exp(-exp(-xbetas))
  a.real = pweibull(-log(0.5), shape=1/alpha.real, scale = b.real^(alpha.real), lower.tail = TRUE, log.p = TRUE)/log(m.real)

  for(b in 1:B){
      u = runif(N[k])
      y = (pweibull(-log(1-u), shape=1/alpha.real, scale = b.real^(alpha.real), lower.tail = TRUE, log.p = FALSE))^(1/a.real)
      y[y>.9999999999999999] <- .9999999999999999
      hist(y)
    results = try(
      optim(start,Q,method = "L-BFGS-B",control=list(fnscale=-1),hessian = T),silent = T)
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
write.table(emv,"emv_covariates_reparametrization_median_delta.txt",col.names = T)
emv.results = read.table("emv_covariates_reparametrization_median_delta.txt")
alpha.real = .25
boxplot(emv.results$emv.b.6~emv.results$N)
abline(h=alpha.real)
H=results$hessian
s.d = sqrt(diag( -solve(H)))
EMV = results$par
EMV -1.96*s.d
EMV +1.96*s.d

