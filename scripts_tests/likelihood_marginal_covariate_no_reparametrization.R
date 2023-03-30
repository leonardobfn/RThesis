

rm(list=ls())
Q=function(theta){
  #theta = start
  Beta = theta[1:3]
  lambda = theta[c(4:5)]
  alpha = theta[6]
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  a = exp(xbeta)
  b = exp(wlambda)
  l=sum( log(a) + log(alpha) + alpha*log(b) + (a-1)*log(y) - log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) - (-b*log(1-y^(a)))^(alpha)   )
  return(l)
}


N=c(1,12,24,36,48)*120
emv = NULL
B=50

betas = as.matrix(c(1,1,.5))
lambdas = as.matrix(c(1,-1))
alpha.real = .75
#start=c(betas,lambdas,alpha.real)
start=c(.1,.1,.1,.1,-.1,.5)
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
a.real = exp(xbetas)
b.real = exp(wlambdas)


  for(b in 1:B){
    u = runif(N[k])
    y=(pweibull(-log(1-u), shape=1/alpha.real, scale = b.real^(alpha.real), lower.tail = TRUE, log.p = FALSE))^(1/a.real)
    y[y>.9999999999999999] <- .9999999999999999
    hist(y)
    results = try(
      optim(start,Q,method = "L-BFGS-B",control=list(fnscale=-1),hessian = T),silent = T)
    if(class(results)=="try-error"){
      results <- list()
      results$par <- rep("NA",length(start))
    }
    emv.b[b,] = results$par
  }
  emvs = data.frame(N=N[k],emv.b = emv.b)
  emv = rbind(emv,emvs)
}
emv$N.f = factor(emv$N,levels = N)
write.table(emv,"emv_covariates_no_parametrization.txt",col.names = T)
emv.results = read.table("emv_covariates_no_parametrization.txt")
alpha.real = .75
boxplot(emv.results$emv.b.6~emv.results$N)
abline(h=alpha.real)
H=results$hessian
s.d = sqrt(diag( -solve(H)))
EMV = results$par
EMV -1.96*s.d
EMV +1.96*s.d

