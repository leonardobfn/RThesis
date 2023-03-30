

rm(list=ls())
Q=function(theta){
  a = theta[1]
  b = theta[2]
  alpha = theta[3]
  l=sum( log(a) + log(alpha) + alpha*log(b) + (a-1)*log(y) - log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) - (-b*log(1-y^(a)))^(alpha)   )
  return(l)
}


a.real = 2
b.real = 2
alpha.real = .65
N=c(120,500,1000,2000,3000)
#N=3000
N.f = factor(N,levels = N)
emv = NULL
B=50
emv.b = matrix(0,B,3)
for(k in 1:length(N)){
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
boxplot(emv$emv.b.3~emv$N)
abline(h=alpha.real)
H=results$hessian
s.d = sqrt(diag( -solve(H)))
EMV = results$par
EMV -1.96*s.d
EMV +1.96*s.d

