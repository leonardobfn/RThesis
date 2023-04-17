log.f.cond.beta <- function(theta,alpha,lambda,y,x,w,ncx,ncv){
  #theta = c(par.cov.start.marginal, alpha.start)
  Beta = theta
  lambda = lambda #theta[c((1+ncx):(ncx+ncv))]
  alpha = alpha
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% Beta
  a = exp(xbeta)
  b = 1/delta

  log.f.cond. = NULL

  gt <- extraDistr::dkumar(y,a,b)
  Gt <- extraDistr::pkumar(y,a,b,lower.tail = F)
  r <- gt/Gt # derivada de -log(Gt)
  LAMBDA <- -log(Gt)
  log.f1 <- -LAMBDA[1]^alpha+log(alpha)+(alpha-1)*log(LAMBDA[1])+log(r[1])
  for(t in 2:length(a)){
    log.f.cond.[t] <- log(r[t]) + (1-alpha)*log(LAMBDA[t-1])+LAMBDA[t-1]^(alpha)-(LAMBDA[t-1]+LAMBDA[t])^alpha+
      (alpha-2)*log(LAMBDA[t-1]+LAMBDA[t])+log(alpha*( LAMBDA[t-1]+LAMBDA[t]  )^alpha+(1-alpha))

  }
  log.f.cond. <- log.f.cond.[is.finite(log.f.cond.)==T]
  l = c(log.f1,sum(log.f.cond.))
  return( sum(l[is.finite(l)==T]))

}
