log.kumar.gbase = function(theta, x, w, y,ncx,ncv) {
  #theta = par.covariates.start;x=cov_a;w=cov_delta;alpha=.5

  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1/delta
  a = exp(xbeta)


  l = extraDistr::dkumar(
    x = y,
    a = a,
    b = b,
    log = T
  )

  l = l [is.finite(l)==T]
  return(sum(l))
}
