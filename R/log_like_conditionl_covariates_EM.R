log.like_conditionl_covariates_EM = function(theta, y, x, w, Etil1,ncx,ncv) {
  # x: betas covariates
  # w: lambdas covariates
  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)] # real
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% Beta
  a = exp(xbeta)
  Gb = pkumar(y,
              a,
              b = 1 / delta,
              lower.tail = FALSE,
              log.p = FALSE)
  d = dkumar(y, a, b = 1 / delta, log = FALSE)
  r  = d / Gb
  l = Etil1 * log(Gb)
  l = l[is.finite(l)]
  r. = r[is.finite(r)]
  L = sum(l + log(r.), na.rm = T)

  return(L)
}
