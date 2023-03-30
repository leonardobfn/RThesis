rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

##log_like_gbase-----

log.like.kumar = function(theta, x, w, y) {
  #theta = start

  Beta = theta[1:3]
  lambda = theta[c(4:5)]
  ar = theta[6]
  ma = theta[7]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta
  a = mediana = NULL
  a[1] = exp(xbeta[1])
  mediana[1] = log(a[1]) - log(-log(1 - .5 ^ (b[1])))
  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1]) - xbetas[t - 1])) +
                 ma * (-log(-log(y[t - 1])) - log(a[t - 1])))
    #ma*(-log(-log(y[t-1])) - mediana[t-1]))
    mediana[t] <-
      log(a[t]) - log(-log(1 - .5 ^ (1 / b[t])))
  }

  l = sum(extraDistr::dkumar(
    x = y,
    a = a,
    b = b,
    log = T
  ))
  return(l)
}

## log_like_marginal------

log.like.marginal_alpha = function(theta, x, w, y,thetas_covariates) {
  #theta = start
  Beta = thetas_covariates[1:3] # real
  lambda = thetas_covariates[4:5] # real
  ar = thetas_covariates[6]
  ma = thetas_covariates[7]
  alpha = theta
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta
  a = NULL

  a[1] = exp(xbeta[1])
  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) + ma *
                 (-log(-log(y[t - 1])) - log(a[t - 1])))

  }


  l = sum(
    log(a) + log(alpha) +
      alpha * log(b) + (a - 1) * log(y) -
      log(1 - y ^ (a)) + (alpha - 1) * log(-log(1 - y ^ (a))) -
      (-b * log(1 - y ^ (a))) ^ (alpha)
  )
  return(l)
}

## log_like_conditional----

log.like_conditionl_covariates = function(theta, y, x, w, Etil1) {
  # x: betas covariates
  # w: lambdas covariates
  beta <- c(theta[1], theta[2], theta[3])
  lambda <- c(theta[4:5])
  ar <- theta[6]
  ma <- theta[7]
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% beta
  a <- Gb <- d <- r <- NULL
  a[1] <- exp(xbeta[1])
  Gb[1] <- pkumar(y[1],
                  a[1],
                  b = 1 / delta[1],
                  lower.tail = FALSE,
                  log.p = FALSE)
  d[1] <- dkumar(y[1], a[1], b = 1 / delta[1], log = FALSE)
  r[1]  = d[1] / Gb[1]
  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) + ma *
                 (-log(-log(y[t - 1])) - log(a[t - 1])))
    Gb[t] = pkumar(y[t],
                   a[t],
                   b = 1 / delta[t],
                   lower.tail = FALSE,
                   log.p = FALSE)
    d[t]  = dkumar(y[t], a[t], b = 1 / delta[t], log = FALSE)
    r[t]  = d[t] / Gb[t]
  }
  L = sum(Etil1 * log(Gb[1:N[k]]) + log(r[1:N[k]]),na.rm = T)
  return(L)
}

## v function-----------

V = function(theta, x, w, y) {
  #ok
  # x: betas covariates
  # w: lambdas covariates
  #theta <-  par.covariates.start
  beta <- c(theta[1], theta[2], theta[3])
  lambda <- c(theta[4:5])
  ar <- theta[6]
  ma <- theta[7]
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% beta
  a <- Gb <- NULL
  a[1] <- exp(xbeta[1])
  Gb[1] <- pkumar(y[1],
                 a[1],
                 b = 1 / delta[1],
                 lower.tail = FALSE,
                 log.p = FALSE)

  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) + ma *
                 (-log(-log(y[t - 1])) - log(a[t - 1])))
    Gb[t] = pkumar(y[t],
                   a[t],
                   b = 1 / delta[t],
                   lower.tail = FALSE,
                   log.p = FALSE)
  }
  Gb[Gb==0] <- 1
  v = sum(-log(Gb[1:N[k]]),na.rm = T)
  return(v)
}

# derivate function--------

d_vn = function(N, alpha, v) {
  #ok
  #N=N[k];v=v;alpha=alpha.real

  if (N == 1) {
    dvn = alpha * v ^ (alpha - 1)
    qnj = 1
    f_alpha_v = alpha ^ ((N)) * v ^ (((N)) * alpha - N)
    dvn <- (-1) ^ N * sum(f_alpha_v * qnj) * exp(-v ^ alpha)
    return(list(
      dvn = dvn,
      qnj = qnj,
      f_alpha_v = f_alpha_v
    ))
  }

  qnj <-  matrix(0, nrow = N, ncol = N + 1)
  qnj[, 1]  <- 0

  colnames(qnj) = paste0("j=", seq(0, N))
  rownames(qnj) = paste0("N=", seq(1, N))

  for (n in 1:N) {
    qnj[n, n + 1] <- 1
  }
  for (n in 2:N) {
    for (j_pos in 2:(N)) {
      #j = j_pos-1
      qnj[n, j_pos] <-
        qnj[n - 1, j_pos - 1] + (n - 1 - (j_pos - 1) * alpha) * qnj[n - 1, j_pos]
    }
  }

  f_alpha_v = alpha ^ ((0:N)) * v ^ (((0:N)) * alpha - N)
  dvn <- (-1) ^ N * sum(f_alpha_v * qnj[N, ]) * exp(-v ^ alpha)
  return(list(
    dvn = dvn,
    qnj = qnj,
    f_alpha_v = f_alpha_v
  ))
}



# data genetarion --------

## real par definition --------

number_of_years = c(5, 10, 15,20, 25)
#number_of_years = c(10)
N = number_of_years*12 # length of sample
betas.real = as.matrix(c(4.5, -1, -.5))
lambdas.real = as.matrix(c(-1, -3))
ar.real = .5
ma.real = .5
alpha.real = .85
theta.real <- c(betas.real,
                lambdas.real,
                ar.real,
                ma.real,
                alpha.real)
## y generation ----------

emv.results = NULL # store estimates
B = 50 # repeat B times
emv.b = matrix(0, B, length(theta.real))

for(k in 1:length(number_of_years)){

  #k=1
  j = 1:N[k]
  x1 = rep(1, N[k])
  x2 = sin(2 * pi * j / 12)
  x3 = cos(2 * pi * j / 12)
  w1 = rep(1, N[k])
  w2 = rep(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), N[k] / 12)
  x = cbind(x1, x2, x3)
  w = cbind(w1, w2)
  xbetas.real = x %*% betas.real
  wlambdas.real = w %*% lambdas.real
  delta.real = exp(wlambdas.real)
  b.real = 1 / delta.real

  for(b in 1:B){

    a.real = NULL
    y = NULL
    a.real[1] = exp(xbetas.real)[1]

    y[1] = (pweibull(
      -log(1 - runif(1)),
      shape = 1 / alpha.real,
      scale = b.real[1] ^ (alpha.real),
      lower.tail = TRUE,
      log.p = FALSE
    )) ^ (1 / a.real[1])

    for (t in 2:N[k]) {
      a.real[t] = exp(xbetas.real[t] +
                        ar.real * (-log(-log(y[t - 1])) - xbetas.real[t - 1]) +
                        ma.real * (-log(-log(y[t - 1])) - log(a.real[t - 1]))
      )
      repeat {
        u = runif(1)
        y[t] = (
          pweibull(
            -log(1 - u),
            shape = 1 / alpha.real,
            scale = b.real[t] ^ (alpha.real),
            lower.tail = TRUE,
            log.p = FALSE
          )
        ) ^ (1 / a.real[t])

        if (format(y[t], scientific = FALSE) != 1) {
          break
        }
      }
    }
    hist(y)
    plot.ts(y)

    # EM--------

    theta.start <-
      c(betas.real, lambdas.real, ar.real, ma.real, alpha.real)
    erro = 10^(-4)

    repeat {
      # betas.start <- theta_start[1:3]
      # lambdas.start <- theta_start[4:5]
      # ar.start <- theta_start[6]
      # ma.start <- theta_start[7]

      par.covariates.start <- theta.start[1:7] # covariates parameter
      alpha.start <- theta.start[8] # alpha parameter


      # Step M----

      alpha.up <-
        try(
          optim(
            par = alpha.start,
            fn = log.like.marginal_alpha,
            control = list(fnscale = -1),
            method = "L-BFGS-B",
            lower = .25,
            upper = .95,
            x = x,
            w = w,
            y = y,
            thetas_covariates = par.covariates.start
          ),
          silent = T
        )

      if(class(alpha.up)=="try-error"){
        emv <- rep("NA",length(theta.real))
        break
      }

      # Step E-----

      v <- V(theta = par.covariates.start,
             x = x,
             w = w,
             y = y)


      derivate_numerator <- d_vn(N = N[k],
                                 v = v,
                                 alpha = alpha.up$par)$dvn

      derivate_denominator <- d_vn(N = N[k] - 1,
                                   v = v,
                                   alpha = alpha.up$par)$dvn


      Esp.z <- -(derivate_numerator / derivate_denominator)
      if(is.finite(Esp.z)==F){
        emv <- rep("NA",length(theta.real))
        break
      }
      # Step M----

      par.covariates.up <-
        try(
          optim(
            par = par.covariates.start,
            fn = log.like_conditionl_covariates,
            control = list(fnscale = -1),
            method = "BFGS",
            x = x,
            w = w,
            y = y,
            Etil1 = Esp.z
          ),
          silent = T
        )

      if(class(par.covariates.up)=="try-error"){
        emv <- rep("NA",length(theta.real))
        break
      }

      theta.up <-
        c(par.covariates.up$par,
          alpha.up$par)

      crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)

      if (crit < erro) {
        emv <- theta.up
        break
      }
      else{
        theta.start <- theta.up
      }
    }
    emv.b[b,] <- emv
  }
  emv.results.aux <- data.frame(N=N[k],
                                emv.b=emv.b)
  emv.results <- rbind(emv.results,
                       emv.results.aux)
}
emv.results$N.f = factor(emv.results$N,
                         levels = N)

write.table(emv.results,"likelihood_marginal_covariate_reparametrization_delta_time_ar_ma_alpha_estimation_EM_covariates_teste1.txt",col.names = T)
emv.results = read.table("likelihood_marginal_covariate_reparametrization_delta_time_ar_ma_alpha_estimation_EM_covariates_teste1.txt")

boxplot(emv.results$emv.b.8~emv.results$N.f)
which(is.finite(emv.results$emv.b.8)==F)

emv.results$N.f <- factor(emv.results$N.f)
teste = pivot_longer(emv.results,cols = c(paste0("emv.b.",1:8)))

real.values = data.frame(name=factor(paste0("emv.b.",1:8),
                                         levels = paste0("emv.b.",1:8)),
                         value.real=theta.real)
ggplot(teste,aes(N.f,value)) +
  geom_boxplot()+
  facet_wrap(~name,scales = "free_y")+
  geom_hline(data = real.values,aes(yintercept = value.real))

