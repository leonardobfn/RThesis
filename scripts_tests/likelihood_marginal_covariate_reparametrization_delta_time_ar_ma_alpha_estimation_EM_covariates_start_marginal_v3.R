rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

##log_like_gbase----_-

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

log.like.marginal = function(theta, x, w, y) {
  #theta = start
  Beta = theta[1:3] # real
  lambda = theta[4:5] # real
  ar = theta[6]
  ma = theta[7]
  alpha = theta[8]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta
  a = phi.median = NULL

  a[1] = exp(xbeta[1])
  median <- (1- exp(-  b[1]*(-log(.5))^1/alpha  ))^(1/a[1])
  phi.median [1] <- -log(-log(median))

  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) + ma *
                 (-log(-log(y[t - 1])) - phi.median[t-1]  ))

    median <- (1- exp(-  b[t]*(-log(.5))^1/alpha  ))^(1/a[t])
    phi.median [t] <- -log(-log(median))
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

log.like_conditionl_covariates = function(theta, y, x, w, Etil1,alpha) {
  # x: betas covariates
  # w: lambdas covariates
  beta <- c(theta[1], theta[2], theta[3])
  lambda <- c(theta[4:5])
  ar <- theta[6]
  ma <- theta[7]
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  b=1/delta
  xbeta <- x %*% beta
  a <- Gb <- d <- r <- phi.median <- NULL
  a[1] <- exp(xbeta[1])
  Gb[1] <- pkumar(y[1],
                  a[1],
                  b = 1 / delta[1],
                  lower.tail = FALSE,
                  log.p = FALSE)
  d[1] <- dkumar(y[1], a[1], b = 1 / delta[1], log = FALSE)
  r[1]  = d[1] / Gb[1]
  #median <- (1- exp(-  b[1]*(-log(.5))^1/alpha  ))^(1/a[1])
  #phi.median [1] <- -log(-log(median))
  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] + ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) + ma *
                 (-log(-log(y[t - 1])) - log(a[t-1] ) )
    )
    #median <- (1- exp(-  b[t]*(-log(.5))^1/alpha  ))^(1/a[t])
    #phi.median [t] <- -log(-log(median))
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
  b=1/delta
  xbeta <- x %*% beta
  a <- Gb <- phi.mediana <- NULL
  a[1] <- exp(xbeta[1])
  Gb[1] <- pkumar(y[1],
                  a[1],
                  b = 1 / delta[1],
                  lower.tail = FALSE,
                  log.p = FALSE)
  # phi.mediana[1] <- # median gbase
  #   log(a[1]) - log(-log(1 - .5 ^ (1 / b[1])))
  for (t in 2:N[k]) {
    a[t] = exp(xbeta[t] +
                 ar * (-log(-log(y[t - 1])) - xbeta[t - 1]) +
                 ma *(-log(-log(y[t - 1])) - log(a[t - 1]) ))
                #ma*(-log(-log(y[t-1])) - phi.mediana[t-1]))

    #phi.mediana[t] <- log(a[t]) - log(-log(1 - .5 ^ (1 / b[t])))
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

#number_of_years = c(5, 10, 15,20, 25)
number_of_years = c(5,10)
N = number_of_years*12 # length of sample
betas.real = as.matrix(c(1, -1, -.5))
lambdas.real = as.matrix(c(.5, -1))
ar.real = .7
ma.real = -.5
alpha.real = .75
theta.real <- c(betas.real,
                lambdas.real,
                ar.real,
                ma.real,
                alpha.real)
## y generation ----------

emv.results = NULL # store estimates
B = 100 # repeat B times
emv.b = matrix(0, B, length(theta.real))

for(k in 1:length(number_of_years)){
#k=2
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

    a.real = phi.median.real=NULL
    y = NULL
    a.real[1] = exp(xbetas.real)[1]
    median.real <- (1- exp(-  b.real[1]*(-log(.5))^1/alpha.real  ))^(1/a.real[1])
    phi.median.real[1] <- -log(-log(median.real))
    z <- stabledist::rstable(N[k],alpha = alpha.real,
                             beta = 1,
                             gamma = (cos(pi*alpha.real/2))^(1/alpha.real),
                             delta=0,
                             pm=1)

    y[1] = extraDistr::rkumar(1,a=a.real[1],b= z[1]*b.real[1])

    for (t in 2:N[k]) {

      a.real[t] = exp(xbetas.real[t] +
                        ar.real * (-log(-log(y[t - 1])) - xbetas.real[t - 1]) +
                        ma.real * (-log(-log(y[t - 1])) - phi.median.real [t-1])
      )
      median.real <- (1- exp(-  b.real[t]*(-log(.5))^1/alpha.real  ))^(1/a.real[t])
      phi.median.real[t] <- -log(-log(median.real))
      repeat {
        y[t] =
          extraDistr::rkumar(
            1,
            a = a.real[t],
            b = z[t]*b.real[t]
          )


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
    # theta.start <-
    #   c(c(1,-1,-1), c(-1,-1),.5, -.5, .5)
    erro = 10^(-4)

    repeat {
      # betas.start <- theta_start[1:3]
      # lambdas.start <- theta_start[4:5]
      # ar.start <- theta_start[6]
      # ma.start <- theta_start[7]

      par.covariates.start <- theta.start[1:7] # covariates parameter
      alpha.start <- theta.start[8] # alpha parameter


      # Step M----

      estimates.aux <-
        try(
          optim(
            par = c(par.covariates.start,alpha.start),
            fn = log.like.marginal,
            control = list(fnscale = -1),
            method = "BFGS",
            #method = "L-BFGS-B",
            #lower = c(rep(-Inf,7),0.25),
            #upper = c(rep(Inf,7),0.95),
            x = x,
            w = w,
            y = y
          ),
          silent = T
        )

      if(class(estimates.aux)=="try-error"){

        emv <- rep("NA",length(theta.real))
        break
      }

      if(estimates.aux$par[8]>=1){

        emv <- rep("NA",length(theta.real))
        break
      }


      alpha.up <- estimates.aux$par[8]
      par.covariates.start <- estimates.aux$par[-8]
      # Step E-----

      v <- V(theta = par.covariates.start,
             x = x,
             w = w,
             y = y)


      derivate_numerator <- d_vn(N = N[k],
                                 v = v,
                                 alpha = alpha.up)$dvn

      derivate_denominator <- d_vn(N = N[k] - 1,
                                   v = v,
                                   alpha = alpha.up)$dvn


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
            Etil1 = Esp.z,
            alpha=alpha.up
          ),
          silent = T
        )

      if(class(par.covariates.up)=="try-error"){
        emv <- rep("NA",length(theta.real))
        break
      }

      theta.up <-
        c(par.covariates.up$par,
          alpha.up)

      crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)
      cat(c(crit,b,k),"\r")
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

write.table(emv.results,"likelihood_marginal_covariate_reparametrization_delta_time_ar_ma_alpha_estimation_EM_covariates_start_marginal_v3_teste3.txt",col.names = T)
#Analysis-------
rm(list = ls())
number_of_years = c(5, 10, 15,20, 25)
#number_of_years = c(15)
N = number_of_years*12 # length of sample
betas.real = as.matrix(c(1, -1, -.5))
lambdas.real = as.matrix(c(.5, -1))
ar.real = .7
ma.real = -.5
alpha.real = .75
theta.real <- c(betas.real,
                lambdas.real,
                ar.real,
                ma.real,
                alpha.real)
emv.results = read.table("likelihood_marginal_covariate_reparametrization_delta_time_ar_ma_alpha_estimation_EM_covariates_start_marginal_v3_teste3.txt")
colnames(emv.results)[-c(1,10)] <- c(
  paste0("beta",0:2),
  paste0("lambda",0:1),
  "ar",
  "ma",
  "alpha"
)

names.par = factor(colnames(emv.results) [c(-1,-10)],colnames(emv.results) [c(-1,-10)])

boxplot(emv.results$alpha~emv.results$N)
which(is.finite(emv.results$alpha)==F)

emv.results$N.f <- factor(emv.results$N)
teste = pivot_longer(emv.results,
                     cols = names.par) %>%
  mutate(name=factor(name,levels =names.par))

estimates = teste %>% select(-N) %>%
  group_by(N.f,name) %>%
  summarise(mean = mean(value,na.rm=T),
            md = median(value,na.rm=T),
            q.025 = quantile(value,c(0.025),na.rm=T),
            q.975 = quantile(value,c(0.975),na.rm=T),
            sdd = sd(value,na.rm = T),
            NAA = length(which(is.na(value)==T))) %>% data.frame()

estimates$real.values <- rep(theta.real,2)

real.values = data.frame(name=names.par,
                         value.real=theta.real
)

ggplot(teste,aes(N.f,value)) +
  geom_boxplot()+
  facet_wrap(~name,scales = "free_y")+
  geom_hline(data = real.values,aes(yintercept = value.real))

