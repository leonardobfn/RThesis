# Funções-----------
#log.f.cond---------
log.f.cond <- function(theta,y,x,w,ncx,ncv){
  #theta = c(par.cov.start.marginal, alpha.start)
  # Beta = theta[1:ncx]
  # lambda = theta[c((1+ncx):(ncx+ncv))]

  alpha = theta[-c((1):(ncx+ncv))]
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
  return( log.f1  + sum(log.f.cond.))

}

#snowfall.simulation--------
snowfall.simulation = function(steps,
                               model,
                               alpha.value,
                               formula,
                               wd_sample,
                               wd,
                               wd_covariates,
                               n)
{
  # model=1
  #  steps = 650
  #  alpha.value = "alpha95"
  #print(MC)

  # setwd(wd_sample)
  # path.sample = paste0("Model_", model,"/simulations/",alpha.value,"/")
  # #file.sample = paste0( path.sample, data.labels[steps])
  # file.sample = paste0( path.sample,"data",steps,".txt")
  # data.sample = read.table(file.sample)
  # wd_covariates=wd.covariates
  #wd_sample = wd.sample
  # formula = FORMULA
  #n = 216
  setwd(wd_covariates)
  covi = read.table("covariates.txt")[1:n,]

  covi$semester <- as.factor(covi$semester)

  data = data.frame(covi, RH = rnorm(nrow(covi)))
  mf <- model.frame(Formula::Formula(formula), data = data)
  #y <- model.response(mf)
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)

  ncx <- ncol(cov_a)

  ncv <- ncol(cov_delta)
  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")

  par_real = read.table(paste0("real_par_model_", model, ".txt"))
  beta.real <- par_real[c(1:ncx), 1]
  lambda.real <-
    par_real[c((1 + ncx):(ncx + ncv)), 1] %>% as.matrix()
  alpha_ger <-
    paste0("0.", stringr::str_sub(alpha.value, 6, 7)) %>% as.numeric()

  a = exp(cov_a %*% beta.real)
  delta = exp(cov_delta %*% lambda.real)
  b = 1 / delta

  # geração dos dados inicio----------
  nn = n
  yt.y1 = matrix(0, nn, 1)

  at1 = a[1]
  bt1 = b[1]
  #u1 = runif(1)
  repeat {
    u1 = runif(1)
    #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
    y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
    yt.y1[1] <-  y1
    if (format(y1) != 1 & format(y1) != 0)
      break
  }

  for (i in 2:nn) {
    #i = 84
    at1 = a[i - 1]
    bt1 = b[i - 1]
    at = a[i]
    bt = b[i]

    # repeat{
    #   u2 = runif(1,0,1)
    #   if(abs(u2-0)>0.02 & abs(u2-1)>0.02)
    #     break
    # }

    u2 = runif(1, 0, 1)

    int.yt = try (uniroot(
      p.cond,
      interval = c(0, 1),
      u = 1 - u2,
      y.t.1 = yt.y1[i - 1],
      at1 = at1,
      at = at,
      bt1 = bt1,
      bt = bt,
      alpha = alpha_ger
    ),
    silent = T)

    if (class(int.yt) == "try-error" & i == 2) {
      repeat {
        repeat {
          u1 = runif(1)
          #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
          y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
          yt.y1[1] <-  y1
          if (format(y1) != 1 & format(y1) != 0)
            break
        }

        u2 = runif(1, 0, 1)

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 1],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }
      }
    }
    #int.yt = list()
    #int.yt$root = 0.000000e+00
    if (class(int.yt) == "try-error") {
      int.yt = list()
      int.yt$root = 0
      test2 = "try-error"
    } else{
      test2 = "no-error"
    }

    if (class(int.yt) != "try-error" &
        i == 2 &
        (format(int.yt$root) == 0 | format(int.yt$root) == 1)) {
      test = "try-error"
      while (test == "try-error") {
        repeat {
          u1 = runif(1)
          #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
          y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
          yt.y1[1] <-  y1
          if (format(y1) != 1 & format(y1) != 0)
            break
        }

        u2 = runif(1, 0, 1)

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 1],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)
        test = class(int.yt)
        if (test != "try-error") {
          if (format(int.yt$root) != 0 & format(int.yt$root) != 1)
            break
          else{
            test = "try-error"
          }
        }
        # if (class(int.yt) != "try-error") {
        #   break
        # }
      }
    }

    if ((class(int.yt) == "try-error" |
         test2 == "try-error") & i > 2) {
      repeat {
        u2 = runif(1, 0, 1)
        at1 = a[i - 2]
        bt1 = b[i - 2]
        at = a[i - 1]
        bt = b[i - 1]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 2],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)
        if (class(int.yt) == "try-error") {
          int.yt = list()
          int.yt$root <- NA

        }

        yt.y1[i - 1] <- int.yt$root

        u2 = runif(1, 0, 1)
        at1 = a[i - 1]
        bt1 = b[i - 1]
        at = a[i]
        bt = b[i]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 1],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }
      }
    }

    if (class(int.yt) != "try-error" &
        i > 2 &
        (format(int.yt$root) == 0 | format(int.yt$root) == 1)) {
      test = "try-error"
      while (test == "try-error") {
        u2 = runif(1, 0, 1)
        at1 = a[i - 2]
        bt1 = b[i - 2]
        at = a[i - 1]
        bt = b[i - 1]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 2],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) == "try-error") {
          int.yt = list()
          int.yt$root <- NA
        }

        yt.y1[i - 1] <- int.yt$root

        u2 = runif(1, 0, 1)
        at1 = a[i - 1]
        bt1 = b[i - 1]
        at = a[i]
        bt = b[i]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 1],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha_ger
        ),
        silent = T)

        test = class(int.yt)

        if (test != "try-error") {
          if (format(int.yt$root) != 0 & format(int.yt$root) != 1)
            break
          else{
            test = "try-error"
          }
        }

        # if (class(int.yt) != "try-error") {
        #   break
        # }

      }
    }

    yt.y1[i] <- int.yt$root
  }

  y = yt.y1
  setwd(wd)
  setwd(wd_sample)
  file.sample  = paste0(alpha.value, "/data", steps, ".txt")

  y.data = data.frame(
    y = y,
    alpha = rep(alpha_ger, length(y)),
    Model = rep(paste0("Model ", model), length(y))
  )
  #saveRDS(y.data, file.sample)
  write.table(y.data, file.sample)
  setwd(wd)
}




## p. cond ----------

p.cond <- function(u,y.t.1, y.t, at1, at, bt1, bt, alpha) {
  G.base.t.1 = extraDistr::pkumar(y.t.1, a = at1, b = bt1,lower.tail = F)
  G.base.t = extraDistr::pkumar(y.t, a = at, b = bt,lower.tail = F)
  LAMBDA.t.1 = -log(G.base.t.1)
  LAMBDA.t = -log(G.base.t)
  g.base.t = extraDistr::dkumar(y.t, a = at, b = bt)
  lambda.t = g.base.t / G.base.t
  return(
    1-u- (LAMBDA.t.1 ^ (1-alpha) *
            exp(LAMBDA.t.1 ^ alpha - (LAMBDA.t.1 + LAMBDA.t) ^ alpha) *
            (LAMBDA.t.1 + LAMBDA.t) ^ (alpha -1 )
          )
  )
}

## f.cond--------
f.cond <- function(y.t.1, y.t, at1, at, bt1, bt, alpha) {
  G.base.t.1 = extraDistr::pkumar(y.t.1, a = at1, b = bt1,lower.tail = F)
  G.base.t = extraDistr::pkumar(y.t, a = at, b = bt,lower.tail = F)
  LAMBDA.t.1 = -log(G.base.t.1)
  LAMBDA.t = -log(G.base.t)
  g.base.t = extraDistr::dkumar(y.t, a = at, b = bt)
  lambda.t = g.base.t / G.base.t
  return(
    lambda.t * LAMBDA.t.1 ^ (1 - alpha) *
      exp(LAMBDA.t.1 ^ alpha - (LAMBDA.t.1 + LAMBDA.t) ^ alpha) *
      (LAMBDA.t.1 + LAMBDA.t) ^ (alpha - 2) *
      (alpha * (LAMBDA.t.1 + LAMBDA.t) ^ alpha + (1 - alpha))
  )
}
## integral da f.cond-----------

int = function(l, up, ...) {
  return(integrate(f.cond,
                   lower = l,
                   upper = up, ...)$value)

}

##log.f.cond.beta ----------
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

##log.f.cond.lambda ----------
log.f.cond.lambda <- function(theta,alpha,beta,y,x,w,ncx,ncv){
  #theta = c(par.cov.start.marginal, alpha.start)
  Beta = beta
  lambda = theta
  alpha = alpha
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  #delta <- exp(wlambda)
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
##log.f.cond.lambda.alpha  ----------

log.f.cond.lambda.alpha <- function(theta,beta,y,x,w,ncx,ncv){
  #theta = c(par.cov.start.marginal, alpha.start)
  #Beta = beta
  #lambda = theta[c(1:ncv)]
  alpha = theta[-c(1:ncv)]
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  #delta <- exp(wlambda)
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
  return( log.f1  + sum(log.f.cond.))

}

##log.f.cond.alpha ----------

log.f.cond.alpha <- function(alpha,theta,y,x,w,ncx,ncv){
  #theta = c(par.cov.start.marginal, alpha.start)
  #Beta = beta;lambda = lambda
  #x = cov_a;w = cov_delta
  Beta = theta[1:ncx]
  lambda = theta[c((1+ncx):(ncx+ncv))]

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
  #return(  log.f1  + sum(  log.f.cond.)  )
  l = c(log.f1,sum(  log.f.cond. ))
  return(  sum(l[is.finite(l)==T])  )

}

##log_like_gbase----

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

##log_like_margina_alpha.ind-------

log.like.marginal.alpha.ind=function(alpha,theta,x,w,y){
  #theta = par.covariates.start.aux;x=cov_a;w=cov_delta;alpha=.5
  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)]
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  delta = exp(wlambda)
  b = c(1/delta)

  a = c(exp(xbeta))

  l = sum( log(a) + log(alpha) +
             alpha*log(b) + (a-1)*log(y) -
             log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) -
             (-b*log(1-y^(a)))^(alpha),na.rm = T)
  return(l)
}

## log_like_marginal.ind------

log.like.marginal.ind = function(theta, x, w, y) {
  #theta = c(par.covariates.start, alpha.start);x=cov_a;w=cov_delta;alpha=.5
  Beta = theta[1:ncx]
  lambda = theta[c((1+ncx):(ncx+ncv))]
  alpha = theta[-c((1):(ncx+ncv))]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta

  a = exp(xbeta)

  l = sum(
    log(a) + log(alpha) +
      alpha * log(b) + (a - 1) * log(y) -
      log(1 - y ^ (a)) + (alpha - 1) * log(-log(1 - y ^ (a))) -
      (-b * log(1 - y ^ (a))) ^ (alpha),
    na.rm = T)
  return(l)
}

## log_like_conditional----

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

## v function-----------

V = function(theta, x, w, y,ncx,ncv) {
  #ok
  # x: betas covariates
  # w: lambdas covariates
  #theta <-  par.covariates.start
  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)]# real
  wlambda <- w %*% lambda
  delta <- exp(wlambda)
  xbeta <- x %*% Beta
  a = exp(xbeta)
  Gb = pkumar(y,
              a,
              b = 1 / delta,
              lower.tail = FALSE,
              log.p = FALSE)
  Gb[Gb == 0] <- 1
  v = sum(-log(Gb), na.rm = T)
  return(v)
}


# derivate function--------

d_vn = function(N,alpha,v){
  #N=1;v=2;alpha=.5

  if(N==1){
    dvn = alpha*v^(alpha-1)
    qnj=1
    f_alpha_v = alpha^((N))*v^(((N))*alpha-N)
    dvn <- (-1)^N*sum(f_alpha_v*qnj)*exp(-v^alpha)
    return(list(dvn=dvn,qnj=qnj,f_alpha_v=f_alpha_v))
  }

  qnj <-  matrix(0, nrow = N, ncol = N+1)
  qnj[,1]  <- 0

  colnames(qnj) = paste0("j=",seq(0,N))
  rownames(qnj) = paste0("N=",seq(1,N))

  for(n in 1:N){
    qnj[n,n+1] <- 1
  }

  for(n in 2:N){

    for(j_pos in 2:(N)){

      #j = j_pos-1
      qnj[n,j_pos] <- qnj[n-1,j_pos-1] + (n-1-(j_pos-1)*alpha)*qnj[n-1,j_pos]
    }
  }

  f_alpha_v = alpha^((0:N))*v^(((0:N))*alpha-N)
  dvn <- (-1)^N*sum(f_alpha_v*qnj[N,])*exp(-v^alpha)
  return(list(dvn=dvn,qnj=qnj,f_alpha_v=f_alpha_v))
}

snowfall.estimates_Method_1 = function(steps,
                                       model,
                                       alpha.value,
                                       erro = 10 ^ (-4),
                                       formula,
                                       wd_sample,
                                       wd_results_emv,
                                       wd_results_hessian,
                                       wd,
                                       wd_covariates,
                                       n
)
{
  # model=1
  #  steps = 650
  #
  #print(MC)

  # setwd(wd_sample)
  # path.sample = paste0("Model_", model,"/simulations/",alpha.value,"/")
  # #file.sample = paste0( path.sample, data.labels[steps])
  # file.sample = paste0( path.sample,"data",steps,".txt")
  # data.sample = read.table(file.sample)
  # erro = 10^(-4)
  # wd_covariates=wd.covariates
  # wd_results_emv = wd.results_emv
  # wd_results_hessian = wd.results_emv
  # wd_sample = wd.sample
  # wd= wd.
  # model = 1
  # alpha.value = "alpha50"
  # formula = FORMULA
  # steps = 1025
  # n = n.

  setwd(wd)
  setwd(wd_sample)
  file.sample = paste0(alpha.value,"/data",steps,".txt")
  data.sample = read.table(file.sample)
  setwd(wd)
  setwd(wd_covariates)
  covi = read.table("covariates.txt")[1:n,]

  covi$semester <- as.factor(covi$semester)

  data = data.frame(covi, RH = data.sample$y)
  mf <- model.frame(Formula::Formula(formula), data = data)
  y <- model.response(mf)
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)

  ncx <- ncol(cov_a)

  ncv <- ncol(cov_delta)
  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")
  par_real = read.table(paste0("real_par_model_", model, ".txt"))


  start_aux_betas = try(coef(lm(-log(-log(RH)) ~ sent + cost, data = data)))

  if (class(start_aux_betas) != "try-error") {
    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux_betas, rep(-.1, ncv)),
        fn = log.kumar.gbase,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,6),0.25),
        # upper = c(rep(Inf,6),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        ncx=ncx,
        ncv=ncv
      ),
      silent = T)

    if (class(par.cov.Gbase) != "try-error" & par.cov.Gbase$value==0) {
      theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
    }

    if (class(par.cov.Gbase) == "try-error") {
      theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
    }else{
      theta.start = par.cov.Gbase$par
    }
  } else{
    theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
  }


  theta.start = c(par.cov.Gbase$par,0.5)
  value.start = theta.start

  repeat {

    beta <- theta.start[c(1:ncx)]
    lambda <- theta.start[c((1 + ncx):(ncx + ncv))]
    alpha <- theta.start[-c((1):(ncx + ncv))]

    emv.beta <-
      try(optim(
        par = beta,
        fn = log.f.cond.beta,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,ncx+ncv),0.25),
        # upper = c(rep(Inf,ncx+ncv),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        lambda = lambda,
        alpha = alpha,
        ncx = ncx,
        ncv = ncv,
        hessian = T
      ),
      silent = T)

    emv.lambda <-
      try(optim(
        par = lambda,
        fn = log.f.cond.lambda,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,ncx+ncv),0.25),
        # upper = c(rep(Inf,ncx+ncv),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        beta = emv.beta$par,
        alpha = alpha,
        hessian = T
      ),
      silent = T)


    emv.alpha <-
      try(optim(
        par = alpha,
        fn = log.f.cond.alpha,
        control = list(fnscale = -1),
        #method = "BFGS",
        method = "L-BFGS-B",
        lower = c(0.1),
        upper = c(.99),
        x = cov_a,
        w = cov_delta,
        y = y,
        ncx = ncx,
        ncv = ncv,
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    # if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
    #   estimates.aux = list()
    #   estimates.aux$par = rep(NA, ncx + ncv + 1)
    #   break
    # }
    #
    # if (class(emv.alpha) == "try-error") {
    #   estimates.aux = list()
    #   estimates.aux$par = rep(NA, ncx + ncv + 1)
    #   break
    # } else{
    #   theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    # }

    if (class(emv.alpha) == "try-error") {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)


    crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)

    if (crit < erro) {
      estimates.aux = list()
      estimates.aux$par = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  if (length(which(is.na(estimates.aux$par) == T)) == 0) {
    emv = list()
    emv <- round(estimates.aux$par, 4)
    emv.BETA = list()
    emv.BETA$hessian = round(emv.beta$hessian, 4)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = round(emv.lambda$hessian, 4)
    emv.ALPHA = list()
    emv.ALPHA$hessian = round(emv.alpha$hessian, 4)
  } else{
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = matrix(0, ncx, ncx)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = matrix(0, ncv, ncv)
    emv.ALPHA = list()
    emv.ALPHA$hessian = 0
  }
  setwd(wd)
  setwd(wd_results_emv)
  #steps=1
  file.estimates.emv  = paste0("estimates_",alpha.value,"_n",n, ".txt")

  estimates = data.frame(steps = rep(steps,ncx+ncv+1),
                         alpha = rep(data.sample$alpha[1],ncx+ncv+1),
                         par_name = par_names,
                         par_real = c(par_real[,1],data.sample$alpha[1]),
                         emv,
                         Method = rep("Method 1",ncx+ncv+1),
                         n = rep(n,ncx+ncv+1),
                         Model = rep(paste0("Model ",model),ncx+ncv+1),
                         value.start = round( value.start,4)
  )

  write.table(
    estimates,
    file.estimates.emv,
    col.names = F,
    row.names = F,
    append = T,
    quote = T,
  )
  setwd(wd)

  setwd(wd.results_hessian)

  file.hessian.alpha = paste0("hessian_",alpha.value,"_n",n, ".txt")
  file.hessian.covBeta = paste0("hessian_Beta_",alpha.value,"_n",n,".txt")
  file.hessian.covGamma = paste0("hessian_Gamma_",alpha.value,"_n",n,".txt")

  hessianAlpha = cbind(steps, emv.ALPHA$hessian)
  hessianCovBeta = cbind(steps, emv.BETA$hessian)
  hessianCovLambda = cbind(steps, emv.LAMBDA$hessian)

  write.table(
    hessianAlpha,
    file.hessian.alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  write.table(
    hessianCovBeta,
    file.hessian.covBeta,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  write.table(
    hessianCovLambda,
    file.hessian.covGamma,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  setwd(wd)
}

snowfall.estimates_Method_2 = function(steps,
                                       model,
                                       alpha.value,
                                       erro = 10 ^ (-4),
                                       formula,
                                       wd_sample,
                                       wd_results_emv,
                                       wd_results_hessian,
                                       wd,
                                       wd_covariates,
                                       n)
{
  # model=1
  #  steps = 650
  #
  #print(MC)


  # erro = 10^(-4)
  # wd_covariates=wd.covariates
  # wd_results_emv = wd.results_emv
  # wd_results_hessian = wd.results_emv
  # wd_sample = wd.sample
  # wd= wd.
  # model = 1
  # alpha.value = "alpha50"
  # formula = FORMULA
  # steps = 1025
  # n = n.

  setwd(wd)
  setwd(wd_sample)
  file.sample = paste0(alpha.value, "/data", steps, ".txt")
  data.sample = read.table(file.sample)
  setwd(wd)
  setwd(wd_covariates)
  covi = read.table("covariates.txt")[1:n, ]

  covi$semester <- as.factor(covi$semester)

  data = data.frame(covi, RH = data.sample$y)
  mf <- model.frame(Formula::Formula(formula), data = data)
  y <- model.response(mf)
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)

  ncx <- ncol(cov_a)

  ncv <- ncol(cov_delta)
  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")
  par_real = read.table(paste0("real_par_model_", model, ".txt"))


  start_aux_betas = try(coef(lm(-log(-log(RH)) ~ sent + cost, data = data)))

  if (class(start_aux_betas) != "try-error") {
    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux_betas, rep(-.1, ncv)),
        fn = log.kumar.gbase,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,6),0.25),
        # upper = c(rep(Inf,6),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        ncx = ncx,
        ncv = ncv
      ),
      silent = T)

    if (class(par.cov.Gbase) != "try-error" &
        par.cov.Gbase$value == 0) {
      theta.start = c(0.80 , 0.20, -0.02 , rep(-.1, ncv))
    }

    if (class(par.cov.Gbase) == "try-error") {
      theta.start = c(0.80 , 0.20, -0.02 , rep(-.1, ncv))
    } else{
      theta.start = par.cov.Gbase$par
    }
  } else{
    theta.start = c(0.80 , 0.20, -0.02 , rep(-.1, ncv))
  }


  theta.start = c(par.cov.Gbase$par, 0.5)


  repeat {
    beta <- theta.start[c(1:ncx)]
    lambda <- theta.start[c((1 + ncx):(ncx + ncv))]
    alpha <- theta.start[-c((1):(ncx + ncv))]

    emv.beta <-
      try(optim(
        par = beta,
        fn = log.f.cond.beta,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,ncx+ncv),0.25),
        # upper = c(rep(Inf,ncx+ncv),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        lambda = lambda,
        alpha = alpha,
        ncx = ncx,
        ncv = ncv,
        hessian = T
      ),
      silent = T)

    emv.lambda <-
      try(optim(
        par = lambda,
        fn = log.f.cond.lambda,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,ncx+ncv),0.25),
        # upper = c(rep(Inf,ncx+ncv),0.95),
        x = cov_a,
        w = cov_delta,
        y = y,
        beta = emv.beta$par,
        alpha = alpha,
        hessian = T
      ),
      silent = T)


    emv.alpha <-
      try(optim(
        par = alpha,
        fn = log.f.cond.alpha,
        control = list(fnscale = -1),
        #method = "BFGS",
        method = "L-BFGS-B",
        lower = c(0.1),
        upper = c(.99),
        x = cov_a,
        w = cov_delta,
        y = y,
        ncx = ncx,
        ncv = ncv,
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)


    if (class(emv.alpha) == "try-error") {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)

    crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)

    if (crit < erro) {
      estimates.aux = list()
      estimates.aux$par = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  alpha.up <- estimates.aux$par[-c(1:(ncx + ncv))]
  par.covariates.start <- estimates.aux$par[c(1:(ncx + ncv))]
  value.start = c(estimates.aux$par[c(1:(ncx + ncv))],0.5)

  repeat {
    if (length(which(is.na(estimates.aux$par) == T)) > 0) {
      par.covariates.up = list()
      par.covariates.up$hessian = matrix(0,ncx+ncv,ncx+ncv)
      emv.alpha = list()
      emv.alpha$hessian = c(0)
      emv <- rep(NA, ncx + ncv + 1)
      break
    }
    #       # Step E-----
    v <- V(
      theta = par.covariates.start,
      x = cov_a,
      w = cov_delta,
      y = y,
      ncx=ncx,
      ncv=ncv
    )


    derivate_numerator <- try(d_vn(N = length(y) + 1,
                                   v = v,
                                   alpha = alpha.up)$dvn)

    if (class(derivate_numerator) == "try-error" |
        is.finite(derivate_numerator) == F |
        derivate_numerator == 0) {
      par.covariates.up = list()
      par.covariates.up$hessian = matrix(0,ncx+ncv,ncx+ncv)
      emv.alpha = list()
      emv.alpha$hessian = c(0)
      emv <- rep(NA, ncx + ncv + 1)
      break
    }

    derivate_denominator <- try(d_vn(N = length(y),
                                     v = v,
                                     alpha = alpha.up)$dvn)

    if (class(derivate_denominator) == "try-error" |
        is.finite(derivate_denominator) == F |
        derivate_denominator == 0) {
      par.covariates.up = list()
      par.covariates.up$hessian = matrix(0,ncx+ncv,ncx+ncv)
      emv.alpha = list()
      emv.alpha$hessian = c(0)
      emv <- rep(NA, ncx + ncv + 1)
      break
    }

    Esp.z <- -(derivate_numerator / derivate_denominator)

    # Step M----

    par.covariates.up <-
      try(optim(
        par = par.covariates.start,
        fn = log.like_conditionl_covariates_EM,
        control = list(fnscale = -1),
        method = "BFGS",
        x = cov_a,
        w = cov_delta,
        y = y,
        Etil1 = Esp.z,
        hessian = T,
        ncx = ncx,
        ncv=ncv
      ),
      silent = T)

    if (class(par.covariates.up) == "try-error") {
      par.covariates.up = list()
      par.covariates.up$hessian = matrix(0,ncx+ncv,ncx+ncv)
      emv.alpha = list()
      emv.alpha$hessian = c(0)
      emv <- rep(NA, ncx + ncv + 1)
      break
    }

    if (class(par.covariates.up) != "try-error" &
        par.covariates.up$value == 0) {
      par.covariates.up = list()
      par.covariates.up$hessian = matrix(0,ncx+ncv,ncx+ncv)
      emv.alpha = list()
      emv.alpha$hessian = c(0)
      emv <- rep(NA, ncx + ncv + 1)
      break
    }

    crit <-
      sum(((par.covariates.up$par - par.covariates.start) / par.covariates.start
      ) ^ 2)
    if (crit < erro) {
      emv <- c(par.covariates.up$par, alpha.up)
      break
    }
    else{
      par.covariates.start <- par.covariates.up$par
    }
  }




  setwd(wd)
  setwd(wd_results_emv)
  #steps=1
  file.estimates.emv  = paste0("estimates_", alpha.value, "_n", n, ".txt")

  estimates = data.frame(
    steps = rep(steps, ncx + ncv + 1),
    alpha = rep(data.sample$alpha[1], ncx + ncv + 1),
    par_name = par_names,
    par_real = c(par_real[, 1], data.sample$alpha[1]),
    emv,
    Method = rep("Method 2", ncx + ncv + 1),
    n = rep(n, ncx + ncv + 1),
    Model = rep(paste0("Model ", model), ncx + ncv +
                  1),
    value.start = round(value.start, 4)
  )

  write.table(
    estimates,
    file.estimates.emv,
    col.names = F,
    row.names = F,
    append = T,
    quote = T,
  )
  setwd(wd)

  setwd(wd.results_hessian)

  file.hessian.alpha = paste0("hessian_", alpha.value, "_n", n, ".txt")
  file.hessian.cov = paste0("hessian_cov_", alpha.value, "_n", n, ".txt")

  hessianAlpha = cbind(steps, emv.alpha$hessian)
  hessianCov = cbind(steps, par.covariates.up$hessian)

  write.table(
    hessianAlpha,
    file.hessian.alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  write.table(
    hessianCov,
    file.hessian.cov,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  setwd(wd)
}


