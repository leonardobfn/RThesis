# Funções-----------
#log.f.cond---------
log.f.cond <- function(theta,y,x,w){
  #theta = c(par.cov.start.marginal, alpha.start)
  Beta = theta[1:ncx]
  lambda = theta[c((1+ncx):(ncx+ncv))]

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

snowfall.simulation <- function(nn, a, b, alpha, model, MC) {
  #id=1
  # alpha = 0.5
  # alpha = 0.8

  alpha.value <- switch (
    as.character(alpha),
    "0.35" = "alpha35",
    "0.5" = "alpha50",
    "0.65" = "alpha65",
    "0.8" = "alpha80",
    "0.95" = "alpha95"
  )

  path.sample <-
    paste0("Data_simulation/Model_",
           model,
           "/simulations/",
           alpha.value,
           "/data",
           MC,
           ".txt")

  yt.y1 = matrix(0, nn, 1)

  at1 = a[1]
  bt1 = b[1]
  u1 = runif(1)
  G = (1 - .8 ^ (1 / b[1])) ^ a[1]
  y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha))) ^ (1 / at1)
  yt.y1[1] <-  y1

  y.aux <- data.frame(y.sim = yt.y1[1],
                      t = 1,
                      alpha = alpha)

  write.table(
    y.aux,
    path.sample,
    sep = " ",
    append = T,
    quote = T,
    row.names = F,
    col.names = F
  )

  for (i in 2:nn) {
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
      alpha = alpha
    ),
    silent = T)

    if (class(int.yt) == "try-error" & i == 2) {
      repeat {
        u1 = runif(1)
        y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha))) ^ (1 / at1)
        yt.y1[1] <-  y1
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
          alpha = alpha
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }
      }
    }

    if (class(int.yt) == "try-error" & i>2) {
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
          alpha = alpha
        ),
        silent = T)

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
          alpha = alpha
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }

      }
    }

    yt.y1[i] <- int.yt$root

    y.aux <- data.frame(y.sim = yt.y1[i],
                        t = i,
                        alpha = alpha)

    write.table(
      y.aux,
      path.sample,
      sep = " ",
      append = T,
      quote = T,
      row.names = F,
      col.names = F
    )
  }

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
log.f.cond.beta <- function(theta,alpha,lambda,y,x,w){
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
  return( log.f1  + sum(log.f.cond.))

}

##log.f.cond.lambda ----------
log.f.cond.lambda <- function(theta,alpha,beta,y,x,w){
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
  return( log.f1  + sum(log.f.cond.))

}
##log.f.cond.lambda.alpha  ----------

log.f.cond.lambda.alpha <- function(theta,beta,y,x,w){
  #theta = c(par.cov.start.marginal, alpha.start)
  Beta = beta
  lambda = theta[c(1:ncv)]
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

log.f.cond.alpha <- function(alpha,theta,y,x,w){
  #theta = c(par.cov.start.marginal, alpha.start)
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
  return( log.f1  + sum(log.f.cond.))

}

##log_like_gbase----

log.kumar.gbase = function(theta, x, w, y) {
  #theta = par.covariates.start;x=cov_a;w=cov_delta;alpha=.5

  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1/delta
  a = exp(xbeta)



  l = sum(extraDistr::dkumar(
    x = y,
    a = a,
    b = b,
    log = T
  ),na.rm = T)
  return(l)
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

log.like_conditionl_covariates_EM = function(theta, y, x, w, Etil1) {
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

  L = sum(Etil1 * log(Gb) + log(r), na.rm = T)

  return(L)
}

## v function-----------

V = function(theta, x, w, y) {
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

d_vn = function(N, alpha, v) {
  #ok
  #N=TT;v=v;alpha=0.4753165

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

snowfall.estimates_Method_1 = function(steps, model, alpha,erro = 10 ^ (-4)){
  # model=1
  #  steps = 1
  #   alpha = 0.35
  #  print(MC)

  alpha.value <- switch (
    as.character(alpha),
    "0.35" = "alpha35",
    "0.5" = "alpha50",
    "0.65" = "alpha65",
    "0.8" = "alpha80",
    "0.95" = "alpha95"
  )
  path.sample <- paste0("Data_simulation/Model_",
                        model,
                        "/simulations/",
                        alpha.value,
                        "/data",
                        steps,
                        ".txt")


  s = read.table(path.sample)

  covi = read.table("Data_simulation/covariates.txt")
  covi$semester <- as.factor(covi$semester)
  data = data.frame(covi, RH = s$V1)
  formula <- RH ~ sent + cost | semester
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
  start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost, data = data))

  start_aux <-
    c(start_aux_betas, -.1, -.1)

  par.cov.Gbase <-
    try(optim(
      par =  c(start_aux),
      fn = log.kumar.gbase,
      control = list(fnscale = -1),
      method = "BFGS",
      # method = "L-BFGS-B",
      # lower = c(rep(-Inf,6),0.25),
      # upper = c(rep(Inf,6),0.95),
      x = cov_a,
      w = cov_delta,
      y = y
    ),
    silent = T)

  if (class(par.cov.Gbase) == "try-error") {
    start_aux_lambdas = try(betareg::betareg(
      formula = formula,
      link = "loglog",
      link.phi = "log",
      data = data
    )$coefficients$precision,silent = T)
    if(class(start_aux_lambdas)=="try-error"){
      start_aux = c(start_aux_betas,-2,-2)
    }else{
      start_aux <-
        c(start_aux_betas, start_aux_lambdas)
    }

    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux),
        fn = log.kumar.gbase,
        control = list(fnscale = -1),
        method = "BFGS",
        # method = "L-BFGS-B",
        # lower = c(rep(-Inf,6),0.25),
        # upper = c(rep(Inf,6),0.95),
        x = cov_a,
        w = cov_delta,
        y = y
      ),
      silent = T)

    if (class(par.cov.Gbase) == "try-error") {
      par.cov.Gbase = list()
      par.cov.Gbase$par <- rep(NA,ncx+ncv+1)
    }

  }

  # estimates.aux <-
  #   try(optim(
  #     par = c(par.cov.Gbase$par,.5),
  #     fn = log.f.cond,
  #     control = list(fnscale = -1),
  #     method = "BFGS",
  #     # method = "L-BFGS-B",
  #     #  lower = c(rep(-Inf,ncx+ncv),0.01),
  #     #  upper = c(rep(Inf,ncx+ncv),1),
  #     x = cov_a,
  #     w = cov_delta,
  #     y = y
  #   ),
  #   silent = T)
  #theta.start = c(start_aux,.5)

  theta.start = c(par.cov.Gbase$par, .5)


  beta <- theta.start[c(1:ncx)]
  lambda <- theta.start[c((1 + ncx):(ncx + ncv))]
  alpha <- theta.start[-c((1):(ncx + ncv))]

  emv.alpha <-
    try(optim(
      par = alpha,
      fn = log.f.cond.alpha,
      control = list(fnscale = -1),
      #method = "BFGS",
      method = "L-BFGS-B",
      lower = c(0.1),
      upper = c(1),
      x = cov_a,
      w = cov_delta,
      y = y,
      theta = c(beta, lambda),
      hessian = T
    ),
    silent = T)

  if(class(emv.alpha)=="try-error"|emv.alpha$value==0){
    estimates.aux = list()
    estimates.aux$par = rep(NA,ncx+ncv+1)
  }else{
    estimates.aux = list()
    estimates.aux$par = c(beta,lambda,emv.alpha$par)
  }


  alpha.up <- estimates.aux$par[-c(1:(ncx + ncv))]
  par.covariates.start <- estimates.aux$par[c(1:(ncx + ncv))]
  #print(estimates.aux)
  repeat {

    if(length(which(is.na(estimates.aux$par)==T))>0){
      par.covariates.up = list()
      par.covariates.up$value = 0
      emv <- rep(NA,ncx+ncv+1)
      break
    }
    #       # Step E-----
    v <- try(V(
      theta = par.covariates.start,
      x = cov_a,
      w = cov_delta,
      y = y
    ),silent = T)

    if(class(v)=="try-error"){
      v = NA
    }


    derivate_numerator <- try(d_vn(N = length(y) + 1,
                                   v = v,
                                   alpha = alpha.up)$dvn)

    if(class(derivate_numerator)=="try-error" | is.finite(derivate_numerator)==F|
       derivate_numerator==0){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
      break
    }

    derivate_denominator <- try(d_vn(N = length(y),
                                     v = v,
                                     alpha = alpha.up)$dvn)

    if(class(derivate_denominator)=="try-error" | is.finite(derivate_denominator)==F|
       derivate_denominator==0){

      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
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
        Etil1 =  Esp.z,
        hessian = T
      ),
      silent = T)

    if(class(par.covariates.up)=="try-error"){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
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

  if((class(par.covariates.up)!="try-error" & par.covariates.up$value ==0))
  {
    emv <- rep(NA,ncx+ncv+1)
  }


  path_estimates = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_1/estimates/",
    "estimates",
    alpha.value,
    ".txt"
  )

  estimates = data.frame(steps, emv, par_names)
  write.table(
    estimates,
    path_estimates,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  # if(length(which(is.na(emv)==T)==0)){
  # colnames(par.covariates.up$hessian) = par_names[-6]
  # rownames(par.covariates.up$hessian) = par_names[-6]
  # }

  path_hessian_cov = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_1/hessian/",
    "hessianCov_",
    alpha.value,
    ".txt"
  )
  if(length(which(is.na(emv)==T))>0){
    par.covariates.up = list()
    par.covariates.up$hessian <- matrix(0,ncx+ncv,ncx+ncv)
    hessianCov = cbind(steps, par.covariates.up$hessian)
  }else{
    hessianCov = cbind(steps, par.covariates.up$hessian)
  }


  write.table(
    hessianCov,
    path_hessian_cov,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  path_hessian_alpha = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_1/hessian/",
    "hessian_",
    alpha.value,
    ".txt"
  )
  if(length(which(is.na(emv)==T))>0){
    emv.alpha = list()
    emv.alpha$hessian <- matrix(0,1,1)
    hessianAlpha = cbind(steps, emv.alpha$hessian)
  }else{
    hessianAlpha = cbind(steps, emv.alpha$hessian)
  }

  write.table(
    hessianAlpha,
    path_hessian_alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )
}



