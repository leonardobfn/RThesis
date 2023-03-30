rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções ---------

## f.cond ---------

f.cond <- function(theta,y,x,w){
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
##log_like_gbase-------

log.like.kumar = function(theta, x, w, y) {
  #theta = par.covariates.start;x=cov_a;w=cov_delta;alpha=.5

  Beta = theta[1:ncx]
  lambda = theta[-c(1:ncx)]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta
  a = exp(xbeta)


  l = sum(extraDistr::dkumar(
    x = y,
    a = a,
    b = b,
    log = T
  ),na.rm = T)
  return(l)
}
##log_like_margina_alpha-------
log.like.marginal.alpha=function(alpha,theta,x,w,y){
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

## log_like_marginal------

log.like.marginal = function(theta, x, w, y) {
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

log.like_conditionl_covariates = function(theta, y, x, w, Etil1) {
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
number_of_years = c(14)

N = number_of_years * 12 # length of sample
beta0.real = 0.44875154
betas.real = as.matrix(c(beta0.real,0.2191953 ,-0.1158410))
lambda0.real = -3.0767779
lambdas.real = as.matrix(c(lambda0.real, -0.5170042))
alpha.real = .85
theta.real <- c(betas.real,
                lambdas.real,
                alpha.real)
par.names = c(
  paste0("beta",0:2),
  paste0("lambda",0:1),
  "alpha"
)
## y generation ----------

emv.results = NULL # store estimates
B = 50 # repeat B times
emv.b = matrix(0, B, length(theta.real))
tic <- tictoc::tic()
for(h in 1:100){
z = stabledist::rstable(
  1,
  alpha.real,
  beta = 1,
  gamma = (cos((pi*alpha.real)/2))^(1/alpha.real),
  delta = 0,
  pm = 1
)

for (k in 1:length(number_of_years)) {
  #  k=1
  j = 1:N[k]
  x1 = rep(1, N[k])
  x2 = sin(2 * pi * j / 12)#rnorm(N[k])#
  x3 = cos(2 * pi * j / 12)#rnorm(N[k])#cos(2 * pi * j / 12)
  w1 = rep(1, N[k])
  w2 = rep(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), N[k] / 12)#rnorm(N[k])#rep(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), N[k] / 12)
  x = cbind(x1,x2, x3)
  w = cbind(w1, w2)
  xbetas.real = x %*% betas.real
  wlambdas.real = w %*% lambdas.real
  delta.real = exp(wlambdas.real)
  b.real = 1 / delta.real

  for (b in 1:B) {


    a.real = NULL
    y = NULL
    a.real[1] = exp(xbetas.real)[1]
    #%>% median()
    y[1] = (
      extraDistr::rkumar(
        1,
        a = a.real[1],
        b = z[1]*b.real[1]
      )
    )

    for (t in 2:N[k]) {

      a.real[t] = exp(xbetas.real[t])

      repeat {
        y[t] = (
          extraDistr::rkumar(
            1,
            a = a.real[t],
            b = z*b.real[t]
          )
        )

        if (format(y[t], scientific = FALSE) != 1) {
          break
        }
      }
    }


    hist(y)
    #plot.ts(y)
    dados <- data.frame(cbind(y,x2,x3,w2) )
    formula <- y ~  x2 + x3|w2
    mf <- model.frame(Formula::Formula(formula), data = dados)
    y <- model.response(mf)
    cov_a <- model.matrix(Formula::Formula(formula), data = dados, rhs = 1)
    cov_delta <- model.matrix(Formula::Formula(formula), data = dados, rhs = 2)
    ncx <- ncol(cov_a)
    ncv <- ncol(cov_delta)
    # EM--------
    aux.start.betas <- coef(lm(-log(-log(y))~x2+x3))
    aux.start.lambdas <- betareg::betareg(formula = formula,data=dados)$coefficients$precision

    par.cov.start.base.line <- c(aux.start.betas,aux.start.lambdas)

    par.cov.start.marginal = try(
      optim(
        par = par.cov.start.base.line,
        fn = log.like.kumar,
        method = "BFGS",
        control = list(fnscale = -1),
        y = y,
        x = cov_a,
        w = cov_delta,
      )$par,
      silent = T
    )



    if(class(par.cov.start.marginal)=="try-erro"){
      par.cov.start.marginal <- par.cov.start.base.line
    }

    alpha.start <- .5

    estimates.aux <-
      try(optim(
        par = c(par.cov.start.marginal, alpha.start),
        fn = f.cond ,
        control = list(fnscale = -1),
        method = "BFGS",
        #lower = c(rep(-Inf,7),0.25),
        #upper = c(rep(Inf,7),0.95),
        x = x,
        w = w,
        y = y
      ),
      silent = T
      )

    if (class(estimates.aux) == "try-error") {
      next
    }

    alpha.up <- estimates.aux$par[-c(1:(ncx+ncv))]


    if (  alpha.up >= 1) {
      emv <- rep("NA", length(theta.real))
      #estimates.aux$par[5] <- 1/estimates.aux$par[5]
      next
    }


    par.cov.cond.start <- c(  estimates.aux$par[c(1:(ncx+ncv))] )

    erro = 10 ^ (-4)
    iter = 1
    repeat {
      cat("Start_iter=",iter,"\n")
      # betas.start <- theta_start[1:3]
      # lambdas.start <- theta_start[4:5]
      # ar.start <- theta_start[6]
      # ma.start <- theta_start[7]


      # Step E-----

      v <- V(
        theta = par.cov.cond.start,
        x = x,
        w = w,
        y = y
      )


      derivate_numerator <- d_vn(N = N[k],
                                 v = v,
                                 alpha = alpha.up)$dvn

      derivate_denominator <- d_vn(N = N[k]-1,
                                   v = v,
                                   alpha = alpha.up)$dvn


      Esp.z <- -(derivate_numerator / derivate_denominator)
      if (is.finite(Esp.z) == F) {
        break
      }
      # Step M----

      par.cov.cond.up <-
        try(optim(
          par = par.cov.cond.start,
          fn = log.like_conditionl_covariates,
          control = list(fnscale = -1),
          method = "BFGS",
          x = x,
          w = w,
          y = y,
          Etil1 = Esp.z
        ),
        silent = T)

      if (class(par.cov.cond.up) == "try-error") {
        break
      }

      crit <- sum(((par.cov.cond.up$par - par.cov.cond.start) / par.cov.cond.start) ^ 2)
      cat("h=",h,"\n")
      cat("k=",k,"\n")
      cat("b=",b,"\n")
      cat("crit=",crit,"\n")
      cat("Estimate=",c(par.cov.cond.up$par,alpha.up),"\n")

      if (crit < erro) {
        emv <- c(par.cov.cond.up$par,alpha.up)
        EMV <- data.frame(EMV.real=theta.real,EMV.M=estimates.aux$par,EMV=emv,par.names=par.names,k=N[k])
        path_EMV <- "scripts_tests/model_time_v2/Block_1/sample_conditional/EMV.txt"
        write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
        break
      }
      else{
        iter = iter + 1
        par.cov.cond.start <- par.cov.cond.up$par
      }
      if(iter==100){
        break
      }
    }
  }

}
}
toc <- tictoc::toc()
#700.84 sec elapsed
#553.04 sec elapsed
#675.61 sec elapsed
#840.16 sec elapsed
#2912.6 sec elapsed
#Analysis of par -------
#rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
par.names = c(
  paste0("beta",0:2),
  paste0("lambda",0:1),
  "alpha"
)

emv.results = read.table("scripts_tests/model_time_v2/Block_1/sample_conditional/EMV.txt")
colnames(emv.results) = c("id","par.real","est.marg","est.cond","par","N")
emv.results$N <- as.factor(emv.results$N)
emv.results$par <- factor(emv.results$par,levels=par.names)
emv.results <- pivot_longer(emv.results,c( "est.marg","est.cond")) %>% arrange(name)

estimates = emv.results %>%
  group_by(name,N,par) %>%
  summarise(par.real = mean(par.real),
    mean = mean(value,na.rm=T),
            md = median(value,na.rm=T),
            q.025 = quantile(value,c(0.025),na.rm=T),
            q.975 = quantile(value,c(0.975),na.rm=T),
            sdd = sd(value,na.rm = T)) %>% data.frame() %>%
  select(name,N,par,par.real,
          everything())


emv.results %>%
  ggplot(.,aes(N,value,color=name)) +
  geom_boxplot()+
  facet_wrap(~par ,scales = "free_y")+
  geom_hline(aes(yintercept = par.real))

emv.results %>% pivot_wider(c(par))
#Analysis of predicted values
theta <- estimates %>% filter(N==168) %>% select(  md ) %>% unlist()
Beta = theta[1:ncx] # real
lambda = theta[(ncx+1):(ncv+ncx)] # real
alpha = theta[-c(1:(ncv+ncx))][1]
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:N[k]) {
  a[t] = exp(xbeta[t])

}

mediana.mar = ( 1-exp( - delta*(-log(0.5))^(1/.95) ) )^(1/a)
mean((mediana.mar-y)^2)
mean(abs(mediana.bl-y)/y)
# mediana.bl = ( 1-.5^delta )^(1/a)
# plot.ts(y)
# lines(mediana.mar,col=2)
# lines(mediana.bl,col=3)


l = matrix(0,N[k],50)
for(o in 1:50){
  a = NULL
  y.est = NULL
  a[1] = exp(xbeta)[1]

  y.est[1] = (
    pweibull(
      -log(1 - runif(1)),
      shape = 1 / alpha,
      scale = b[1] ^ (alpha),
      lower.tail = TRUE,
      log.p = FALSE
    )
  ) ^ (1 / a.real[1])

  for (t in 2:N[k]) {

    a[t] = exp(xbeta[t])

    repeat {
      u = runif(1)
      y.est[t] = (
        pweibull(
          -log(1 - u),
          shape = 1 / alpha,
          scale = b[t] ^ (alpha),
          lower.tail = TRUE,
          log.p = FALSE
        )
      ) ^ (1 / a[t])

      if (format(y.est[t], scientific = FALSE) != 1) {
        break
      }
    }
  }
  l[,o] <- y.est
  plot.ts(y)
  lines(y.est,col=2)

}
plot.ts(y)
mm = rowMeans(l)
lines(mm,col=2)












al = emv.results %>% filter(V3=="alpha" )
lam = emv.results %>% filter(V3=="lambda0")
require(tidyverse)

cor(matrix(c(al$V2,lam$V2),400,2))




