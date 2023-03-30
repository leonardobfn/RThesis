rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

##log_like_gbase----

log.like.kumar = function(theta, x, w, y) {
  #theta = start

  Beta = theta[1:2]
  lambda = theta[c(3:4)]
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
  ))
  return(l)
}
##log_like_margina_alpha-------
log.like.marginal.alpha=function(alpha,theta,x,w,y){
  #theta = start
  Beta = theta[1:2] # real
  lambda = theta[3:4] # real
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  delta = exp(wlambda)
  b = 1/delta

  a = exp(xbeta)

  l = sum( log(a) + log(alpha) +
             alpha*log(b) + (a-1)*log(y) -
             log(1-y^(a)) + (alpha-1)*log(-log(1-y^(a))) -
             (-b*log(1-y^(a)))^(alpha))
  return(l)
}

## log_like_marginal------

log.like.marginal = function(theta, x, w, y) {
  #theta = start
  Beta = theta[1:2] # real
  lambda = theta[3:4] # real
  alpha = theta[5]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta

  a = exp(xbeta)

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
  Beta = theta[1:2] # real
  lambda = theta[3:4] # real
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
  Beta = theta[1:2] # real
  lambda = theta[3:4] # real
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
  #N=N[k];v=v;alpha=.75

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
number_of_years = c(5,10,12)

N = number_of_years * 12 # length of sample
betas.real = as.matrix(c(0.2191953 ,-0.1158410))
lambdas.real = as.matrix(c(-3.0767779, -0.5170042))
alpha.real = .85
theta.real <- c(betas.real,
                lambdas.real,
                alpha.real)
par.names = c(
  paste0("beta",1:2),
  paste0("lambda",0:1),
  "alpha"
)
## y generation ----------

emv.results = NULL # store estimates
B = 400 # repeat B times
emv.b = matrix(0, B, length(theta.real))
tic <- tictoc::tic()
for (k in 1:length(number_of_years)) {
  k=1
  j = 1:N[k]
  x1 = rep(1, N[k])
  x2 = sin(2 * pi * j / 12)#rnorm(N[k])#
  x3 = cos(2 * pi * j / 12)#rnorm(N[k])#cos(2 * pi * j / 12)
  w1 = rep(1, N[k])
  w2 = rep(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), N[k] / 12)#rnorm(N[k])#rep(c(1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0), N[k] / 12)
  x = cbind(x2, x3)
  w = cbind(w1, w2)
  xbetas.real = x %*% betas.real
  wlambdas.real = w %*% lambdas.real
  delta.real = exp(wlambdas.real)
  b.real = 1 / delta.real

  for (b in 1:B) {
    a.real = NULL
    y = NULL
    a.real[1] = exp(xbetas.real)[1]

    y[1] = (
      pweibull(
        -log(1 - runif(1)),
        shape = 1 / alpha.real,
        scale = b.real[1] ^ (alpha.real),
        lower.tail = TRUE,
        log.p = FALSE
      )
    ) ^ (1 / a.real[1])

    for (t in 2:N[k]) {

      a.real[t] = exp(xbetas.real[t])

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
    #plot.ts(y)

    # EM--------


    theta.start.base.line <-  c(c(1,1), c(1,1))


    par.covariates.base.line = try(optim(par=theta.start.base.line,fn = log.like.kumar,method = "BFGS",
                                         control = list(fnscale=-1),y=y,x=x,w=w)$par,silent = T)

    if (class(par.covariates.base.line) == "try-erro") {
      emv <- rep("NA", length(theta.real))
      next
    }

    alpha.up <-
      try(optim(
        par = .5,
        fn = log.like.marginal.alpha,
        control = list(fnscale = -1),
        method = "BFGS",
        #lower = c(rep(-Inf,7),0.25),
        #upper = c(rep(Inf,7),0.95),
        x = x,
        w = w,
        y = y,
        theta = par.covariates.base.line
      )$par,
      silent = T)

    if (class(alpha.up) == "try-error") {
      emv <- rep("NA", length(theta.real))
      next
    }

    if (alpha.up >= 1) {
      emv <- rep("NA", length(theta.real))
      next
    }


    erro = 10 ^ (-4)
    iter = 1
    theta.start <- par.covariates.base.line
    repeat {
      cat("Start_iter=",iter,"\n")

      par.covariates.start <- theta.start[1:4] # covariates parameter

      # Step E----


      v <- V(
        theta = par.covariates.start,
        x = x,
        w = w,
        y = y
      )


      derivate_numerator <- d_vn(N = N[k]+1,
                                 v = v,
                                 alpha = alpha.up)$dvn

      derivate_denominator <- d_vn(N = N[k],
                                   v = v,
                                   alpha = alpha.up)$dvn


      Esp.z <- -(derivate_numerator / derivate_denominator)
      if (is.finite(Esp.z) == F) {
        emv <- rep("NA", length(theta.real))
        break
      }
      # Step M----

      par.covariates.up <-
        try(optim(
          par = par.covariates.start,
          fn = log.like_conditionl_covariates,
          control = list(fnscale = -1),
          method = "BFGS",
          x = x,
          w = w,
          y = y,
          Etil1 = Esp.z
        ),
        silent = T)

      if (class(par.covariates.up) == "try-error") {
        emv <- rep("NA", length(theta.real))
        break
      }

      theta.up <-
        c(par.covariates.up$par)

      crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)
      cat("k=",k,"\n")
      cat("b=",b,"\n")
      cat("crit=",crit,"\n")
      cat("Estimate=",c(theta.up,alpha.up),"\n")

      if (crit < erro) {
        emv <- c(theta.up,alpha.up)
        EMV <- data.frame(EMV=emv,par.names=par.names,k=N[k])
        path_EMV <- "scripts_tests/model_time_v2/Block_2/EMV.txt"
        write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
        break
      }
      else{
        iter = iter + 1
        theta.start <- theta.up
      }
      if(iter==200){
        break
      }
    }
  }
}
toc <- tictoc::toc()

#Analysis of par -------
#rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
par.names = c(
  paste0("beta",1:2),
  paste0("lambda",0:1),
  "alpha"
)
betas.real = as.matrix(c(0.2191953 ,-0.1158410))
lambdas.real = as.matrix(c(-3.0767779, -0.5170042))
alpha.real = .85
theta.real <- c(betas.real,
                lambdas.real,
                alpha.real)
emv.results = read.table("scripts_tests/model_time_v2/Block_2/EMV.txt")
real.values=rep(theta.real,3) %>% as.numeric()
estimates = emv.results %>%
  mutate(V3=factor(V3,levels =par.names)) %>%
  group_by(V4,V3) %>%
  summarise(mean = mean(V2,na.rm=T),
            md = median(V2,na.rm=T),
            q.025 = quantile(V2,c(0.025),na.rm=T),
            q.975 = quantile(V2,c(0.975),na.rm=T),
            sdd = sd(V2,na.rm = T)) %>% data.frame() %>%
  mutate(theta.real = real.values) %>%
  select( V4,  V3,theta.real,
          everything())


real.values = data.frame(V3=par.names,
                         value.real=theta.real) %>%
  mutate(V3=factor(V3,levels = par.names))
emv.results %>% mutate(V3=factor(V3,levels = par.names)) %>%
  ggplot(.,aes(as.factor(V4),V2)) +
  geom_boxplot()+
  facet_wrap(~V3,scales = "free_y")+
  geom_hline(data = real.values,aes(yintercept = value.real))


#Analysis of predicted values
theta <- estimates %>% filter(V4==120) %>% select(  theta.real ) %>% unlist()
Beta = theta[1:2] # real
lambda = theta[3:4] # real
alpha = theta[5]
xbeta = x %*% Beta
wlambda = w %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:N[k]) {
  a[t] = exp(xbeta[t])

}
#mediana.real = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
mediana.est = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
plot.ts(y)
#lines(mediana.real,col=2)
lines(mediana.est,col=4)
mean((mediana.est-y)^2)
mean(abs(mediana.est-y)/y)
phi
sub = emv.results %>% filter(V3=="alpha",V4=="120")
nrow(sub)
#y.s = sort(y)
#y.s = seq(0,1,by=.01)
F.mar = exp( - ( -b[10]*log(1-y.s^(a[10])) )^(alpha)  )
plot(y.s,1-F.mar,type="o")

log.like.kumar = function(theta, x, w, y) {
  #theta = start

  Beta = theta[1:2]
  lambda = theta[c(3:4)]
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
  ))
  return(l)
}

par.base = optim(theta.real[-5],log.like.kumar,method = "BFGS",x=x,y=y,w=w,control = list(fnscale=-1))
Beta = par.base$par[1:2] # real
lambda = par.base$par[3:4] # real
alpha = par.base$par[5]
xbeta = x %*% Beta
wlambda = w %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:N[k]) {
  a[t] = exp(xbeta[t])

}
y.s = seq(0,1,by=0.01)
q.base = (1-(1-.5)^delta)^(1/a)
F.base = pkumar(y.s,a[10],b[10])
lines(y.s,F.base,col=2)
mean((q.base-y)^2)

