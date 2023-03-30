rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

##log_like_gbase----_-

log.like.kumar = function(theta, x, w, y) {
  #theta =start_aux

  Beta = theta[1:3]
  lambda = theta[c(4:6)]
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

## log_like_marginal------

log.like.marginal = function(theta, x, w, y) {
  #theta = start
  Beta = theta[1:3]
  lambda = theta[c(4:6)]
  # real
  alpha = theta[7]
  xbeta = x %*% Beta
  wlambda = w %*% lambda
  delta = exp(wlambda)
  b = 1 / delta
  a = NULL

  a = exp(xbeta)

  l = sum(
    log(a) + log(alpha) +
      alpha * log(b) + (a - 1) * log(y) -
      log(1 - y ^ (a)) + (alpha - 1) * log(-log(1 - y ^ (a))) -
      (-b * log(1 - y ^ (a))) ^ (alpha)
  ,na.rm = T)
  return(l)
}

## log_like_conditional----

log.like_conditionl_covariates = function(theta, y, x, w, Etil1) {
  # x: betas covariates
  # w: lambdas covariates
  Beta = theta[1:3]
  lambda = theta[c(4:6)]
  # real
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
  Beta = theta[1:3]
  lambda = theta[c(4:6)]
  # real
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
number_of_years = c(5,10,12)

N = number_of_years * 12 # length of sample
betas.real = as.matrix(c(0.8551271  ,0.1957253, -0.0827679))
lambdas.real = as.matrix(c(-7.8704209 ,0.4902060,0.4205015))
alpha.real = 0.8970391
theta.real <- c(betas.real,
                lambdas.real,
                alpha.real)
x <- as.matrix(read.table("scripts_tests/Block_1/cov_a_manaus_m2s1.txt"))[,-1]
w <- as.matrix(read.table("scripts_tests/Block_1/cov_delta_manaus_m2s1.txt"))
par.names = c(
  paste0("beta",1:3),
  paste0("lambda",0:2),
  "alpha"
)
## y generation ----------

emv.results = NULL # store estimates
B = 400 # repeat B times
emv.b = matrix(0, B, length(theta.real))
tic <- tictoc::tic()
for (k in 1:length(number_of_years)) {
  k=1

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
    start_aux <-
      c(0.30835622,0.21981277,-0.04973522,-0.10000000,-0.10000000,-0.10000000 )

    par.covariates.start.aux = try(optim(par=start_aux,fn = log.like.kumar,method = "BFGS",
                                         control = list(fnscale=-1),y=y,x=x,w=w)$par,silent = T)
    if(class(  par.covariates.start.aux)=="try-erro"){
      par.covariates.start.aux <- start_aux
    }
    theta.start <-
      c(  par.covariates.start.aux, .5)

    erro = 10 ^ (-4)
    iter = 1
    repeat {
      cat("Start_iter=",iter,"\n")
      # betas.start <- theta_start[1:3]
      # lambdas.start <- theta_start[4:5]
      # ar.start <- theta_start[6]
      # ma.start <- theta_start[7]

      par.covariates.start <-
        theta.start[1:6] # covariates parameter
      alpha.start <- theta.start[7] # alpha parameter



      # Step M----

      estimates.aux <-
        try(optim(
          par = c(par.covariates.start, alpha.start),
          fn = log.like.marginal,
          control = list(fnscale = -1),
          method = "BFGS",
          #lower = c(rep(-Inf,7),0.25),
          #upper = c(rep(Inf,7),0.95),
          x = x,
          w = w,
          y = y
        ),
        silent = T)

      if (class(estimates.aux) == "try-error") {
        emv <- rep("NA", length(theta.real))
        break
      }

      if (estimates.aux$par[7] >= 1) {
        emv <- rep("NA", length(theta.real))
        #estimates.aux$par[5] <- 1/estimates.aux$par[5]
        break
      }


      alpha.up <- estimates.aux$par[7]
      par.covariates.start <- estimates.aux$par[-7]
      # Step E-----

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
        c(par.covariates.up$par,
          alpha.up)

      crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)
      cat("k=",k,"\n")
      cat("b=",b,"\n")
      cat("crit=",crit,"\n")
      cat("Estimate=",theta.up,"\n")

      if (crit < erro) {
        emv <- theta.up
        EMV <- data.frame(EMV=emv,par.names=par.names,k=N[k])
        path_EMV <- "scripts_tests/Block_1/EMV.txt"
        write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
        break
      }
      else{
        iter = iter + 1
        theta.start <- theta.up
      }
      if(iter==100){
        break
      }
    }
  }

}
toc <- tictoc::toc()
#3509.04 sec elapsed
#Analysis of par -------
#rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
par.names = c(
  paste0("beta",1:3),
  paste0("lambda",0:2),
  "alpha"
)
betas.real = as.matrix(c(0.8551271  ,0.1957253, -0.0827679))
lambdas.real = as.matrix(c(-7.8704209 ,0.4902060,0.4205015))
alpha.real = 0.8970391
theta.real <- c(betas.real,
                lambdas.real,
                alpha.real)
emv.results = read.table("scripts_tests/Block_1/m2s1.txt")
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
n = 60
theta <- estimates %>% filter(V4==n) %>% select(  mean ) %>% unlist()
Beta = theta[1:3] # real
lambda = theta[4:6] # real
alpha = theta[7]
xbeta = x %*% Beta
wlambda = w %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:n) {
  a[t] = exp(xbeta[t])

}

mediana = ( 1-exp( - delta[1:n]*(-log(0.5))^(1/alpha) ) )^(1/a[1:n])
mediana = ( 1-.5^delta[1:n] )^(1/a[1:n])
plot.ts(y)
lines(mediana,col=2)
mean((mediana-y)^2)
mean(abs(mediana-y)/y)

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












al = emv.results %>% filter(V3=="lambda0" & V4=="60")
boxplot(al$V2)
abline(h=-3.0767779)
median(al$V2)
mean(al$V2)
nrow(al)
hist(al$V2)
