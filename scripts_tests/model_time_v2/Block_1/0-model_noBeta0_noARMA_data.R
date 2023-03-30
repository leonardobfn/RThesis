rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

##log_like_gbase----

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

# Loading database------
devtools::load_all() # loading my functions

data("data_1")
head(data_1)
## Parse database------
month_names <- factor(month.abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
database <- data_1 %>%
  filter(Station == "1") %>%
  mutate(
    semester= rep(c(1,2),each=6) %>% rep(21) %>% rep(14) %>% as.factor(),
    month_names = rep(month_names, 21) %>% rep(14),
    date = paste0(Year, "/", month_names),
    t = seq(1, 252) %>% rep(14),
    Group = rep(paste0("Group ", c(1, 2, 3, 3, 2, 1, 1, 3, 2, 3, 3, 3, 1, 1)), each = 252),
    cost = cos(2 * pi * as.numeric(Month) / 12),
    sent = sin(2 * pi * as.numeric(Month) / 12),
    lles = (  log(10^((7.5*TBS)/(237.3+TBS)))  )
  )

head(database, 13)
citys <- database$City %>% unique()
results <- matrix(0,length(citys),8)
TT=144
y.data <- matrix(0,length(citys),TT)
for(p in 1:length(citys)){
data <- database %>% filter(City==citys[p]) %>% slice(1:TT)
data$precp[data$precp==0] <- 1
#p=10
formula <- RH ~  lles + sent+cost|semester+log(precp)-1
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)
# write.table(cov_a,"scripts_tests/Block_1/cov_a_manaus_m4s1.txt")
# write.table(cov_delta,"scripts_tests/Block_1/cov_delta_manaus_m4s1.txt")
ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)

tic <- tictoc::tic()
# EM--------
start_aux_lambdas = betareg::betareg(formula = formula,
                                     link = "loglog",
                                     link.phi = "log",
                                     data=data)
start_aux_betas = coef(lm(-log(-log(RH)) ~ lles +sent+cost,data=data ))
start_aux <-
  c(start_aux_betas,-.1,-.1,-.1)

theta.start <-
  c(start_aux, .5)

erro = 10 ^ (-4)
iter = 1
    repeat {
      cat("Start_iter=",iter,"\n")
      # betas.start <- theta_start[1:3]
      # lambdas.start <- theta_start[4:5]
      # ar.start <- theta_start[6]
      # ma.start <- theta_start[7]

      par.covariates.start <-
        theta.start[1:(ncx+ncv)] # covariates parameter
      alpha.start <- theta.start[-c(1:(ncx+ncv))] # alpha parameter


      # Step M----

      estimates.aux <-
        try(optim(
          par = c(par.covariates.start, alpha.start),
          fn = log.like.marginal,
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

      if (class(estimates.aux) == "try-error") {
        emv <- rep("NA", length(theta.start))
        break
      }

      # if (estimates.aux$par[5] >= 1) {
      #   #emv <- rep("NA", length(theta.start))
      #   estimates.aux$par[5] <- 1/ estimates.aux$par[5]
      #   #break
      # }


      alpha.up <- estimates.aux$par[-c(1:(ncx+ncv))]
      par.covariates.start <- estimates.aux$par[c(1:(ncx+ncv))]
      # Step E-----

      v <- V(
        theta = par.covariates.start,
        x = cov_a,
        w = cov_delta,
        y = y
      )


      derivate_numerator <- d_vn(N = TT+1,
                                 v = v,
                                 alpha = alpha.up)$dvn

      derivate_denominator <- d_vn(N = TT,
                                   v = v,
                                   alpha = alpha.up)$dvn


      Esp.z <- -(derivate_numerator / derivate_denominator)
      if (is.finite(Esp.z) == F) {
        emv <- rep("NA", length(theta.start))
        break
      }
      # Step M----

      par.covariates.up <-
        try(optim(
          par = par.covariates.start,
          fn = log.like_conditionl_covariates,
          control = list(fnscale = -1),
          method = "BFGS",
          x = cov_a,
          w = cov_delta,
          y = y,
          Etil1 = Esp.z
        ),
        silent = T)

      if (class(par.covariates.up) == "try-error") {
        emv <- rep("NA", length(theta.start))
        break
      }

      theta.up <-
        c(par.covariates.up$par,
          alpha.up)

      crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)
      # cat("k=",k,"\n")
      # cat("b=",b,"\n")
      cat("crit=",crit,"\n")
      cat("Estimate=",theta.up,"\n")

      if (crit < erro) {
        emv <- theta.up
        #EMV <- data.frame(EMV=emv,par.names=par.names,k=N[k])
        #path_EMV <- "scripts_tests/Block_1/EMV.txt"
        #write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
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
results[p,] <- emv
y.data[p,] <- y
}
rownames(results) <- citys
#Analysis of par -------
#rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
#Analysis of predicted values
id = 14
theta <- results[id,] %>% as.numeric()
y.id <- y.data[id,]
Beta = theta[1:ncx] # real
lambda = theta[(ncx+1):c(ncx+ncv)] # real
alpha = theta[-c(1:c(ncx+ncv))]
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:TT) {
  a[t] = exp(xbeta[t])

}

mediana = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
plot.ts(y.id)
lines(mediana,col=2)
mean((mediana-y.id)^2)
mean(abs(mediana-y.id)/y)
phi.m <- -log(-log(mediana))
phi.y <- -log(-log(y.id))
plot(phi.m,phi.y)
acf(mediana-y.id)

acf(y.id)










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

teta = optim(start_aux,log.like.kumar,x=cov_a,w=cov_delta,control = list(fnscale=-1),y=y.id)$par

Beta.b = theta[1:ncx] # real
lambda.b = theta[(ncx+1):c(ncx+ncv)] # real
xbeta.b = cov_a %*% Beta.b
wlambda.b = cov_delta %*% lambda.b
delta.b = exp(wlambda.b)
b.b = 1 / delta.b
a.b = NULL

a.b[1] = exp(xbeta[1])
for (t in 2:TT) {
  a.b[t] = exp(xbeta.b[t])

}

mediana.base = (1-.5^delta.b)^(1/a.b)
plot.ts(y.id)
lines(mediana.base,col=2)
mean((mediana.base-y.id)^2)
mean(abs(mediana.base-y.id)/y)
lines(mediana,col=3)







