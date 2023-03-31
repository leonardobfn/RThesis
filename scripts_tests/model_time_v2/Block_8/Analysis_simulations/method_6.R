rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
# Funções-----------

p.cond <- function(y.t.1, y.t, at1, at, bt1, bt, alpha) {
  G.base.t.1 = extraDistr::pkumar(y.t.1, a = at1, b = bt1,lower.tail = F)
  G.base.t = extraDistr::pkumar(y.t, a = at, b = bt,lower.tail = F)
  LAMBDA.t.1 = -log(G.base.t.1)
  LAMBDA.t = -log(G.base.t)
  g.base.t = extraDistr::dkumar(y.t, a = at, b = bt)
  lambda.t = g.base.t / G.base.t
  return(
    (LAMBDA.t.1 ^ (1-alpha) *
            exp(LAMBDA.t.1 ^ alpha - (LAMBDA.t.1 + LAMBDA.t) ^ alpha) *
            (LAMBDA.t.1 + LAMBDA.t) ^ (alpha -1 )
    )
  )
}
##f.cond----------
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
f.cond.alpha <- function(alpha,theta,y,x,w){
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

# Loading database------
#devtools::load_all() # loading my functions

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
results <- results.mar <- matrix(0,length(citys),8)
TT=143
y.data <- matrix(0,length(citys),TT+1)
tic <- tictoc::tic()
#for(p in 1:length(citys)){
p=10
data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
data$precp[data$precp==0] <- 1
BB = 1000
model=1
emvs = matrix(0,BB,6)
for(ii in 1:BB){
  path <-
    paste0("Data_simulation/Model_",model,"/simulations/","alpha50","/data",ii,".txt")
yy =read.table(path)
data$RH <- yy %>% group_by(V2) %>% summarise(y = sample (V1,1,replace=T,prob=NULL)) %>%
  pull(y)
# data$RH <- yy %>% group_by(V2) %>% summarise(y = median (V1)) %>%
#   pull(y)

 # yy =read.table("scripts_tests/model_time_v2/Block_8/simulation_1/alpha2_m1.txt")
 # data$RH <- yy %>% group_by(V2) %>% summarise(y = sample (V1,1,replace=T,prob=NULL)) %>%
 #   pull(y)



 #p=10
#formula <- RH ~ lles|lles-1
formula <- RH ~ sent + cost|semester
#p=10
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)
# write.table(cov_a,"scripts_tests/Block_1/cov_a_manaus_m4s1.txt")
# write.table(cov_delta,"scripts_tests/Block_1/cov_delta_manaus_m4s1.txt")
ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)
par.names <- c(paste0(colnames(cov_a),"_a"),
               paste0(colnames(cov_delta),"_delta"),
               "alpha")
par.names <- factor(par.names,levels=par.names)
# start values------
start_aux_lambdas = betareg::betareg(formula = formula,
                                     link = "loglog",
                                     link.phi = "log",
                                     data=data)$coefficients$precision

start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost,data=data ))
#start_aux_betas = coef(lm(-log(-log(RH)) ~ lles,data=data ))
start_aux <-
  c(start_aux_betas,-.1,-1)

# start_aux <-
#   c(start_aux_betas,start_aux_lambdas)

par.cov.Gbase <-
  try(optim(
    par =  c(start_aux_betas,-.1,-1),
    fn = log.like.kumar,
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
  par.cov.Gbase <- list()
  par.cov.Gbase$par <- c(start_aux_betas,-.1,-1)
  #break
}


theta.start <- c(start_aux,.5) #c(par.cov.Gbase$par,.45)#
theta.start <- c(par.cov.Gbase$par,.5)#
#theta.start <- c(3.06857482, -0.06830916, 1.08626420, -4.89575000, 0.56244984)
erro = 10 ^ (-4)
iter = 1
# Estimates-------------



emv <-
  try(optim(
    par =theta.start,
    fn = f.cond,
    control = list(fnscale = -1),
    method = "BFGS",
    # method = "L-BFGS-B",
    #  lower = c(rep(-Inf,ncx+ncv),0.50),
    # upper = c(rep(Inf,ncx+ncv),0.95),
    x = cov_a,
    w = cov_delta,
    y = y,
    hessian=T
  ),
  silent = T)

emvs[ii,] <- emv$par
}

colMeans(emvs)
conf = .05
t(round(apply(emvs[1:(ii-1),],2,quantile,c(conf/2,0.5,1-(conf/2))),3))
# if (class(emv) == "try-error") {
#   emv <- rep("NA", ncv+ncx+1)
#   break
# }
#
# # hessian <-data.frame(h=emv$hessian,city=citys[p])
# # path_hessian <- "scripts_tests/model_time_v2/Block_6/hessian.txt"
# # write.table(hessian, file = path_hessian, append = T, quote = T, col.names = F)
# # EMV <- data.frame(city=citys[p],emv.cond=emv$par,par.names=par.names)
# # path_EMV <- "scripts_tests/model_time_v2/Block_6/EMV.txt"
# # write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
# results[p,] <- emv$par
# #results.mar[p,] <- c(par.cov.Gbase$par,alpha.up)
# # y.data[p,] <- y
# }
# tic <- tictoc::toc()
#
# rownames(results) <- citys
# #Analysis of par -------
# #rm(list = ls())
# require(tidyr)
# require(dplyr)
# require(extraDistr)
# require(ggplot2)
# #Analysis of predicted values
# #id = 10
# #p=10
# data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
# data$precp[data$precp==0] <- 1
# #p=10
# #formula <- RH ~  lles+|semester
# mf <- model.frame(Formula::Formula(formula), data = data)
# y <- model.response(mf)
# cov_a <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
# cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)
# theta <- emv$par#results[id,] %>% as.numeric()
# #theta <- results.mar[id,] %>% as.numeric()
# y.id <- y
# Beta = theta[1:ncx] # real
# lambda = theta[(ncx+1):c(ncx+ncv)] # real
# alpha = theta[-c(1:c(ncx+ncv))]
# xbeta = cov_a %*% Beta
# wlambda = cov_delta %*% lambda
# delta = exp(wlambda)
# b = 1 / delta
# a = NULL
#
# a[1] = exp(xbeta[1])
# for (t in 2:(TT+1)) {
#   a[t] = exp(xbeta[t])
#
# }
#
# mediana = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
# mean((mediana-y.id)^2)
# mean(abs(mediana-y.id)/y)
# plot.ts(y.id)
# lines(mediana,col=2)
# phi.m <- -log(-log(mediana))
# phi.y <- -log(-log(y.id))
# plot(phi.m,phi.y)
# acf(mediana-y.id)
# yt.y1 = y2.dado.y1 = NULL
# u1 = runif(1000)
# y1 = (1-exp( - (1/b[1])*(-log(u1))^(1/alpha) ) )^(1/a[1])
# plot.ts(y1)
# abline(h=y[1])
# yt.y1[1] <- median(y1)
# qyt <- seq(0.001,1,by=0.001)
# for(i in 1:200){
#   cat(i,"\n")
#   u2 <- runif(1,0,1)
#   # int.yt=sapply(qyt,int,
#   #               l = 0,
#   #               y.t.1 = yt.y1[1],
#   #               at1 = a[1],
#   #               at = a[2],
#   #               bt1 = 1/b[1],
#   #               bt = 1/b[2],
#   #               alpha = alpha
#   # )
#   for(j in 1:length(qyt)){
#     #cat(j,"\n")
#     #j=1
#     int.yt=int(up=qyt[j],
#                l = 0,
#                y.t.1 = yt.y1[1],
#                at1 = a[1],
#                at = a[2],
#                bt1 = 1/b[1],
#                bt = 1/b[2],
#                alpha = alpha
#     )
#     if(abs(u2-int.yt)<10^(-1)){
#       y2.dado.y1[i] <- qyt[j]# qyt[which.min(dife)]
#       break
#
#     }
#
#
#   }
#   #int.x.y = sapply(qx,int,l=0,y=y)
#   # dife = abs(u2-int.yt)
#   # y2.dado.y1[i] <- qyt[which.min(dife)]
# }
# #hist(y2.dado.y1)
# plot.ts(y2.dado.y1)
# abline(h=y[2])
#







# G.est <- ((pkumar(y ,a,b,lower.tail = F)))
# F.mar.bar = exp(-(  sort( -log(G.est) ) )^alpha)
# log.lam1 <- (-log(F.mar.bar))# weibul
#
# Fy.emp = (ecdf(y))
# Fy.marg = 1-F.mar.bar
#
# plot(Fy.emp)
# lines(sort(y),Fy.marg,col=2)
#
# plot(sort(y),log.lam1)
# plot(log(sort(y)),log(log.lam1))
#
# st <- 1-knots(Fy.emp)
# eta.g <- -log(-log(G.est[1:64]))
# eta.s <- -log(-log(st))[1:64]
# m = matrix(eta.g)
# m1 = matrix(eta.s)
# f = function(al){
#   return(sum((m1-al*m)^2))
# }
# optim(par = .8,fn = f)
# m = matrix(rnorm(50))
# m1 = matrix(4*m)
#
# lm.fit(m,m1)
# log.lam1.emp <- log(lam.1.emp)# weibul
#
# lm(log.lam1.emp~log(-log(G.est[1:64]))-1)
#
# lam.1.emp <- -log(1-knots(Fy.emp))
# log.lam1.emp <- log(lam.1.emp)# weibul
#
# lm(lam.1.emp~G.est[1:length(lam.1.emp)]-1)
# 1/1.831
#
# alpha
#
#
#
#
# log.like.kumar = function(theta, x, w, y) {
#   #theta = par.covariates.start;x=cov_a;w=cov_delta;alpha=.5
#
#   Beta = theta[1:ncx]
#   lambda = theta[-c(1:ncx)]
#   xbeta = x %*% Beta
#   wlambda = w %*% lambda
#   delta = exp(wlambda)
#   b = 1 / delta
#   a = exp(xbeta)
#
#
#
#   l = sum(extraDistr::dkumar(
#     x = y,
#     a = a,
#     b = b,
#     log = T
#   ),na.rm = T)
#   return(l)
# }
#
# teta = optim(start_aux,log.like.kumar,x=cov_a,w=cov_delta,control = list(fnscale=-1),y=y.id)$par
#
# Beta.b = teta[1:ncx] # real
# lambda.b = teta[(ncx+1):c(ncx+ncv)] # real
# xbeta.b = cov_a %*% Beta.b
# wlambda.b = cov_delta %*% lambda.b
# delta.b = exp(wlambda.b)
# b.b = 1 / delta.b
# a.b = NULL
#
# a.b[1] = exp(xbeta.b[1])
# for (t in 2:(TT+1)) {
#   a.b[t] = exp(xbeta.b[t])
#
# }
#
# mediana.base = (1-.5^delta.b)^(1/a.b)
# plot.ts(y.id)
# lines(mediana.base,col=2)
# mean((mediana.base-y.id)^2)
# mean(abs(mediana.base-y.id)/y)
# lines(mediana,col=3)
# acf(mediana.base-y.id)
#
#
diag(solve(-emv$hessian))

sdd = sqrt(diag(solve(-emv$hessian)))
emv$par-2.06*sdd
emv$par+2.06*sdd
emv


