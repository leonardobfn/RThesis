rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)

compiler::enableJIT(3)

# Funções-----------

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

##log.f.cond ----------
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

##log.f.cond ----------
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

log.like.kumar = function(theta, x, w, y) {
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
#ff = read.table(file = "scripts_tests/model_time_v2/Block_7/sample_int_alpha_s_1_m2.txt",sep = " ")
#ff = read.table(file = "scripts_tests/model_time_v2/Block_7/sample_int_alpha_1_m2_est2.txt",sep = " ")
ff = read.table(file = "scripts_tests/model_time_v2/Block_7/sample_int_alpha_1_m2_est1.txt",sep = " ")

yy = ff %>% group_by(V2) %>% summarise(m=sample(na.omit(V1),1,prob = NULL,replace=T)) %>%
  pull(m)
plot.ts(yy)
# write.table(yy,"scripts_tests/model_time_v2/Block_7/
#             bom_resultado.txt")
data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
data$precp[data$precp==0] <- 1
data$RH <- yy
#p=10
formula <- RH ~ lles+sent+cost|semester+lles
formula <- RH ~ sent+cost|semester
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
start_aux_lambdas = - betareg::betareg(formula = formula,
                                       link = "loglog",
                                       link.phi = "log",
                                       data=data)$coefficients$precision
fit <- coef(VGAM::vglm(RH ~ lles+sent+cost,family = "kumar", data = data, trace = F))
start_aux_betas <- fit[c(1,3,5,7)]
start_aux_betas = coef(lm(-log(-log(RH)) ~ sent+cost,data=data ))
summary(fit)
start_aux <-
  c(start_aux_betas,-.1,-.1)
emv
# start_aux <-
#   c(start_aux_betas,start_aux_lambdas)
plot.ts(y)
par.cov.Gbase <-
  try(optim(
    par =  c(start_aux),
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
  par.cov.Gbase$par <- c(start_aux_betas,-.1,-.1,-.1)
  break
}


theta.start <-c(start_aux,.5) #c(par.cov.Gbase$par,.5)#
theta.start <-c(par.cov.Gbase$par,.5)#c(start_aux,.5) #
theta.start <-c(emv1$par[-6],0.5) #c(par.cov.Gbase$par,.5)#
#theta.start <-c(emv1$par) #c(par.cov.Gbase$par,.5)#
#theta.start <- c(2.97,0.25,-0.001,-8.36,-.3,.67)
erro = 10 ^ (-4)
iter = 1
# Estimates-------------

beta <- theta.start[c(1:ncx)]
lambda <- theta.start[c((1+ncx):(ncx+ncv))]
alpha <- theta.start[-c((1):(ncx+ncv))]

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
    alpha=alpha,
    hessian=T
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
    alpha=alpha,
    hessian=T
  ),
  silent = T)

emv.lambda.alpha <-
  try(optim(
    par = c(lambda,alpha),
    fn = log.f.cond.lambda.alpha,
    control = list(fnscale = -1),
    method = "BFGS",
    # method = "L-BFGS-B",
    # lower = c(rep(-Inf,ncx+ncv),0.25),
    # upper = c(rep(Inf,ncx+ncv),0.95),
    x = cov_a,
    w = cov_delta,
    y = y,
    beta = emv.beta$par,
    hessian=T
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
    upper = c(1),
    x = cov_a,
    w = cov_delta,
    y = y,
    theta = c(emv.beta$par,emv.lambda$par),
    hessian=T
  ),
  silent = T)

al = seq(0,1,by=0.01)
l=NULL
for(i in 1:length(al)){
  l[i]=log.f.cond.alpha(al[i],c(emv.beta$par,emv.lambda$par),y=y,x = cov_a,
                        w = cov_delta)}
plot(al,l)
# al[which.max(l)]
# if (class(emv.alpha) == "try-error") {
#   emv <- rep("NA", ncv+ncx+1)
#   break
# }

emv1 = list()
emv1$par <- c(emv.beta$par,emv.lambda$par,emv.alpha$par)
#alpha.up <- emv.alpha$par
emv1
#emv <- emv1
acf(y.id)
gg = emv.lambda
ff = sqrt(diag(solve(-gg$hessian)))
gg$par+2.05*ff
gg$par-2.05*ff
emv1
emv
qnorm(0.99)
# hessian <-data.frame(h=emv$hessian,city=citys[p])
# path_hessian <- "scripts_tests/model_time_v2/Block_6/hessian.txt"
# write.table(hessian, file = path_hessian, append = T, quote = T, col.names = F)
# EMV <- data.frame(city=citys[p],emv.cond=emv$par,par.names=par.names)
# path_EMV <- "scripts_tests/model_time_v2/Block_6/EMV.txt"
# write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
results[p,] <- emv$par
#results.mar[p,] <- c(par.cov.Gbase$par,alpha.up)
# y.data[p,] <- y
#}
tic <- tictoc::toc()

rownames(results) <- citys
#Analysis of predicted values -------
data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
data$precp[data$precp==0] <- 1
mf <- model.frame(Formula::Formula(formula), data = data)
#y <- model.response(mf)
cov_a <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)
theta <- emv$par#results[id,] %>% as.numeric()
y.id <- y
Beta = theta[1:ncx] # real
lambda = theta[(ncx+1):c(ncx+ncv)] # real
alpha = theta[-c(1:c(ncx+ncv))]
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL
a[1] = exp(xbeta[1])
for (t in 2:(TT+1)) {
  a[t] = exp(xbeta[t])

}

mediana = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
mean((mediana-y.id)^2)
mean(abs(mediana-y.id)/y)
plot.ts(y.id)
lines(mediana,col=3)
phi.m <- -log(-log(mediana))
phi.y <- -log(-log(y.id))
plot(phi.m,phi.y)
acf(mediana-y.id)


# geração yt dado yt-1  ----------


yt.y1 = y2.dado.y1 = NULL
n = 100
nn=100 # n da unif
tt = 2
at1 = a[tt]
bt1 = b[tt]
at = a[tt+1]
bt = b[tt+1]

u1 = runif(nn)
y1 = (1-exp( - (1/bt1)*(-log(u1))^(1/.85) ) )^(1/at1)
plot.ts(y1)
abline(h=median(y1))
abline(h=y[tt])
abline(h=mediana[tt])
yt.y1[1] <- sample(y1,1,prob = NULL,replace = T)
qyt <- seq(0.001,1,by=0.001)




for(i in 1:n){
  cat(i,"\n")
  u2 = runif(1)
  #u2 <- runif(1,0,u1)
  # int.yt=sapply(qyt,int,
  #               l = 0,
  #               y.t.1 = yt.y1[1],
  #               at1 = a[1],
  #               at = a[2],
  #               bt1 = b[1],
  #               bt = b[2],
  #               alpha = alpha
  # )
  for(j in 1:length(qyt)){
    #cat(j,"\n")
    #j=1
    int.yt = try(int(
      up = qyt[j],
      l = 0,
      y.t.1 = yt.y1[1],
      at1 = at1,
      at = at,
      bt1 = bt1,
      bt = bt,
      alpha = alpha
    ),
    silent = T
    )
    if(class(int.yt)=="try-error"){
      next
    }

    if(abs(u2-int.yt)<10^(-3)){
      y2.dado.y1[i] <- qyt[j]# qyt[which.min(dife)]
      break
    }


  }
  #int.x.y = sapply(qx,int,l=0,y=y)
  # dife = abs(u2-int.yt)
  # y2.dado.y1[i] <- qyt[which.min(dife)]
}
#hist(y2.dado.y1)
plot.ts(na.omit(y2.dado.y1))
abline(h=median(y2.dado.y1,na.rm = T))
abline(h=mediana[tt+1])
abline(h=y[tt+1])


# geração forma geral -----
#alphas = c(0.6770567,.35,.5,.85)
alphas = c(0.7258732 )
n=144;N=200
tic <- tictoc::tic()
for(al in 1:length(alphas)){
  #al=1

  alpha <- alphas[al]
  cat("alpha=",alpha,"\n")
  at1 = a[1]
  bt1 = b[1]
  u1 = runif(N)
  y1 = (1-exp( - (1/bt1)*(-log(u1))^(1/alpha) ) )^(1/at1)
  path.sample.int <- paste0("scripts_tests/model_time_v2/Block_7/sample_int_alpha_",al,"_m2_est3",".txt")
  y.aux <- data.frame(y.sim = y1,t=1,alpha=alpha)
  write.table(y.aux,path.sample.int,sep=" ",append = T,quote = T,row.names = F,col.names = F)
  yt.y1 <- array(c(0),dim=c(n,N))
  yt.y1[1,] <-  y1
  qyt <- seq(0.001,1,by=0.001)
  #seq.i = seq(3,5,2)
  for(i in 2:n){
    #  i=2
    cat(i,"\n")
    at1 = a[i-1]
    bt1 = b[i-1]
    at = a[i]
    bt = b[i]
    #u2 <- runif(1,0,1)
    for(k in 1:N){
      u2 <- runif(1,0,1)
      # int.yt=sapply(qyt,int,
      #               l = 0,
      #               y.t.1 = mean(yt.y1[i-1,]),
      #               at1 = a[i-1],
      #               at = a[i],
      #               bt1 = 1/b[i-1],
      #               bt = 1/b[i],
      #               alpha = alpha
      # )
      for(j in 1:length(qyt)){
        #cat(j,"\n")
        #j=1
        int.yt = try(int(
          up = qyt[j],
          l = 0,
          #y.t.1 = sample(na.omit(yt.y1[i-1,]),1,replace = T,prob = NULL),
          y.t.1 = median(yt.y1[i-1,],na.rm = T),
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha
        ),
        silent = T
        )
        if(class(int.yt)=="try-error"){
          yt.y1[i,k] <- NA
          next
        }

        if(abs(u2-int.yt)<10^(-3)){
          yt.y1[i,k] <- qyt[j]# qyt[which.min(dife)]
          #write.table(EMV, file = path_EMV, append = T, quote = T, col.names = F)
          break
        }


      }

    }

    path.sample.int <- paste0("scripts_tests/model_time_v2/Block_7/sample_int_alpha_",al,"_m2_est3",".txt")
    y.aux <- data.frame(y.sim = yt.y1[i,],t=i,alpha=alpha)
    write.table(y.aux,path.sample.int,sep=" ",append = T,quote = T,row.names = F,col.names = F)
  }
}
toc <- tictoc::toc()

dd = 4
plot.ts(na.omit(yt.y1[dd,]))
abline(h=median(yt.y1[dd,],na.rm = T))
abline(h=mediana[dd])
abline(h=y[dd])
ff = read.table(file = "scripts_tests/model_time_v2/Block_7/sample_int.txt",sep = " ")
t.i = ff %>% filter(V2==144) %>% pull(V1) %>% na.omit()
plot.ts(t.i)
abline(h=median(t.i,na.rm = T))
abline(h=mediana[144])
abline(h=y[144])
ff[nrow(ff),]
41100/300
plot.ts( c(ff[1,]) )
mean(c(ff[1,]),na.rm=T)

ff = read.table(file = "scripts_tests/model_time_v2/Block_7/sample_int.txt",sep = " ")
yy = ff %>% group_by(V2) %>% summarise(m=sample(na.omit(V1),1,prob = NULL,replace=T)) %>%
  pull(m)
plot.ts(yy)
#y.real = y
lines(yy,col=2)
acf(abs(yy))
acf(y.real)


