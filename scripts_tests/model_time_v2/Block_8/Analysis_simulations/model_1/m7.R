rm(list = ls())
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
source("scripts_tests/model_time_v2/Block_8/auxiliary_functions.R")
compiler::enableJIT(3)

# Loading database------
#devtools::load_all() # loading my functions

#data("data_1")
load("/cloud/project/data/data_1.rda")
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
p=10
data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
data$precp[data$precp==0] <- 1
# yy =read.table("model_time_v2/Block_8/simulation_1/alpha3_m1.txt")
# data$RH <- yy %>% group_by(V2) %>% summarise(y = sample (V1,1,replace=T,prob=NULL)) %>%
#   pull(y)

yy =read.table("scripts_tests/model_time_v2/Block_8/simulation_4/alpha1ar1_m2.txt")
# data$RH <- yy %>% group_by(V2) %>% summarise(y = sample (V1,1,replace=T,prob=NULL)) %>%
#   pull(y)
data$RH <- yy %>% group_by(V2) %>% summarise(y = median (V1)) %>%
  pull(y)




#p=10
formula <- RH ~ lles|lles
formula <- RH ~ sent + cost|semester
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
fit <- coef(VGAM::vglm(RH ~ lles,family = "kumar", data = data, trace = F))
start_aux_betas <- fit[c(1,3,5,7)]
start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost,data=data ))
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
theta.start <-c(emv1$par) #c(par.cov.Gbase$par,.5)#
#theta.start <- c(2.97,0.25,-0.001,-8.36,-.3,.67)
#theta.start = c(3.06857482, -0.06830916, 1.08626420, -4.89575000, 0.56244984)
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
emv1
 emv = list()
emv$par=c(2.970735469  ,0.250967153 ,-0.001374638, -8.356448522 ,-0.308273254 , 0.677422743)

#alpha.up <- emv.alpha$par
emv1
#emv <- emv1
acf(y.id)
gg = emv.beta
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
nn = 200 # n da unif
tt = 3
at1 = a[tt]
bt1 = b[tt]
at = a[tt+1]
bt = b[tt+1]

u1 = runif(nn)
y1 = (1-exp( - (1/bt1)*(-log(u1))^(1/alpha) ) )^(1/at1)
plot.ts(y1)
abline(h=median(y1))
abline(h=y[tt])
abline(h=mediana[tt])
yt.y1[1] <- sample(y1,1,prob = NULL,replace = T)
yt.y1[1] <-median(y1)
qyt <- seq(0.001,1,by=0.001)


for(i in 1:nn){
  cat(i,"\n")
  u2 = runif(1)
  int.yt = try (uniroot(
    p.cond,
    interval = c(0, 1),
    u = u2,
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
  else{
    y2.dado.y1[i] <-int.yt$root
  }

}

#hist(y2.dado.y1)
plot.ts(na.omit(y2.dado.y1))
abline(h=median(y2.dado.y1,na.rm = T))
abline(h=mediana[tt+1])
abline(h=y[tt+1])


# geração forma geral -----
alphas = c(0.6770567,.35,.5,.85)
#alphas = c(0.7258732 )
n=144;N=200
tic <- tictoc::tic()
for(al in 1:length(alphas)){
  al=1

  alpha <- alphas[al]
  cat("alpha=",alpha,"\n")
  at1 = a[1]
  bt1 = b[1]
  u1 = runif(N)
  y1 = (1-exp( - (1/bt1)*(-log(u1))^(1/alpha) ) )^(1/at1)
  path.sample.int <- paste0("scripts_tests/model_time_v2/Block_7/sample_root_alpha_",al,"_m2_est1",".txt")
  y.aux <- data.frame(y.sim = y1,t=1,alpha=alpha)
  write.table(y.aux,path.sample.int,sep=" ",append = T,quote = T,row.names = F,col.names = F)
  yt.y1 <- array(c(0),dim=c(n,N))
  yt.y1[1,] <-  y1

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
      int.yt = try (uniroot(
        p.cond,
        interval = c(0, 1),
        u = u2,
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
        next
      }
      else{
        yt.y1[i,k] <-int.yt$root
      }
    }

    path.sample.int <- paste0("scripts_tests/model_time_v2/Block_7/sample_root_alpha_",al,"_m2_est1",".txt")
    y.aux <- data.frame(y.sim = yt.y1[i,],t=i,alpha=alpha)
    write.table(y.aux,path.sample.int,sep=" ",append = T,quote = T,row.names = F,col.names = F)
  }
}
toc <- tictoc::toc()

dd = 144
plot.ts(na.omit(yt.y1[dd,]))
abline(h=median(yt.y1[dd,],na.rm = T))
abline(h=mediana[dd])
abline(h=y[dd])
path.sample.int <- paste0("scripts_tests/model_time_v2/Block_7/sample_root_alpha_",al,"_m2_est1",".txt")
ff = read.table(file = path.sample.int,sep = " ")
t.i = ff %>% filter(V2==144) %>% pull(V1) %>% na.omit()
plot.ts(t.i)
abline(h=median(t.i,na.rm = T))
abline(h=mediana[144])
abline(h=y[144])
ff[nrow(ff),]
41100/300
plot.ts( c(ff[1,]) )
mean(c(ff[1,]),na.rm=T)

ff = read.table(file = path.sample.int,sep = " ")
yy = ff %>% group_by(V2) %>% summarise(m=median(V1)) %>%
  pull(m)
plot.ts(yy)
#y.real = y
lines(yy,col=2)
acf(abs(yy))
acf(y.real)


