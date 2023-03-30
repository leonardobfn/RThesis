rm(list=ls())
require(extraDistr)
require(tidyverse)
devtools::load_all()

#--- link function-----

linkfun = function(mu) -log(-log(mu))

linkinv = function(eta) {
  pmax(pmin(
    exp(-exp(-eta)),
    1 - .Machine$double.eps
  ), .Machine$double.eps)
}
#---------data organization ------
TT <- 108
one <- rep(1,2*108)
tt <- rep(1:12,9) %>% rep(2)
city <- rep(c(1,0),each=108)
idx <- rep(c("A","B"),each=108)
sent <- sin(2*pi*tt/12)
cost <- cos(2*pi*tt/12)

X. <- data.frame(one,sent,cost)
W. <- data.frame(one,sent,cost)

betas <- c(1.50,.10,-0.015)
lambda <- c(-5,0.1,-0.04)
ar.ma <- c(0.57,-.10)

XB <- as.matrix(X.)%*%betas
WL <- as.matrix(W.)%*%lambda
delta <- exp(WL)
quantil  <- .5


n=1000
mt=.20
set.seed(1)
zm = rstable_pos(n,.50)%>%mean(trim=mt)
w <- .5
H <- zm*w


y1 <- rep(0,108)
y1[1] <- linkinv(XB[1])
r1 <- eta_k1t <- rep(0,108)

for(t in 2:108){

  eta_k1t[t] <- XB[t] + ar.ma[1]*(linkfun(y1[t-1])-XB[t-1]) + ar.ma[2]*r1[t-1]
  betaa = rbeta(1,H,1)
  bq_t=1/delta[t]
  kl_t <- linkinv(eta_k1t[t])
  bq_t=1/delta[t]
  phi_ts = log(kl_t)/log(1-(1-quantil)^(delta[t]))
  aq_t = 1/phi_ts
  y1[t] <- median(rkumar(200,aq_t,bq_t*H))
  r1[t] <- linkfun(y1[t]) - eta_k1t[t]

}


y2 <- rep(0,108)
y2[1] <- linkinv(XB[109])
r2 <- eta_k2t <- rep(0,108)

for(t in 2:108){

  eta_k2t[t] <- XB[-c(1:108)][t] + ar.ma[1]*(linkfun(y2[t-1])-XB[-c(1:108)][t]) + ar.ma[2]*r2[t-1]
  bq_t=1/delta[t]
  kl_t <- linkinv(eta_k2t[t])
  bq_t=1/delta[t]
  phi_ts = log(kl_t)/log(1-(1-quantil)^(delta[t]))
  aq_t = 1/phi_ts
  y2[t] <- median(rkumar(200,aq_t,bq_t*H))
  r2[t] <- linkfun(y2[t]) - eta_k2t[t]

}

cor(matrix(c(y1,y2),108,2))
plot(y1,y2)
y = c(y1,y2)

pesos <- rep(w,2*108)
grupos <- rep("G1",2*108)
dados <- data.frame(y=y,X.,W.,idx,grupos,weights=pesos)
data <- dados
head(data)

formula <- y ~ sent + cost    | sent + cost
# data <- estimations$data
grupos <- data$grupos
w <- data$weights
quantil <- .5
N <- 500
erro <- 10^(-4)
link.kl <- "loglog"
link.delta <- "log"
p <- 1 # AR(p)
q <- 1 # MA(q)
idx <- data$idx
tic <- tictoc::tic()
x <- gev_kuma(
  formula = formula,
  p = p,
  q = q,
  grupos = grupos,
  w = w,
  N = N,
  quantil = quantil,
  data = data,
  link.kl = "loglog", link.delta = "log", erro = erro,
  id = data$idx
)
#----------- analise-----

betas <- c(x$kl$Estimate)
lambda <- c(x$delta$Estimate)
ar.ma <- c(x$ar$Estimate)

XB <- as.matrix(X.)%*%betas
WL <- as.matrix(W.)%*%lambda
delta <- exp(WL)
quantil  <- .5


n=1000
mt=.20
#set.seed(1)
zm = rstable_pos(n,x$alpha$Estimate)%>%mean(trim=mt)
w <- .5
H <- zm*w


y1 <- rep(0,108)
y1[1] <- linkinv(XB[1])
r1 <- eta_k1t <- rep(0,108)

for(t in 2:108){

  eta_k1t[t] <- XB[t] + ar.ma[1]*(linkfun(y1[t-1])-XB[t-1]) + ar.ma[2]*r1[t-1]
  betaa = rbeta(1,H,1)
  bq_t=1/delta[t]
  kl_t <- linkinv(eta_k1t[t])
  bq_t=1/delta[t]
  phi_ts = log(kl_t)/log(1-(1-quantil)^(delta[t]))
  aq_t = 1/phi_ts
  y1[t] <- median(rkumar(200,aq_t,bq_t*H))
  r1[t] <- linkfun(y1[t]) - eta_k1t[t]

}


y2 <- rep(0,108)
y2[1] <- linkinv(XB[109])
r2 <- eta_k2t <- rep(0,108)

for(t in 2:108){

  eta_k2t[t] <- XB[-c(1:108)][t] + ar.ma[1]*(linkfun(y2[t-1])-XB[-c(1:108)][t]) + ar.ma[2]*r2[t-1]
  bq_t=1/delta[t]
  kl_t <- linkinv(eta_k2t[t])
  bq_t=1/delta[t]
  phi_ts = log(kl_t)/log(1-(1-quantil)^(delta[t]))
  aq_t = 1/phi_ts
  y2[t] <- rkumar(1,aq_t,bq_t*H)
  r2[t] <- linkfun(y2[t]) - eta_k2t[t]

}

cor(matrix(c(y1,y2),108,2))
plot(y1,y2)
y3 = c(y1,y2)
plot.ts(y)
lines(y3,col=2)
abline(v=108)
mean(abs(y3-y)/y)
