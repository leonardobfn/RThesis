rm(list=ls())
require(extraDistr)
require(tidyverse)
devtools::load_all()


path_estimation <- "estimations_predicts/estimations.rds"
results <- readRDS(path_estimation)

p <- results$p
q <- results$q
database <-results$data %>% arrange(Group)%>%group_by(City) %>% slice(-(1:(max(q,p))))
quantil <- results$quantil
city <- "Manicoré"
estacao <- results$DADOS %>% filter(idx==city)
cov_delta <- estacao[,c(7,8,9)]
cov_kl <- estacao[,c(3,4,5,6)]
link.kl <-results$link.kl$linkfun
link.eta <-results$link.kl$linkinv
link.delta <- results$link.delta$linkinv
betas <- c(results$kl$Estimate)
lambda <- c(results$delta$Estimate)
ar.ma <- results$ar$Estimate
alfas <- results$alpha[,1]


XB.aux <- as.matrix(cov_kl)%*%betas
XB <- rbind(XB.aux[12],XB.aux)
delta <- link.delta(as.matrix(cov_delta)%*%lambda)
y <- results$data %>% filter(City==city) %>% pull(RH)
r <- eta_klt <- rep(0,108)
for(t in 2:108){
  eta_klt[t] <- XB[t] + ar.ma[1]*(link.kl(y[t-1])-XB[t-1]) + ar.ma[2]*r[t-1]
  r[t] <- link.kl(y[t]) - eta_klt[t]

}

n=1000
mt=.25
B=100
a = d=s=matrix(0,107,B)

for(b in 1:B){
  zm = rstable_pos(n,alfas[3])%>%mean(trim=mt)
  w <- estacao$weights[1]
  H <- zm*w
  betaa = rbeta(1,H,1)
  kl_t <- link.eta(eta_klt[-1])
  bq_t=1/delta
  phi_ts = log(kl_t)/log(1-(1-quantil)^(delta))
  aq_t = 1/phi_ts
  yest = (1-betaa^(delta))^(phi_ts)
  #ycond = (1-(1-runif(1))^(1/(bq_t*H)))^(1/aq_t)
  yk <- rkumar(107,aq_t,bq_t*H)
a[,b] <- yk
#d[,b] <- ycond
s[,b] <- yest
}


yk <- apply(a, 1, quantile,c(0.025,0.5,0.975))
yest <- apply(s, 1, quantile,c(0.025,0.5,0.975))
#yk <- apply(a, 1, median)
#ycond <- apply(s, 1, quantile,c(0.025,0.5,0.975))
plot(estacao$y,type="o",col=1)
#lines(estacao$kl_t,col="2",type="o")
#lines(yest,col="3",type="o")
# lines(yest[3,],col="2",type="o")
# lines(yest[1,],col="2",type="o")
#
# lines(ycond,col="4",type="o")
lines(yk[2,],col="54",type="l")
lines(yk[3,],col="5",type="o")
lines(yk[1,],col="5",type="o")

zm1=NULL
par(mfrow=c(2,2))
for(i in 1:100){
  zm1[i] = rstable_pos(n,alfas[3])%>%mean(trim=.20)

}
hist(zm1)
