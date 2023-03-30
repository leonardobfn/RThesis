
rm(list=ls())
require(tidyr)
require(ggplot2)
require(dplyr)

G = function(y,a,b){
  return(  (1-y^(a))^(1/b) )
}

FF = function(y,alpha,G){
  return(exp(-(-log(G))^(alpha)))
}


y =seq(0.01,0.99,by=0.001)
a = c(2)
b = c(2)
alpha = c(.8)
par = expand.grid(a,b,alpha)
results = par %>% group_by(Var1,Var2,Var3) %>%
  summarise(GG = round(G(Var1 ,Var2 ,y=y),5),
            FF.=round(FF(alpha=Var3,y=y,G=GG),5)) %>%
  mutate(a=paste0("a=",Var1),
         b=paste0("b=",Var2),
         al=paste0("al=",Var3),
         y=rep(y,length(FF.)/length(y)),
         dife_GG_FF = abs(GG-FF.))

results
h = exp(-exp(-(results$dife_GG_FF)))
hist(h)
mean(1-h)
hist(results$dife_GG_FF)
summary(results$dife_GG_FF)
id_ponto_inter = which.min(results$dife_GG_FF)
ponto_inter = y[id_ponto_inter]

ggplot(results) +
  geom_line(aes(y,GG),colour="red")+
  geom_line(aes(y,FF.))+
  facet_grid(a+b~al,scales = "free")+
  ylab("")+
  geom_vline(xintercept = ponto_inter)

plot.ts(results$dife_GG_FF)
abline(h=0)
(results%>%data.frame())[which.min(results$dife_GG_FF),]







qb = .58
eta_qm = alpha*(-log(-log(1-qb)))
qm = 1-exp(-exp(-eta_qm))
quantile(y,1-qm)
quantile(y,1-qb)

eta_yb = log(a) + log(-log(1-(1-qb)^b))
yb = exp(-exp(-eta_yb))
eta_ym = log(a) + log(-log(1-exp(-b*(log(1-qm))^(1/alpha)) ))
ym = exp(-exp(-eta_ym))
