
# thetastart <- estimation_mc %>% filter(par=="alpha_3") %>% pull(estimate)
# ahat=0
# thetahat = 0.307989428
# alpha = 0.025
ciBCA <- function(ahat,thetahat,thetastart,alpha){
  nboot <- length(thetastart)
  z0 <- qnorm(sum(thetastart<thetahat)/nboot)
  alpha1 = pnorm(z0+(z0+qnorm(alpha))/(1-ahat*(z0+qnorm(alpha))))
  alpha2 = pnorm(z0+(z0+qnorm(1-alpha))/(1-ahat*(z0+qnorm(1-alpha))))
  confpoint = quantile(thetastart,probs=c(alpha1,alpha2))
  return(confpoint)
}
