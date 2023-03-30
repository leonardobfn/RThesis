

  # rgum <- function(N,mi,sigma,A){
  #   u <- runif(N)
  #   x  <- mi-sigma*(-log(-log(1-(1-u)^(1/A))))
  #
  #
  #   return(x)
  # }

  rgum <- function(N,mi,sigma,A){
    u <- runif(N)
    x  <- mi+sigma*(-log(-log(1-(1-u)^(1/A))))
    return((x))
  }

#
#   N = 1000
#   mi <- 0
#   sigma <- 1
#   A <- 2
#   eps <-(rgum(N,mi=mi,sigma=1,A))
#   mean(eps)
#   hist(eps)
#
#   N = 1000
#   qq=0.5;delta=0.002944610
#   H = 0.9872747
#   mi <- 0 + (log(-log(1-qq^(delta)))) + .70 #log(-log(1-qq^(delta)))#0
#   sigma <- -log(1-qq^(delta))
#   A <- H/delta
#   eps <-(rgum(N,mi=mi,sigma=1,A)) #+ (log(-log(1-qq^(delta))))
#   mean(eps)
#   median(eps)
#   hist(eps)
# #
#
#   N = 1000
#   qq=0.5;delta=0.002944610
#   H = 0.9872747
#   mi <- log(-log(1-qq^(delta)))#0
#   sigma <- 1
#   A <- H/delta
#   eps <-(rgum(N,mi=mi,sigma=1,A))
#   mean(eps)
#   hist(eps)
#
#
#
#
#
#   #
#   l = seq(0,10,1)
#   aux = ((-1)^l*(mi+sigma*exp(-1)+sigma*log(l+1)))/(factorial(l+1)*gamma(A-l))
#   ex = gamma(A+1)*sum((aux))
#   ex
#   abline(v=ex)
