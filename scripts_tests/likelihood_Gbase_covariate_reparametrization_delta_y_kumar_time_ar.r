rm(list=ls())
# Função Q------
log.like.kumar=function(theta){
  #theta = start

  Beta = theta[1:3]
  lambda = theta[c(4:5)]
  ar = theta[6]
  xbeta = x%*%Beta
  wlambda = w%*%lambda
  delta = exp(wlambda)
  b=1/delta
  a=NULL
  a[1]=exp(xbeta[1])
  for(t in 2:N[k]){

    a[t] = exp(
      xbeta[t] + ar*(-log(-log(y[t-1])-xbetas[t-1]))#+ ma*(-log(-log(y[t-1])) + log(-log(a[t-1])))
    )
  }

  l = sum(extraDistr::dkumar(x=y,a=a,b=b,log = T))
  return(l)
}

# valores dos parâmetros

betas = as.matrix(c(-.5,-.06,-.05))
lambdas = as.matrix(c(1,-1))
ar.real=.5
#start=c(1,0.15,-0.1,1,1,.5)
start=c(betas,lambdas,ar.real)
N=c(1,12,24,36,48)*120
#N=c(1)*120
emv = NULL
B=10
emv.b = matrix(0,B,length(start))
for(k in 1:length(N)){
  #  k=2
  j=1:N[k]
  x1=rep(1,N[k])
  x2=sin(2*pi*j/12)
  x3=cos(2*pi*j/12)
  w1=rep(1,N[k])
  w2=rep(c(1,1,1,1,1,1,0,0,0,0,0,0),N[k]/12)
  x=cbind(x1,x2,x3)
  w=cbind(w1,w2)
  xbetas=x%*%betas
  wlambdas = w%*%lambdas
  delta.real = exp(wlambdas)
  b.real = 1/delta.real
  a.real = NULL
  y = NULL
  a.real[1] = exp(xbetas[1] )



  for(b in 1:B){

    y[1] = extraDistr::rkumar(1,a.real[1],b=b.real[1])


    for(t in 2:N[k]){
      a.real[t] = exp(
        xbetas[t] + ar.real*(-log(-log(y[t-1])-xbetas[t-1]))
      )
      repeat{
        y[t]=extraDistr::rkumar(1,a=a.real[t],b=b.real[t])
        if(format(y[t], scientific = FALSE)!=1){
          break
        }
      }
    }
    hist(y)
    plot.ts(y)
    results = try(
      optim(start,log.like.kumar,method = "BFGS",control=list(fnscale=-1)),silent = T)
    if(class(results)=="try-error"){
      results <- list()
      results$par <- rep("NA",length(start))
    }
    results
    emv.b[b,] = results$par
  }
  emvs = data.frame(N=N[k],emv.b = emv.b)
  emv = rbind(emv,emvs)
}
emv$N.f = factor(emv$N,levels = N)
write.table(emv,"likelihood_Gbase_covariate_reparametrization_delta_y_kumar_time_ar.txt",col.names = T)
emv.results = read.table("likelihood_Gbase_covariate_reparametrization_delta_y_kumar_time_ar.txt")

boxplot(emv.results$emv.b.6~emv.results$N)
abline(h=ar.real)
means.est <- colMeans(subset(emv.results,N==5760)[,-c(1,8)])
means.est <- (means.est)
Beta = means.est[1:3]
lambda = means.est[c(4:5)]
ar =
  xbeta = x%*%Beta
wlambda = w%*%lambda
delta = exp(wlambda)
b=1/delta
a=NULL
a[1]=exp(xbeta[1])
for(t in 2:N[k]){

  a[t] = exp(
    xbeta[t] #+ ar*(xbeta[t-1]) #+ ma*(-log(-log(y[t-1])) + log(-log(a[t-1])))
  )
}

mean.kumar = b*beta(1+(1/a),b)
var.kumar = b*beta(1+(2/a),b) - (b*beta(1+(1/a),b))^2
plot.ts(y)
lines(mean.kumar,col=2)
resi.1 = (y-mean.kumar)/sqrt(var.kumar)
resi.2 = abs(y-mean.kumar)/y
acf(resi.2)
hist(resi.1)
hist(resi.2)
al = evir::hill(abs(resi.2))
plot.ts(al$y[al$y<1])
max(al$y)
y.sort = sort(resi.2)
n=N
k=60
i = seq(0,k-1)
num = (y.sort[n-k]/y.sort[n-i])^2
den = y.sort[n-k]/y.sort[n-i]
plot.ts(num/den)
mean(num/den)

n=N
k=N-1
i = seq(0,k-1)
num = 2*sum( (y.sort[n-k]/y.sort[n-i])^2 ) - sum( (y.sort[n-k]/y.sort[n-i])  )
den = sum( (y.sort[n-k]/y.sort[n-i])  ) - sum((y.sort[n-k]/y.sort[n-i])^2 )
num/den
plot.ts(num/den)

# abline(h=alpha.real)
# H=results$hessian
# s.d = sqrt(diag( -solve(H)))
# EMV = results$par
# EMV -1.96*s.d
# EMV +1.96*s.d


a = mean(log(y.sort[n-i])) - log(y.sort[n-k])
1/a
��
