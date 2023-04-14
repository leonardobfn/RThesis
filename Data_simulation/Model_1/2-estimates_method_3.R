
rm(list = ls())

#packages--------------
# require(tidyr)
# require(dplyr)
# require(extraDistr)
 source("auxiliary_functions.R")
#compiler::enableJIT(3) #


snowfall.estimates_Method_3 = function(steps,model,alpha,erro=10^(-4)) {
  model=1
   steps = 998
    alpha = 0.95
   print(MC)

  alpha.value <- switch (
    as.character(alpha),
    "0.35" = "alpha35",
    "0.5" = "alpha50",
    "0.65" = "alpha65",
    "0.8" = "alpha80",
    "0.95" = "alpha95"
  )
  path.sample <- paste0(
    "Data_simulation/Model_",
    model,
    "/simulations/",
    alpha.value,
    "/data",
    steps,
    ".txt"
  )


  s = read.table(path.sample)
  covi = read.table("Data_simulation/covariates.txt")
  covi$semester <- as.factor(covi$semester)
  data = data.frame(covi, RH = s$V1)
  formula <- RH ~ sent + cost | semester
  mf <- model.frame(Formula::Formula(formula), data = data)
  y <- model.response(mf)
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)

  ncx <- ncol(cov_a)

  ncv <- ncol(cov_delta)

  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")

  start_aux_betas = try(coef(lm(-log(-log(RH)) ~ sent + cost, data = data)))

  if (class(start_aux_betas) != "try-error") {
    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux_betas, rep(-.1, ncv)),
        fn = log.kumar.gbase,
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
      theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
    }else{
      theta.start = par.cov.Gbase$par
    }
  } else{
    theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
  }


  theta.start = c(theta.start, .5)
  value.start = theta.start

  repeat {
    beta <- theta.start[c(1:ncx)]
    lambda <- theta.start[c((1 + ncx):(ncx + ncv))]
    alpha <- theta.start[-c((1):(ncx + ncv))]

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
        alpha = alpha,
        hessian = T
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
        alpha = alpha,
        hessian = T
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
        upper = c(.99),
        x = cov_a,
        w = cov_delta,
        y = y,
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    if(class(emv.alpha)!="try-error" & emv.alpha$value==0){
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    if (class(emv.alpha) == "try-error"  ) {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }else{
      theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

    crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)

    if (crit < erro) {
      estimates.aux = list()
      estimates.aux$par = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  if (length(which(is.na(estimates.aux$par) == T)) == 0) {
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = emv.beta$hessian
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = emv.lambda$hessian
    emv.ALPHA = list()
    emv.ALPHA$hessian = emv.alpha$hessian
  }else{
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = matrix(0, ncx, ncx)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = matrix(0, ncv, ncv)
    emv.ALPHA = list()
    emv.ALPHA$hessian = 0
  }



  path_estimates = paste0(
    "Data_simulation/Model_",
    1,
    "/estimates/Method_3/estimates/",
    "estimates",
    "alpha95",
    ".txt"
  )
estimates = data.frame(steps,emv,par_names,value.start)
  write.table(
    estimates,
    path_estimates,
    col.names = F,
    row.names = F,
    append = T,
  )

  path_hessian_covBeta = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_3/hessian/",
    "hessianCovBeta_",
    alpha.value,
    ".txt"
  )


  hessianCovBeta = cbind(steps, emv.BETA$hessian)
  write.table(
    hessianCovBeta,
    path_hessian_covBeta,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )
  path_hessian_covLambda = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_3/hessian/",
    "hessianCovLambda_",
    alpha.value,
    ".txt"
  )


  hessianCovLambda = cbind(steps, emv.LAMBDA$hessian)
  write.table(
    hessianCovLambda,
    path_hessian_covLambda,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  path_hessian_alpha = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_3/hessian/",
    "hessian_",
    alpha.value,
    ".txt"
  )

  hessianAlpha = cbind(steps, emv.ALPHA$hessian)
  write.table(
    hessianAlpha,
    path_hessian_alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )
}

# path_estimates = paste0(
#   "Data_simulation/Model_",
#   1,
#   "/estimates/Method_3/estimates/",
#   "estimates",
#   "alpha95",
#   ".txt"
# )
alphas = c(0.95, 0.80, 0.65, 0.50, 0.35) # values real alpha
model = 1
n = 144 # length of series
MC = 1000
cpus <-5
ncv  = 2
ncx = 3

erro = 10 ^ (-4)
alphas = c(0.95)
require(snowfall)

sfInit(cpus=6,type = "SOCK",parallel=TRUE)
sfExportAll()
#sfLibrary("snowfall")
sfLibrary(tidyr)
sfLibrary(dplyr)
sfLibrary(extraDistr)
sfLibrary(Formula)
#sfClusterSetupRNG(seed=7221)
tic <- tictoc::tic()
for (alpha. in alphas) {
  cat(alpha., "\r")
 # alpha. = 0.95
  # sfClusterApplyLB( 1:MC,
  #                   model = model,
  #                   alpha = alpha.,
  #                   fun=snowfall.estimates_Method_3
  #                 )
  sfLapply( 1:MC,
                    model = model,
                    alpha = alpha.,
                    fun=snowfall.estimates_Method_3)
}
toc <- tictoc::toc()
sfStop()
