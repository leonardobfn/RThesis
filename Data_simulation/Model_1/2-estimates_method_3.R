rm(list = ls())

#packages--------------
  require(tidyr)
require(dplyr)
require(extraDistr)
source("auxiliary_functions.R")
compiler::enableJIT(3) #


snowfall.estimates_Method_3 = function(steps, model, alpha,erro = 10 ^ (-4)){
  # model=1
  #  steps = 1
  #   alpha = 0.35
  #  print(MC)

  alpha.value <- switch (
    as.character(alpha),
    "0.35" = "alpha35",
    "0.5" = "alpha50",
    "0.65" = "alpha65",
    "0.8" = "alpha80",
    "0.95" = "alpha95"
  )
  path.sample <- paste0("Data_simulation/Model_",
                        model,
                        "/simulations/",
                        alpha.value,
                        "/data",
                        steps,
                        ".txt")


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
  start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost, data = data))
  start_aux <-
    c(start_aux_betas, -.1, -.1)

  par.cov.Gbase <-
    try(optim(
      par =  c(start_aux),
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
    start_aux_lambdas = betareg::betareg(
      formula = formula,
      link = "loglog",
      link.phi = "log",
      data = data
    )$coefficients$precision

    start_aux <-
      c(start_aux_betas, start_aux_lambdas)

    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux),
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
      par.cov.Gbase = list()
      par.cov.Gbase$par <- start_aux
    }

  }

  # estimates.aux <-
  #   try(optim(
  #     par = c(par.cov.Gbase$par,.5),
  #     fn = log.f.cond,
  #     control = list(fnscale = -1),
  #     method = "BFGS",
  #     # method = "L-BFGS-B",
  #     #  lower = c(rep(-Inf,ncx+ncv),0.01),
  #     #  upper = c(rep(Inf,ncx+ncv),1),
  #     x = cov_a,
  #     w = cov_delta,
  #     y = y
  #   ),
  #   silent = T)
  #theta.start = c(start_aux,.5)
  theta.start = c(par.cov.Gbase$par, .5)


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
    if(class(emv.beta)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }
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
    if(class(emv.lambda)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }
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
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    if(class(emv.alpha)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }

    theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    crit <- sum(((theta.up - theta.start)/ theta.start)^2)
    if (crit< erro) {
      estimates.aux = list()
      estimates.aux$par <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  if(length(which(is.na( estimates.aux$par )==T))==0){
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = emv.beta$hessian
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = emv.lambda$hessian
    emv.ALPHA = list()
    emv.ALPHA$hessian = emv.alpha$hessian
  }

  if(length(which(is.na(estimates.aux$par )==T))>0){
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = matrix(0,ncx,ncx)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = matrix(0,ncv,ncv)
    emv.ALPHA = list()
    emv.ALPHA$hessian = 0
  }


  path_estimates = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_3/estimates/",
    "estimates",
    alpha.value
  )

  estimates = data.frame(steps, emv, par_names)
  write.table(
    estimates,
    path_estimates,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
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

  hessianAlpha = cbind(steps,emv.ALPHA$hessian)
  write.table(
    hessianAlpha,
    path_hessian_alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )
}


library(snowfall)
alphas = c(0.95,0.80,0.65,0.50,0.35) # values real alpha
model = 1
n = 144 # length of series
MC = 1000
cpus <- 4
ncv  = 2
ncx=3
#alphas = c(0.35)
sfInit(parallel = TRUE, cpus = cpus)
sfExportAll()
sfLibrary(tidyr)
sfLibrary(dplyr)
sfLibrary(extraDistr)
tic <- tictoc::tic()
for (alpha. in alphas) {
  cat(alpha., "\r")
  # alpha.=0.95
  # sfClusterApplyLB( 1:MC,
  #                   model = model,
  #                   fun=snowfall.estimates_Method_3,
  #                   alpha = alpha.)
  sfLapply(
    1:MC,
    model = model,
    fun=snowfall.estimates_Method_3,
    alpha = alpha.
  )
}
toc <- tictoc::toc()
sfStop()
#1603.73 sec elapsed 95
#2851.64 sec 80
#1726.08 + 346.71 sec elapsed sec elapsed 65
# 976 sec al 50
# 493.37 sec elapsed
# estimate - for -------
erro = 10 ^ (-4)
model = 1
alphas = c(0.95,0.80,0.50, 0.65,0.35)
alphas = c(0.65)
tic <- tictoc::tic()
for(alpha in alphas ){
  alpha = 0.50
  steps = 2

  alpha.value <- switch (
    as.character(alpha),
    "0.35" = "alpha35",
    "0.5" = "alpha50",
    "0.65" = "alpha65",
    "0.8" = "alpha80",
    "0.95" = "alpha95"
  )
  path.sample <- paste0("Data_simulation/Model_",
                        model,
                        "/simulations/",
                        alpha.value,
                        "/data",
                        steps,
                        ".txt")


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
  start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost, data = data))
  start_aux <-
    c(start_aux_betas, -.1, -.1)

  par.cov.Gbase <-
    try(optim(
      par =  c(start_aux),
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
    start_aux_lambdas = betareg::betareg(
      formula = formula,
      link = "loglog",
      link.phi = "log",
      data = data
    )$coefficients$precision

    start_aux <-
      c(start_aux_betas, start_aux_lambdas)

    par.cov.Gbase <-
      try(optim(
        par =  c(start_aux),
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
      par.cov.Gbase = list()
      par.cov.Gbase$par <- start_aux
    }

  }

  # estimates.aux <-
  #   try(optim(
  #     par = c(par.cov.Gbase$par,.5),
  #     fn = log.f.cond,
  #     control = list(fnscale = -1),
  #     method = "BFGS",
  #     # method = "L-BFGS-B",
  #     #  lower = c(rep(-Inf,ncx+ncv),0.01),
  #     #  upper = c(rep(Inf,ncx+ncv),1),
  #     x = cov_a,
  #     w = cov_delta,
  #     y = y
  #   ),
  #   silent = T)
  #theta.start = c(start_aux,.5)
  theta.start = c(par.cov.Gbase$par, .5)


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
    if(class(emv.beta)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }
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
    if(class(emv.lambda)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }
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
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    if(class(emv.alpha)=="try-error"){
      estimates.aux = list()
      estimates.aux$par = rep(NA,ncx+ncv+1)
      break
    }

    theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    crit <- sum(((theta.up - theta.start)/ theta.start)^2)
    if (crit< erro) {
      estimates.aux = list()
      estimates.aux$par <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  if(length(which(is.na(estimates.aux$par )==T))==0){
  emv = list()
  emv <- estimates.aux$par
  emv.BETA = list()
  emv.BETA$hessian = emv.beta$hessian
  emv.LAMBDA = list()
  emv.LAMBDA$hessian = emv.lambda$hessian
  emv.ALPHA = list()
  emv.ALPHA$hessian = emv.alpha$hessian
  }

  if(length(which(is.na(estimates.aux$par )==T))>0){
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = matrix(0,ncx,ncx)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = matrix(0,ncv,ncv)
    emv.ALPHA = list()
    emv.ALPHA$hessian = 0
  }


  path_estimates = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_3/estimates/",
    "estimates",
    alpha.value
  )

  estimates = data.frame(steps, emv, par_names)
  write.table(
    estimates,
    path_estimates,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
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

  hessianAlpha = cbind(steps,emv.ALPHA$hessian)
  write.table(
    hessianAlpha,
    path_hessian_alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )
}
}
toc <- tictoc::toc()


