rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
source("auxiliary_functions.R")
compiler::enableJIT(3)



snowfall.estimates_method_2 = function(steps, model, alpha,erro = 10 ^ (-4)){
  # model=1
  # steps = 1
  #  alpha = 0.35
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
                         "/",
                         "data",
                         steps,
                         ".txt")


  s = read.table(path.sample)

  covi = read.table("Data_simulation/covariates.txt")
  data = data.frame(covi, RH = s$V1)
  formula <- RH ~ sent + cost | semester
  mf <- model.frame(Formula::Formula(formula), data = data)
  y <- model.response(mf)
  #print(y[1])
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)
  #print(cov_delta)
  ncx <- ncol(cov_a)
  ncv <- ncol(cov_delta)
  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")

  start_aux_betas = coef(lm(-log(-log(RH)) ~ sent + cost, data = data))
  start_aux <-
    c(start_aux_betas, -.1, -.1)
#print(start_aux)

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

#print("asdasd")

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
      par.cov.Gbase$par <- c(start_aux)
    }

  }
   #print(par.cov.Gbase$par)
  #
  estimates.aux <-
    try(optim(
      par = c(par.cov.Gbase$par,.5),
      fn = log.f.cond,
      control = list(fnscale = -1),
      method = "BFGS",
      # method = "L-BFGS-B",
      #  lower = c(rep(-Inf,ncx+ncv),0.01),
      #  upper = c(rep(Inf,ncx+ncv),1),
      x = cov_a,
      w = cov_delta,
      y = y
    ),
    silent = T)

  # theta.start = c(par.cov.Gbase$par,.5)
  # estimates.aux = list()
  #
  # repeat{
  #
  #   beta <- theta.start[c(1:ncx)]
  #   lambda <- theta.start[c((1+ncx):(ncx+ncv))]
  #   alpha <- theta.start[-c((1):(ncx+ncv))]
  #
  #   emv.beta <-
  #     try(optim(
  #       par = beta,
  #       fn = log.f.cond.beta,
  #       control = list(fnscale = -1),
  #       method = "BFGS",
  #       # method = "L-BFGS-B",
  #       # lower = c(rep(-Inf,ncx+ncv),0.25),
  #       # upper = c(rep(Inf,ncx+ncv),0.95),
  #       x = cov_a,
  #       w = cov_delta,
  #       y = y,
  #       lambda = lambda,
  #       alpha=alpha,
  #     ),
  #     silent = T)
  #   print(emv.beta$par)
  #
  #   emv.lambda <-
  #     try(optim(
  #       par = lambda,
  #       fn = log.f.cond.lambda,
  #       control = list(fnscale = -1),
  #       method = "BFGS",
  #       # method = "L-BFGS-B",
  #       # lower = c(rep(-Inf,ncx+ncv),0.25),
  #       # upper = c(rep(Inf,ncx+ncv),0.95),
  #       x = cov_a,
  #       w = cov_delta,
  #       y = y,
  #       beta = emv.beta$par,
  #       alpha=alpha,
  #     ),
  #     silent = T)

    # emv.alpha <-
    #   try(optim(
    #     par = alpha,
    #     fn = log.f.cond.alpha,
    #     control = list(fnscale = -1),
    #     method = "BFGS",
    #     # method = "L-BFGS-B",
    #     # lower = c(0.1),
    #     # upper = c(1),
    #     x = cov_a,
    #     w = cov_delta,
    #     y = y,
    #     theta = c(emv.beta$par,emv.lambda$par),
    #     hessian=T
    #   ),
    #   silent = T)

    # if(abs(emv.alpha$par-alpha)<erro){
    #   estimates.aux$par <- c(emv.beta$par,emv.lambda$par,emv.alpha$par)
    #   break
    # }else{
    #   theta.start = c(emv.beta$par,emv.lambda$par,emv.alpha$par)
    # }

  #}
  #
  alpha.up <- estimates.aux$par[-c(1:(ncx+ncv))]
  par.covariates.start <- estimates.aux$par[c(1:(ncx+ncv))]
  print(alpha.up)
  #
  # repeat {
  #   # Step E-----
  #   v <- V(theta = par.covariates.start,
  #          x = cov_a,
  #          w = cov_delta,
  #          y = y)
  #
  #
  #   derivate_numerator <- d_vn(N = length(y) + 1,
  #                              v = v,
  #                              alpha = alpha.up)$dvn
  #
  #   derivate_denominator <- d_vn(N = length(y),
  #                                v = v,
  #                                alpha = alpha.up)$dvn
  #
  #
  #   Esp.z <- -(derivate_numerator / derivate_denominator)
  #
  #   # Step M----
  #
  #   par.covariates.up <-
  #     try(optim(
  #       par = par.covariates.start,
  #       fn = log.like_conditionl_covariates_EM,
  #       control = list(fnscale = -1),
  #       method = "BFGS",
  #       x = cov_a,
  #       w = cov_delta,
  #       y = y,
  #       Etil1 = Esp.z,
  #       hessian = T
  #     ),
  #     silent = T)
  #
  #   crit <-
  #     sum(((par.covariates.up$par - par.covariates.start) / par.covariates.start
  #     ) ^ 2)
  #   if (crit < erro) {
  #     emv <- c(par.covariates.up$par, alpha.up)
  #     break
  #   }
  #   else{
  #     par.covariates.start <- par.covariates.up$par
  #   }
  # }
  #
  #
  # path_estimates = paste0(
  #   "Data_simulation/Model_",
  #   model,
  #   "/estimates/Method_2/estimates/",
  #   "estimates",
  #   alpha.value
  # )
  #
  # estimates = data.frame(steps, emv, par_names)
  # write.table(
  #   estimates,
  #   path_estimates,
  #   col.names = F,
  #   row.names = F,
  #   quote = T,
  #   append = T
  # )
  #
  # colnames(par.covariates.up$hessian) = par_names[-6]
  # rownames(par.covariates.up$hessian) = par_names[-6]
  #
  # path_hessian_cov = paste0(
  #   "Data_simulation/Model_",
  #   model,
  #   "/estimates/Method_2/hessian/",
  #   "_hessianCov_",
  #   alpha.value,
  #   ".txt"
  # )
  #
  # hessianCov = cbind(1, par.covariates.up$hessian)
  # write.table(
  #   hessianCov,
  #   path_hessian_cov,
  #   col.names = F,
  #   row.names = F,
  #   quote = T,
  #   append = T
  # )
  #
  # path_hessian_alpha = paste0(
  #   "Data_simulation/Model_",
  #   model,
  #   "/estimates/Method_2/hessian/",
  #   "_hessian_",
  #   alpha.value,
  #   ".txt"
  # )
  #
  # hessianAlpha = cbind(1, emv.alpha$hessian)
  # write.table(
  #   hessianAlpha,
  #   path_hessian_alpha,
  #   col.names = F,
  #   row.names = F,
  #   quote = T,
  #   append = T
  # )
}

library(snowfall)
alphas = c(0.35,0.50, 0.65, 0.80, 0.95) # values real alpha
model = 1
n = 144 # length of series
MC = 5
cpus <- 4
data.labels <- paste0("data", 1:MC, ".txt")

sfInit(parallel = TRUE, cpus = cpus)
sfExportAll()
#tic <- tictoc::tic()
#for (alpha. in alphas) {
  cat(alpha., "\r")
  alpha.=0.95
  sfClusterApplyLB( 1:MC,
                    model = model,
                    fun=snowfall.estimates_method_2,
                    alpha = alpha.)
  sfLapply(
    1:MC,
    model = model,
    fun=snowfall.estimates_method_2,
    alpha = alpha.
  )
}
toc <- tictoc::toc()
sfStop()




# estimate - for -------
model = 1
alphas = c(0.50, 0.65)
for(alpha in alphas ){
alpha.value <- switch (
  as.character(alpha),
  "0.35" = "alpha35",
  "0.5" = "alpha50",
  "0.65" = "alpha65",
  "0.8" = "alpha80",
  "0.95" = "alpha95"
)

id = c(1:1000)
for (MC in id) {
  # MC=1
  path.sample <-

    paste0("Data_simulation/Model_",
           model,
           "/simulations/",
           alpha.value,
           "/data",
           MC,
           ".txt")


  s = read.table(path.sample)

  covi = read.table("Data_simulation/covariates.txt")
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

  # par_names <- c(paste0(colnames(cov_a), "_a"),
  #                paste0(colnames(cov_delta), "_delta"),
  #                "alpha")

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

  theta.start = c(par.cov.Gbase$par, .5)
  estimates.aux = list()

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
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    if (abs(emv.alpha$par - alpha) < erro) {
      estimates.aux$par <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }



  alpha.up <- estimates.aux$par[-c(1:(ncx + ncv))]
  par.covariates.start <- estimates.aux$par[c(1:(ncx + ncv))]

  repeat {
    # Step E-----
    v <- V(
      theta = par.covariates.start,
      x = cov_a,
      w = cov_delta,
      y = y
    )


    derivate_numerator <- d_vn(N = length(y) + 1,
                               v = v,
                               alpha = alpha.up)$dvn

    derivate_denominator <- d_vn(N = length(y),
                                 v = v,
                                 alpha = alpha.up)$dvn


    Esp.z <- -(derivate_numerator / derivate_denominator)

    # Step M----

    par.covariates.up <-
      try(optim(
        par = par.covariates.start,
        fn = log.like_conditionl_covariates_EM,
        control = list(fnscale = -1),
        method = "BFGS",
        x = cov_a,
        w = cov_delta,
        y = y,
        Etil1 = Esp.z,
        hessian = T
      ),
      silent = T)

    crit <-
      sum(((par.covariates.up$par - par.covariates.start) / par.covariates.start
      ) ^ 2)
    if (crit < erro) {
      emv <- c(par.covariates.up$par, alpha.up)
      break
    }
    else{
      par.covariates.start <- par.covariates.up$par
    }
  }


  path_estimates = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_2/estimates/",
    "estimates",
    alpha.value
  )

  estimates = data.frame(MC, emv, par_names)
  write.table(
    estimates,
    path_estimates,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  colnames(par.covariates.up$hessian) = par_names[-6]
  rownames(par.covariates.up$hessian) = par_names[-6]

  path_hessian_cov = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_2/hessian/",
    "hessianCov_",
    alpha.value,
    ".txt"
  )

  hessianCov = cbind(MC, par.covariates.up$hessian)
  write.table(
    hessianCov,
    path_hessian_cov,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  path_hessian_alpha = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_2/hessian/",
    "hessian_",
    alpha.value,
    ".txt"
  )

  hessianAlpha = cbind(MC, emv.alpha$hessian)
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



