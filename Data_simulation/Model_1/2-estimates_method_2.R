rm(list = ls())

packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
source("auxiliary_functions.R")
compiler::enableJIT(3)


snowfall.estimates_method_2 = function(steps, model, alpha,erro = 10 ^ (-4)){
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

    path.sample <-

      paste0("Data_simulation/Model_",
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


      if (abs(emv.alpha$par - alpha) < erro) {
        estimates.aux = list()
        estimates.aux$par <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
        break
      } else{
        theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      }

    }

     alpha.up <- estimates.aux$par[-c(1:(ncx + ncv))]
     par.covariates.start <- estimates.aux$par[c(1:(ncx + ncv))]
      #print(estimates.aux)
    repeat {

      if(length(which(is.na(estimates.aux$par)==T))>0){
        par.covariates.up = list()
        par.covariates.up$value =0
        emv <- rep(NA,ncx+ncv+1)
        break
      }
#       # Step E-----
      v <- V(
        theta = par.covariates.start,
        x = cov_a,
        w = cov_delta,
        y = y
      )


      derivate_numerator <- try(d_vn(N = length(y) + 1,
                                 v = v,
                                 alpha = alpha.up)$dvn)

      if(class(derivate_numerator)=="try-error" | is.finite(derivate_numerator)==F |
         derivate_numerator==0){
        par.covariates.up = list()
        par.covariates.up$value =0
        emv <- rep(NA,ncx+ncv+1)
        break
      }

      derivate_denominator <- try(d_vn(N = length(y),
                                   v = v,
                                   alpha = alpha.up)$dvn)

      if(class(derivate_denominator)=="try-error" | is.finite(derivate_denominator)==F |
         derivate_denominator==0){
        par.covariates.up = list()
        par.covariates.up$value =0
        emv <- rep(NA,ncx+ncv+1)
        break
      }

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

      if(class(par.covariates.up)=="try-error"){
        par.covariates.up = list()
        par.covariates.up$value =0
        emv <- rep(NA,ncx+ncv+1)
        break
      }

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

   if(class(par.covariates.up)!="try-error" & par.covariates.up$value ==0)
   {
      emv <- rep(NA,ncx+ncv+1)
    }


    path_estimates = paste0(
      "Data_simulation/Model_",
      model,
      "/estimates/Method_2/estimates/",
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

    # if(length(which(is.na(emv)==T)==0)){
    # colnames(par.covariates.up$hessian) = par_names[-6]
    # rownames(par.covariates.up$hessian) = par_names[-6]
    # }

    path_hessian_cov = paste0(
      "Data_simulation/Model_",
      model,
      "/estimates/Method_2/hessian/",
      "hessianCov_",
      alpha.value,
      ".txt"
    )
    if(length(which(is.na(emv)==T)>0)){
      par.covariates.up = list()
      par.covariates.up$hessian <- matrix(0,ncx+ncv,ncx+ncv)
    }

    hessianCov = cbind(steps, par.covariates.up$hessian)
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
    if(length(which(is.na(emv)==T)>0)){
      emv.alpha = list()
      emv.alpha$hessian <- matrix(0,1,1)
    }
    hessianAlpha = cbind(steps, emv.alpha$hessian)
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
data.labels <- paste0("data", 1:MC, ".txt")
alphas = c(0.35)
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
  #                   fun=snowfall.estimates_method_2,
  #                   alpha = alpha.)
  sfLapply(
    1:MC,
    model = model,
    fun=snowfall.estimates_method_2,
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
  alpha = 0.65
alpha.value <- switch (
  as.character(alpha),
  "0.35" = "alpha35",
  "0.5" = "alpha50",
  "0.65" = "alpha65",
  "0.8" = "alpha80",
  "0.95" = "alpha95"
)

id = c(1:5)
for (steps in id) {
  steps=84
  path.sample <-

    paste0("Data_simulation/Model_",
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


    if (abs(emv.alpha$par - alpha) < erro) {
      estimates.aux = list()
      estimates.aux$par <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
      break
    } else{
      theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    }

  }

  alpha.up <- estimates.aux$par[-c(1:(ncx + ncv))]
  par.covariates.start <- estimates.aux$par[c(1:(ncx + ncv))]
  #print(estimates.aux)
  repeat {

    if(length(which(is.na(estimates.aux$par)==T))>0){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
      break
    }
    #       # Step E-----
    v <- V(
      theta = par.covariates.start,
      x = cov_a,
      w = cov_delta,
      y = y
    )


    derivate_numerator <- try(d_vn(N = length(y) + 1,
                                   v = v,
                                   alpha = alpha.up)$dvn)

    if(class(derivate_numerator)=="try-error" | is.finite(derivate_numerator)==F |
       derivate_numerator==0){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
      break
    }

    derivate_denominator <- try(d_vn(N = length(y),
                                     v = v,
                                     alpha = alpha.up)$dvn)

    if(class(derivate_denominator)=="try-error" | is.finite(derivate_denominator)==F |
       derivate_denominator==0){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
      break
    }

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

    if(class(par.covariates.up)=="try-error"){
      par.covariates.up = list()
      par.covariates.up$value =0
      emv <- rep(NA,ncx+ncv+1)
      break
    }

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

  if(class(par.covariates.up)!="try-error" & par.covariates.up$value ==0){
    emv <- rep(NA,ncx+ncv+1)
  }


  path_estimates = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_2/estimates/",
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

  # if(length(which(is.na(emv)==T)==0)){
  # colnames(par.covariates.up$hessian) = par_names[-6]
  # rownames(par.covariates.up$hessian) = par_names[-6]
  # }

  path_hessian_cov = paste0(
    "Data_simulation/Model_",
    model,
    "/estimates/Method_2/hessian/",
    "hessianCov_",
    alpha.value,
    ".txt"
  )
  if(length(which(is.na(emv)==T)>1)){
    par.covariates.up = list()
    par.covariates.up$hessian <- matrix(0,ncx+ncv,ncx+ncv)
  }

  hessianCov = cbind(steps, par.covariates.up$hessian)
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
  if(length(which(is.na(emv)==T)>1)){
    emv.alpha = list()
    emv.alpha$hessian <- matrix(0,1,1)
  }
  hessianAlpha = cbind(steps, emv.alpha$hessian)
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


