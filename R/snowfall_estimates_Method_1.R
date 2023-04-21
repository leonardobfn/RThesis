snowfall.estimates_Method_1 = function(steps,
                                       model,
                                       alpha.value,
                                       erro = 10 ^ (-4),
                                       formula,
                                       wd_sample,
                                       wd_results_emv,
                                       wd_results_hessian,
                                       wd,
                                       wd_covariates,
                                       n
)
{
  # model=1
  #  steps = 650
  #
  #print(MC)

  # setwd(wd_sample)
  # path.sample = paste0("Model_", model,"/simulations/",alpha.value,"/")
  # #file.sample = paste0( path.sample, data.labels[steps])
  # file.sample = paste0( path.sample,"data",steps,".txt")
  # data.sample = read.table(file.sample)
  # erro = 10^(-4)
  # wd_covariates=wd.covariates
  # wd_results_emv = wd.results_emv
  # wd_results_hessian = wd.results_emv
  # wd_sample = wd.sample
  # wd= wd.
  # model = 1
  # alpha.value = "alpha50"
  # formula = FORMULA
  # steps = 1025
  # n = n.

  setwd(wd)
  setwd(wd_sample)
  file.sample = paste0(alpha.value,"/data",steps,".txt")
  data.sample = read.table(file.sample)
  setwd(wd)
  setwd(wd_covariates)
  covi = read.table("covariates.txt")[1:n,]

  covi$semester <- as.factor(covi$semester)

  data = data.frame(covi, RH = data.sample$y)
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
  par_real = read.table(paste0("real_par_model_", model, ".txt"))


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
        y = y,
        ncx=ncx,
        ncv=ncv
      ),
      silent = T)

    if (class(par.cov.Gbase) != "try-error" & par.cov.Gbase$value==0) {
      theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
    }

    if (class(par.cov.Gbase) == "try-error") {
      theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
    }else{
      theta.start = par.cov.Gbase$par
    }
  } else{
    theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
  }


  theta.start = c(par.cov.Gbase$par,0.5)
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
        ncx = ncx,
        ncv = ncv,
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
        ncx = ncx,
        ncv = ncv,
        theta = c(emv.beta$par, emv.lambda$par),
        hessian = T
      ),
      silent = T)

    # if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
    #   estimates.aux = list()
    #   estimates.aux$par = rep(NA, ncx + ncv + 1)
    #   break
    # }
    #
    # if (class(emv.alpha) == "try-error") {
    #   estimates.aux = list()
    #   estimates.aux$par = rep(NA, ncx + ncv + 1)
    #   break
    # } else{
    #   theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
    # }

    if (class(emv.alpha) == "try-error") {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
      estimates.aux = list()
      estimates.aux$par = rep(NA, ncx + ncv + 1)
      break
    }

    theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)


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
    emv <- round(estimates.aux$par, 4)
    emv.BETA = list()
    emv.BETA$hessian = round(emv.beta$hessian, 4)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = round(emv.lambda$hessian, 4)
    emv.ALPHA = list()
    emv.ALPHA$hessian = round(emv.alpha$hessian, 4)
  } else{
    emv = list()
    emv <- estimates.aux$par
    emv.BETA = list()
    emv.BETA$hessian = matrix(0, ncx, ncx)
    emv.LAMBDA = list()
    emv.LAMBDA$hessian = matrix(0, ncv, ncv)
    emv.ALPHA = list()
    emv.ALPHA$hessian = 0
  }
  setwd(wd)
  setwd(wd_results_emv)
  #steps=1
  file.estimates.emv  = paste0("estimates_",alpha.value,"_n",n, ".txt")

  estimates = data.frame(steps = rep(steps,ncx+ncv+1),
                         alpha = rep(data.sample$alpha[1],ncx+ncv+1),
                         par_name = par_names,
                         par_real = c(par_real[,1],data.sample$alpha[1]),
                         emv,
                         Method = rep("Method 1",ncx+ncv+1),
                         n = rep(n,ncx+ncv+1),
                         Model = rep(paste0("Model ",model),ncx+ncv+1),
                         value.start = round( value.start,4)
  )

  write.table(
    estimates,
    file.estimates.emv,
    col.names = F,
    row.names = F,
    append = T,
    quote = T,
  )
  setwd(wd)

  setwd(wd.results_hessian)

  file.hessian.alpha = paste0("hessian_",alpha.value,"_n",n, ".txt")
  file.hessian.covBeta = paste0("hessian_Beta_",alpha.value,"_n",n,".txt")
  file.hessian.covGamma = paste0("hessian_Gamma_",alpha.value,"_n",n,".txt")

  hessianAlpha = cbind(steps, emv.ALPHA$hessian)
  hessianCovBeta = cbind(steps, emv.BETA$hessian)
  hessianCovLambda = cbind(steps, emv.LAMBDA$hessian)

  write.table(
    hessianAlpha,
    file.hessian.alpha,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  write.table(
    hessianCovBeta,
    file.hessian.covBeta,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  write.table(
    hessianCovLambda,
    file.hessian.covGamma,
    col.names = F,
    row.names = F,
    quote = T,
    append = T
  )

  setwd(wd)
}
