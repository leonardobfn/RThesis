rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
source("auxiliary_functions.R")
compiler::enableJIT(3)



snowfall.estimates = function(MC, model, alpha){
  # MC = 1
  # alpha = 0.95
  # print(MC)
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
         MC,
         ".txt")


s = read.table(path.sample)


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
  browser()
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

estimates.aux <-
  try(optim(
    par = .50,
    fn = log.f.cond.alpha,
    control = list(fnscale = -1),
    #method = "BFGS",
     method = "L-BFGS-B",
     lower = c(0),
     upper = c(0.99),
    x = cov_a,
    w = cov_delta,
    y = y,
    theta = par.cov.Gbase$par
  ),
  silent = T)

alpha.up <- estimates.aux$par
par.covariates.start <- start_aux
erro = 10 ^ (-4)
repeat {
  # Step E-----
  v <- V(theta = par.covariates.start,
         x = cov_a,
         w = cov_delta,
         y = y)


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
  "/estimates/Method_1/estimates/",
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
  "/estimates/Method_1/hessian/",
  MC,
  "_hessianCov_",
  alpha.value,
  ".txt"
)

hessianCov = cbind(1, par.covariates.up$hessian)
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
  "/estimates/Method_1/hessian/",
  MC,
  "_hessian_",
  alpha.value,
  ".txt"
)

hessianAlpha = cbind(1, estimates.aux$hessian)
write.table(
  hessianAlpha,
  path_hessian_alpha,
  col.names = F,
  row.names = F,
  quote = T,
  append = T
)
}
covi = read.table("Data_simulation/covariates.txt")
library(snowfall)
alphas = c(0.35,0.50, 0.65, 0.80, 0.95) # values real alpha
model = 1
n = 144 # length of series
MC = 2
cpus <- 4
sfInit(parallel = TRUE, cpus = cpus)
sfExportAll()
tic <- tictoc::tic()
for (alpha. in alphas) {
  cat(alpha., "\r")
  alpha.=0.95
  sfClusterApplyLB( 1:MC,
                    model = model,
                    fun=snowfall.estimates,
                    alpha = alpha.)
  sfLapply(
    1:MC,
    model = model,
    fun=snowfall.estimates,
    alpha = alpha.
  )
}
toc <- tictoc::toc()
sfStop()
