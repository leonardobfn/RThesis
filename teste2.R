

rm(list = ls())
devtools::load_all() # meu pacote
devtools::install()
require(RThesis)
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
#source("auxiliary_functions.R")
#compiler::enableJIT(3) #
# snowfall.estimates_Method_3 = function(steps,
#                                        model,
#                                        alpha.value,
#                                        erro = 10 ^ (-4),
#                                        formula,
#                                        wd_results,
#                                        wd_sample,
#                                        wd,
#                                        wd_covariates
# )
#   {
#   # model=1
#   #  steps = 650
#   #  alpha.value = "alpha95"
#   #print(MC)
#
#   # setwd(wd_sample)
#   # path.sample = paste0("Model_", model,"/simulations/",alpha.value,"/")
#   # #file.sample = paste0( path.sample, data.labels[steps])
#   # file.sample = paste0( path.sample,"data",steps,".txt")
#   # data.sample = read.table(file.sample)
#   # wd_covariates=wd.covariates
#   # formula = FORMULA
#   # setwd(wd_covariates)
#   covi = read.table("Data_simulation/covariates.txt")
#   covi$semester <- as.factor(covi$semester)
#
#   data = data.frame(covi, RH = rnorm(nrow(covi)))
#   mf <- model.frame(Formula::Formula(formula), data = data)
#   #y <- model.response(mf)
#   cov_a <-
#     model.matrix(Formula::Formula(formula),
#                  data = data,
#                  rhs = 1)
#   cov_delta <-
#     model.matrix(Formula::Formula(formula),
#                  data = data,
#                  rhs = 2)
#
#   ncx <- ncol(cov_a)
#
#   ncv <- ncol(cov_delta)
#   par_names <- c(paste0(colnames(cov_a), "_a"),
#                  paste0(colnames(cov_delta), "_delta"),
#                  "alpha")
#   par_real = read.table("real_par_model_1.txt")
#   beta.real <- par_real[c(1:ncx),1]
#   lambda.real <- par_real[c((1 + ncx):(ncx + ncv)),1]%>% as.matrix()
#   alpha_ger <- 0.95
#
#   a = exp(cov_a%*%beta.real)
#   delta = exp(cov_delta%*%lambda.real)
#   b = 1/delta
#
#   # geração dos dados inicio----------
#   nn = 144
#   yt.y1 = matrix(0, nn, 1)
#
#   at1 = a[1]
#   bt1 = b[1]
#   #u1 = runif(1)
#   repeat {
#     u1 = runif(1)
#     #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
#     y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
#     yt.y1[1] <-  y1
#     if( format( y1 ) != 1 & format( y1 ) != 0)
#       break
#   }
#
#   for (i in 2:nn) {
#     #i=96
#     at1 = a[i - 1]
#     bt1 = b[i - 1]
#     at = a[i]
#     bt = b[i]
#
#     # repeat{
#     #   u2 = runif(1,0,1)
#     #   if(abs(u2-0)>0.02 & abs(u2-1)>0.02)
#     #     break
#     # }
#
#     u2 = runif(1, 0, 1)
#
#     int.yt = try (uniroot(
#       p.cond,
#       interval = c(0,1),
#       u = 1 - u2,
#       y.t.1 = yt.y1[i - 1],
#       at1 = at1,
#       at = at,
#       bt1 = bt1,
#       bt = bt,
#       alpha = alpha_ger
#     ),
#     silent = T)
#
#     if (class(int.yt) == "try-error" & i == 2) {
#       repeat{
#         repeat {
#           u1 = runif(1)
#           #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
#           y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
#           yt.y1[1] <-  y1
#           if( format( y1 ) != 1 & format( y1 ) != 0)
#             break
#         }
#
#         u2 = runif(1, 0, 1)
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 1],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#
#         if (class(int.yt) != "try-error") {
#           break
#         }
#       }
#     }
#     #int.yt = list()
#     #int.yt$root = 0.000000e+00
#     if (class(int.yt) != "try-error" & i == 2 & (format(int.yt$root)==0|format(int.yt$root)==1)) {
#       test = "try-error"
#       while(test == "try-error") {
#         repeat {
#           u1 = runif(1)
#           #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
#           y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
#           yt.y1[1] <-  y1
#           if( format( y1 ) != 1 & format( y1 ) != 0)
#             break
#         }
#
#         u2 = runif(1, 0, 1)
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 1],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#         test = class( int.yt )
#         if(test!="try-error" ){
#           if(format(int.yt$root)!=0 & format(int.yt$root)!=1)
#             break else{
#               test = "try-error"
#             }
#         }
#         # if (class(int.yt) != "try-error") {
#         #   break
#         # }
#       }
#     }
#
#     if (class(int.yt) == "try-error" & i>2) {
#
#       repeat {
#
#         u2 = runif(1, 0, 1)
#         at1 = a[i - 2]
#         bt1 = b[i - 2]
#         at = a[i - 1]
#         bt = b[i - 1]
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 2],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#         if(class(int.yt)=="try-error"){
#           int.yt = list()
#           int.yt$root <- NA
#         }
#         yt.y1[i - 1] <- int.yt$root
#
#         u2 = runif(1, 0, 1)
#         at1 = a[i - 1]
#         bt1 = b[i - 1]
#         at = a[i]
#         bt = b[i]
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 1],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#
#         if (class(int.yt) != "try-error") {
#           break
#         }
#
#       }
#     }
#
#     if (class(int.yt) != "try-error" & i>2 & (format(int.yt$root)==0|format(int.yt$root)==1)) {
#       test = "try-error"
#       while(test == "try-error") {
#
#         u2 = runif(1, 0, 1)
#         at1 = a[i - 2]
#         bt1 = b[i - 2]
#         at = a[i - 1]
#         bt = b[i - 1]
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 2],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#
#         if(class(int.yt)=="try-error"){
#           int.yt = list()
#           int.yt$root <- NA
#         }
#
#         yt.y1[i - 1] <- int.yt$root
#
#         u2 = runif(1, 0, 1)
#         at1 = a[i - 1]
#         bt1 = b[i - 1]
#         at = a[i]
#         bt = b[i]
#
#         int.yt = try (uniroot(
#           p.cond,
#           interval = c(0,1),
#           u = 1 - u2,
#           y.t.1 = yt.y1[i - 1],
#           at1 = at1,
#           at = at,
#           bt1 = bt1,
#           bt = bt,
#           alpha = alpha_ger
#         ),
#         silent = T)
#
#         test = class( int.yt )
#
#         if(test!="try-error" ){
#           if(format(int.yt$root)!=0 & format(int.yt$root)!=1)
#             break else{
#               test = "try-error"
#             }
#         }
#
#         # if (class(int.yt) != "try-error") {
#         #   break
#         # }
#
#       }
#     }
#
#     yt.y1[i] <- int.yt$root
#   }
#   data$RH = yt.y1
#   y = yt.y1
#   # geração dos dados fim ----------
#
#
#
#
#
#   start_aux_betas = try(coef(lm(-log(-log(RH)) ~ sent + cost, data = data)))
#
#   # if (class(start_aux_betas) != "try-error") {
#   #   par.cov.Gbase <-
#   #     try(optim(
#   #       par =  c(start_aux_betas, rep(-.1, ncv)),
#   #       fn = log.kumar.gbase,
#   #       control = list(fnscale = -1),
#   #       method = "BFGS",
#   #       # method = "L-BFGS-B",
#   #       # lower = c(rep(-Inf,6),0.25),
#   #       # upper = c(rep(Inf,6),0.95),
#   #       x = cov_a,
#   #       w = cov_delta,
#   #       y = y
#   #     ),
#   #     silent = T)
#   #
#   #   if (class(par.cov.Gbase) == "try-error") {
#   #     theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
#   #   }else{
#   #     theta.start = par.cov.Gbase$par
#   #   }
#   # } else{
#   #   theta.start = c(0.80 ,0.20,-0.02 , rep(-.1, ncv))
#   # }
#
#
#   theta.start = c(round(start_aux_betas, 4), -.1, -.1, 0.5)
#   value.start = theta.start
#
#   repeat {
#
#     beta <- theta.start[c(1:ncx)]
#     lambda <- theta.start[c((1 + ncx):(ncx + ncv))]
#     alpha <- theta.start[-c((1):(ncx + ncv))]
#
#     emv.beta <-
#       try(optim(
#         par = beta,
#         fn = log.f.cond.beta,
#         control = list(fnscale = -1),
#         method = "BFGS",
#         # method = "L-BFGS-B",
#         # lower = c(rep(-Inf,ncx+ncv),0.25),
#         # upper = c(rep(Inf,ncx+ncv),0.95),
#         x = cov_a,
#         w = cov_delta,
#         y = y,
#         lambda = lambda,
#         alpha = alpha,
#         ncx = ncx,
#         ncv = ncv,
#         hessian = T
#       ),
#       silent = T)
#
#     emv.lambda <-
#       try(optim(
#         par = lambda,
#         fn = log.f.cond.lambda,
#         control = list(fnscale = -1),
#         method = "BFGS",
#         # method = "L-BFGS-B",
#         # lower = c(rep(-Inf,ncx+ncv),0.25),
#         # upper = c(rep(Inf,ncx+ncv),0.95),
#         x = cov_a,
#         w = cov_delta,
#         y = y,
#         beta = emv.beta$par,
#         alpha = alpha,
#         hessian = T
#       ),
#       silent = T)
#
#
#     emv.alpha <-
#       try(optim(
#         par = alpha,
#         fn = log.f.cond.alpha,
#         control = list(fnscale = -1),
#         #method = "BFGS",
#         method = "L-BFGS-B",
#         lower = c(0.1),
#         upper = c(.99),
#         x = cov_a,
#         w = cov_delta,
#         y = y,
#         ncx = ncx,
#         ncv = ncv,
#         theta = c(emv.beta$par, emv.lambda$par),
#         hessian = T
#       ),
#       silent = T)
#
#     if (class(emv.alpha) != "try-error" & emv.alpha$value == 0) {
#       estimates.aux = list()
#       estimates.aux$par = rep(NA, ncx + ncv + 1)
#       break
#     }
#
#     if (class(emv.alpha) == "try-error") {
#       estimates.aux = list()
#       estimates.aux$par = rep(NA, ncx + ncv + 1)
#       break
#     } else{
#       theta.up <- c(emv.beta$par, emv.lambda$par, emv.alpha$par)
#     }
#
#     crit <- sum(((theta.up - theta.start) / theta.start) ^ 2)
#
#     if (crit < erro) {
#       estimates.aux = list()
#       estimates.aux$par = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
#       break
#     } else{
#       theta.start = c(emv.beta$par, emv.lambda$par, emv.alpha$par)
#     }
#
#   }
#
#   if (length(which(is.na(estimates.aux$par) == T)) == 0) {
#     emv = list()
#     emv <- round(estimates.aux$par, 4)
#     emv.BETA = list()
#     emv.BETA$hessian = round(emv.beta$hessian, 4)
#     emv.LAMBDA = list()
#     emv.LAMBDA$hessian = round(emv.lambda$hessian, 4)
#     emv.ALPHA = list()
#     emv.ALPHA$hessian = round(emv.alpha$hessian, 4)
#   } else{
#     emv = list()
#     emv <- estimates.aux$par
#     emv.BETA = list()
#     emv.BETA$hessian = matrix(0, ncx, ncx)
#     emv.LAMBDA = list()
#     emv.LAMBDA$hessian = matrix(0, ncv, ncv)
#     emv.ALPHA = list()
#     emv.ALPHA$hessian = 0
#   }
#
#   #setwd(wd_results)
#   estimates = data.frame(steps, emv, par_names)
#   file.sample  = paste0("resultados/estimates_",alpha.value,".txt")
#   write.table(
#     estimates,
#     file.sample,
#     col.names = F,
#     row.names = F,
#     append = T,
#     quote = T,
#   )
#   #setwd(wd)
#
#   # hessianCovBeta = cbind(steps, emv.BETA$hessian)
#   # write.table(
#   #   hessianCovBeta,
#   #   path_hessian_covBeta,
#   #   col.names = F,
#   #   row.names = F,
#   #   quote = T,
#   #   append = T
#   # )
#   #
#   #
#   # hessianCovLambda = cbind(steps, emv.LAMBDA$hessian)
#   # write.table(
#   #   hessianCovLambda,
#   #   path_hessian_covLambda,
#   #   col.names = F,
#   #   row.names = F,
#   #   quote = T,
#   #   append = T
#   # )
#   #
#   #
#   # hessianAlpha = cbind(steps, emv.ALPHA$hessian)
#   # write.table(
#   #   hessianAlpha,
#   #   path_hessian_alpha,
#   #   col.names = F,
#   #   row.names = F,
#   #   quote = T,
#   #   append = T
#   # )
# }


MODEL = 1
n = 144 # length of series
MC = 1000
cpus <- 5
#ncv  = 2
#ncx = 3
FORMULA <- RH ~ sent + cost | semester

alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")

#data.labels = paste0("/data",1:MC,".txt")
alpha_value = alpha_value[5]
wd. = getwd()
wd.sample = paste0(wd.,"/Data_simulation")
#wd_results = paste0(wd,"Data_simulation/Model_1/estimates/Method_3/")
wd.results = paste0(wd.,"/resultados")
wd.covariates = paste0(wd.,"/Data_simulation/")

# path_hessian_alpha = paste0(
#   "Data_simulation/Model_",
#   model,
#   "/estimates/Method_3/hessian/",
#   "hessian_",
#   alpha.value,
#   ".txt"
# )
# path_hessian_covLambda = paste0(
#   "Data_simulation/Model_",
#   model,
#   "/estimates/Method_3/hessian/",
#   "hessianCovLambda_",
#   alpha.value,
#   ".txt"
# )
# path_hessian_covBeta = paste0(
#   "Data_simulation/Model_",
#   model,
#   "/estimates/Method_3/hessian/",
#   "hessianCovBeta_",
#   alpha.value,
#   ".txt"
# )

erro = 10 ^ (-4)
require(snowfall)

sfInit(cpus = cpus,
       type = "SOCK",
       parallel = TRUE)
sfExportAll()
sfLibrary(RThesis)
sfLibrary(tidyr)
sfLibrary(dplyr)
sfLibrary(extraDistr)
sfLibrary(Formula)
#sfClusterSetupRNG(seed=7221)
tic <- tictoc::tic()
#for (alpha. in alphas) {
#cat(alpha., "\r")
# alpha. = 0.95

# sfClusterApplyLB(1:MC,fun = snowfall.estimates_Method_3,
#          model = MODEL,
#          alpha.value = alpha_value,
#          formula = FORMULA,
#          wd_results = wd.results,
#          wd_sample = wd.sample,
#          wd=wd,
#          wd_covariates = wd.covariates
#)

sfLapply(1:MC,fun = snowfall.estimates_Method_3,
         model = MODEL,
         alpha.value = alpha_value,
         formula = FORMULA,
         wd_results = wd.results,
         wd_sample = wd.sample,
         wd=wd.,
         wd_covariates = wd.covariates
)

#}
toc <- tictoc::toc()
sfStop()
