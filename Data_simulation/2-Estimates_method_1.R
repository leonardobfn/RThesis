

rm(list = ls())
#devtools::load_all() # meu pacote
#devtools::install()
require(RThesis)
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)

METHOD = 1
MODEL = 1
n. = 144 # length of series
MC = 1200
cpus <- 10
#ncv  = 2
#ncx = 3
FORMULA <- RH ~ sent + cost | semester

alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")

#data.labels = paste0("/data",1:MC,".txt")
alpha_value = alpha_value[2]
wd. = getwd()
#wd.sample = paste0("G:/Meu Drive/data/Model_1/n216/")
wd.sample = paste0(wd.,"/Data_simulation/Model_",MODEL,"/simulations/n",n.)
#wd_results = paste0(wd,"Data_simulation/Model_1/estimates/Method_3/")
wd.results_emv = paste0(wd.,"/Data_simulation/Model_",MODEL,"/estimates/Method_",METHOD,"/estimates/")
wd.results_hessian = paste0(wd.,"/Data_simulation/Model_",MODEL,"/estimates/Method_",METHOD,"/hessian/")
wd.covariates = paste0(wd.,"/Data_simulation/")



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
for (alpha. in alpha_value) {
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

  sfLapply(1:MC,fun = snowfall.estimates_Method_1,
           model = MODEL,
           alpha.value = alpha.,
           formula = FORMULA,
           wd_sample = wd.sample,
           wd=wd.,
           wd_covariates = wd.covariates,
           n=n.,
           wd_results_emv=wd.results_emv,
           wd_results_hessian=wd.results_hessian
  )

}
toc <- tictoc::toc()
sfStop()
#399.16 sec elapsed
