

rm(list = ls())
#devtools::load_all() # meu pacote
#devtools::install()
#require(RThesis)
#source("auxiliary_functions.R")
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)

METHOD = 2
MODEL = 1
n. = 216 # length of series
MC = 1200
cpus <- 10
FORMULA <- RH ~ sent + cost | semester

alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")

#data.labels = paste0("/data",1:MC,".txt")
alpha_value = alpha_value[2]
wd. = getwd()
wd.sample = paste0(wd.,"/Data_simulation/Model_",MODEL,"/simulations/n",n.)
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

  sfLapply(1:MC,fun = snowfall.estimates_Method_2.mpfr,
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
