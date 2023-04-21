

rm(list = ls())
devtools::load_all() # meu pacote
#devtools::install()
require(RThesis)
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)


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

#alpha_value = alpha_value[5]
wd. = getwd()
wd.sample = paste0(wd.,"/Data_simulation/Model_",MODEL,"/simulations/n",n.)
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

  sfLapply(1:MC,fun = snowfall.simulation,
           model = MODEL,
           alpha.value = alpha.,
           formula = FORMULA,
           wd_sample = wd.sample,
           wd=wd.,
           wd_covariates = wd.covariates,
           n=n.
  )

}
toc <- tictoc::toc()
sfStop()
#399.16 sec elapsed
