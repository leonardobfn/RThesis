


rm(list = ls())
devtools::load_all() # meu pacote
#devtools::install()
require(RThesis)
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
n. = c(84, 144, 168)
logs = NULL
MODEL = 1
FORMULA <- RH ~ sent + cost | semester
alpha_value <- c ("alpha35",
                  "alpha50",
                  "alpha65",
                  "alpha80",
                  "alpha95")
alpha_value = alpha_value[c(1,2)]
MC = 5
cpus <- 5
wd. = getwd()
require(snowfall)
sfInit(cpus = cpus,
       type = "SOCK",
       parallel = TRUE)
sfExportAll()
#sfLibrary(RThesis)
sfLibrary(tidyr)
sfLibrary(dplyr)
sfLibrary(extraDistr)
sfLibrary(Formula)
for (k in n.) {
  #n. = N.[k] # length of series

  for (alpha. in alpha_value) {
    tic <- tictoc::tic()
    sfLapply(
      1:MC,
      fun = snowfall.simulation,
      model = MODEL,
      alpha.value = alpha.,
      formula = FORMULA,
      wd = wd.,
      n = k
    )
    toc <- tictoc::toc()
    logs.aux = data.frame(
      cpus = cpus,
      n = k,
      alpha = alpha.,
      Model = paste0("Model ", MODEL),
      time = toc$callback_msg
    )
    logs = rbind(logs, logs.aux)
  }
  setwd(wd.)
}
sfStop()
write.table(logs, "Data_simulation/logs.txt")
read.table("Data_simulation/logs.txt")
#399.16 sec elapsed

# b = read.table("Data_simulation/Model_1/simulations/n144/alpha35/data5.txt")
# plot.ts(b$y)
