rm(list = ls())

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
source("auxiliary_functions.R")
compiler::enableJIT(3)

# Loading database------
#devtools::load_all() # loading my functions

data("data_1")
head(data_1)

#Parse database------

month_names <-
  factor(
    month.abb,
    levels = c(
      "Jan",
      "Feb",
      "Mar",
      "Apr",
      "May",
      "Jun",
      "Jul",
      "Aug",
      "Sep",
      "Oct",
      "Nov",
      "Dec"
    )
  )

database <- data_1 %>%
  filter(Station == "1") %>%
  mutate(
    semester = rep(c(1, 2), each = 6) %>% rep(21) %>% rep(14) %>% as.factor(),
    month_names = rep(month_names, 21) %>% rep(14),
    date = paste0(Year, "/", month_names),
    t = seq(1, 252) %>% rep(14),
    Group = rep(paste0("Group ", c(
      1, 2, 3, 3, 2, 1, 1, 3, 2, 3, 3, 3, 1, 1
    )), each = 252),
    cost = cos(2 * pi * as.numeric(Month) / 12),
    sent = sin(2 * pi * as.numeric(Month) / 12),
    lles = (log(10 ^ ((7.5 * TBS) / (237.3 + TBS)
    )))
  )

head(database, 13)
citys <- database$City %>% unique()
results <- results.mar <- matrix(0, length(citys), 8)
TT = 143
y.data <- matrix(0, length(citys), TT + 1)
tic <- tictoc::tic()
p = 10 # manaus
data <-
  database %>% filter(City == citys[p]) %>% slice((252 - TT):252)
data$precp[data$precp == 0] <- 1
#formula <- RH ~ lles | lles - 1
formula <- RH ~ sent + cost | semester
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 2)
ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)
par_names <- c(paste0(colnames(cov_a),"_a"),
               paste0(colnames(cov_delta),"_delta"))
# values real paramenters betas and lambdas -----------

emv <- list()
emv$par = c(2.973018248,
            0.250886046,
            -0.001380885,
            -8.373751993,
            -0.309405318)
# emvs <- data.frame(real_values=emv$par,par_names = par_names)
# write.table(emvs,"Data_simulation/real_par_model_1.txt",col.names = T,sep=" ")
# read.table("Data_simulation/real_par_model_1.txt")
theta = emv$par
Beta = theta[1:ncx] # real
lambda = theta[(ncx + 1):c(ncx + ncv)] # real
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = exp(xbeta)

# simulations-----


library(snowfall)
alphas = c(0.35,0.50, 0.65, 0.80, 0.95) # values real alpha
model = 1
n = 144 # length of series
MC = 1000
cpus <- 4
sfInit(parallel = TRUE, cpus = cpus)
sfExportAll()
tic <- tictoc::tic()
for (alpha. in alphas) {
  cat(alpha., "\r")
  #sfClusterApplyLB(1:MC,model=model,snowfall.simulation,nn=n,a=a,b=b,alpha=alpha.)
  sfLapply(
    1:MC,
    model = model,
    snowfall.simulation,
    nn = n,
    a = a,
    b = b,
    alpha = alpha.
  )
}
toc <- tictoc::toc()
sfStop()
#934.99 sec elapsed
