rm(list = ls())

#
dir.create("Data_simulation/Model_1/simulations/alpha50")
dir.create("Data_simulation/Model_1/simulations/alpha65")
dir.create("Data_simulation/Model_1/simulations/alpha80")
dir.create("Data_simulation/Model_1/simulations/alpha95")

# packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(ggplot2)
source("scripts_tests/model_time_v2/Block_8/auxiliary_functions.R")
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

# values real paramenters betas and lambdas -----------

emv <- list()
emv$par = c(2.973018248,
            0.250886046,
            -0.001380885,
            -8.373751993,
            -0.309405318)
theta = emv$par
Beta = theta[1:ncx] # real
lambda = theta[(ncx + 1):c(ncx + ncv)] # real
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = exp(xbeta)

# simulations-----


snowfall.simulation <- function(nn, a, b, alpha, model, MC) {
  #id=1
  # alpha = 0.5
  # alpha = 0.8

  alpha.value <- switch (
    as.character(alpha),
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

  yt.y1 = matrix(0, nn, 1)

  at1 = a[1]
  bt1 = b[1]
  u1 = runif(1)
  G = (1-.8^(1/b[1]))^a[1]
  y1.. = exp(-  (-log(G))^alpha )
  u1 = runif(1,0,y1..)
  y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha))) ^ (1 / at1)
  yt.y1[1] <-  y1

  y.aux <- data.frame(y.sim = yt.y1[1],
                      t = 1,
                      alpha = alpha)

  write.table(
    y.aux,
    path.sample,
    sep = " ",
    append = T,
    quote = T,
    row.names = F,
    col.names = F
  )

  for (i in 2:nn) {
    at1 = a[i - 1]
    bt1 = b[i - 1]
    at = a[i]
    bt = b[i]

    # repeat{
    #   u2 = runif(1,0,1)
    #   if(abs(u2-0)>0.02 & abs(u2-1)>0.02)
    #     break
    # }

    u2 = runif(1, 0, 1)

    int.yt = try (uniroot(
      p.cond,
      interval = c(0, 1),
      u = 1 - u2,
      y.t.1 = yt.y1[i - 1],
      at1 = at1,
      at = at,
      bt1 = bt1,
      bt = bt,
      alpha = alpha
    ),
    silent = T)

    if (class(int.yt) == "try-error") {
      repeat {
        u2 = runif(1, 0, 1)
        at1 = a[i - 2]
        bt1 = b[i - 2]
        at = a[i - 1]
        bt = b[i - 1]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 2],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha
        ),
        silent = T)

        yt.y1[i - 1] <- int.yt$root

        u2 = runif(1, 0, 1)
        at1 = a[i - 1]
        bt1 = b[i - 1]
        at = a[i]
        bt = b[i]

        int.yt = try (uniroot(
          p.cond,
          interval = c(0, 1),
          u = 1 - u2,
          y.t.1 = yt.y1[i - 1],
          at1 = at1,
          at = at,
          bt1 = bt1,
          bt = bt,
          alpha = alpha
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }

      }
    }

    yt.y1[i] <- int.yt$root

    y.aux <- data.frame(y.sim = yt.y1[i],
                        t = i,
                        alpha = alpha)

    write.table(
      y.aux,
      path.sample,
      sep = " ",
      append = T,
      quote = T,
      row.names = F,
      col.names = F
    )
  }

}

library(snowfall)
alphas = c(0.50, 0.65, 0.80, 0.95) # values real alpha
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
    alpha = 0.5
  )
}
toc <- tictoc::toc()
sfStop()
