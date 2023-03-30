rm(list = ls())

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
data <- database %>% filter(City == citys[p]) %>% slice((252 - TT):252)
data$precp[data$precp == 0] <- 1
formula <- RH ~ lles|lles
#formula <- RH ~ sent + cost|semester
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <-
  model.matrix(Formula::Formula(formula), data = data, rhs = 2)
ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)
par.names <- c(paste0(colnames(cov_a), "_a"),
               paste0(colnames(cov_delta), "_delta"),
               "alpha")

# values paramenters -----------
emv <- list()
emv$par = c(3.06857482, -0.06830916, 1.08626420, -4.89575000, 0.56244984) # m1
#emv$par = c(2.973018248,0.250886046,-0.001380885,-8.373751993,-0.309405318,0.676654152 )#m2
theta = emv$par
Beta = theta[1:ncx] # real
lambda = theta[(ncx + 1):c(ncx + ncv)] # real
alpha = theta[-c(1:c(ncx + ncv))] # real
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL
a[1] = exp(xbeta[1])

# geração forma geral -----
alphas = c(.35,.5,.65,.80,.95)
model=1
n = 144
N = 1
a.r = c(.10,.20,.30);
tic <- tictoc::tic()
for (aa.r in 1:length(a.r)) {
  arr = aa.r
  ar. = a.r[aa.r]
  for (al in 1:length(alphas)) {
    alpha <- alphas[al]
    cat("alpha=", alpha, "\n")
    yt.y1 <- array(c(0), dim = c(n, N))
    yt.y1[1,] = (pweibull(
      -log(1 - runif(N)),
      shape = 1 / alpha,
      scale = b[1] ^ (alpha),
      lower.tail = TRUE,
      log.p = FALSE
    )) ^ (1 / a[1])

    path.sample <-
      paste0(
        "scripts_tests/model_time_v2/Block_8/simulation_3/alpha",
        al,
        "ar",
        arr,
        "_m",
        model,
        ".txt"
      )

    y.aux <- data.frame(
      y.sim = yt.y1[1],
      t = 1,
      alpha = alpha,
      ar = ar.
    )
    write.table(
      y.aux,
      path.sample,
      sep = " ",
      append = T,
      quote = T,
      row.names = F,
      col.names = F
    )

    for (i in 2:n) {
      cat("i=", i, "\n")

      a[i] = exp(xbeta[i] +
                   ar. * (-log(-log(yt.y1[i - 1]))))

      for (k in 1:N) {
        repeat {
          u = runif(1)
          yt.y1[i, k] = (pweibull(
            -log(1 - u),
            shape = 1 / alpha,
            scale = b[i] ^ (alpha),
            lower.tail = TRUE,
            log.p = FALSE
          )) ^ (1 / a[i])

          if (format(yt.y1[i], scientific = FALSE) != 1) {
            break
          }
        }
      }

      path.sample <-
        paste0(
          "scripts_tests/model_time_v2/Block_8/simulation_3/alpha",
          al,
          "ar",
          arr,
          "_m",
          model,
          ".txt"
        )

      y.aux <- data.frame(
        y.sim = yt.y1[i, ],
        t = i,
        alpha = alpha,
        ar = ar.
      )
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
}
toc <- tictoc::toc()

path.sample = "scripts_tests/model_time_v2/Block_8/simulation_3/alpha2ar1_m1.txt"
ff = read.table(file = path.sample, sep = " ")
t.i = ff %>% pull(V1)
plot.ts(y)
abline(h = median(t.i, na.rm = T))
lines(t.i,col=3)
