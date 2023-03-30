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
formula <- RH ~ sent + cost|semester
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
#emv$par = c(3.06857482, -0.06830916, 1.08626420, -4.89575000, 0.56244984) # m1
emv$par = c(2.973018248,0.250886046,-0.001380885,-8.373751993,-0.309405318,0.676654152 )#m2
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
for (t in 2:(TT + 1)) {
  a[t] = exp(xbeta[t])
}

# geração yt dado yt-1  ----------

yt.y1 = y2.dado.y1 = NULL
nn = 200 # n da unif
tt = 3
at1 = a[tt]
bt1 = b[tt]
at = a[tt + 1]
bt = b[tt + 1]

u1 = runif(nn)
y1 = (1 - exp(-(1 / bt1) * (-log(u1)) ^ (1 / alpha))) ^ (1 / at1)
plot.ts(y1)
abline(h = median(y1))
abline(h = y[tt])
yt.y1[1] <- sample(y1, 1, prob = NULL, replace = T)
yt.y1[1] <- median(y1)

for (i in 1:nn) {
  cat(i, "\n")
  u2 = runif(1)
  int.yt = try (uniroot(
    p.cond,
    interval = c(0, 1),
    u = u2,
    y.t.1 = yt.y1[1],
    at1 = at1,
    at = at,
    bt1 = bt1,
    bt = bt,
    alpha = alpha
  ),
  silent = T)
  if (class(int.yt) == "try-error") {
    next
  }
  else{
    y2.dado.y1[i] <- int.yt$root
  }

}


plot.ts(na.omit(y2.dado.y1))
abline(h = median(y2.dado.y1, na.rm = T))
abline(h = y[tt + 1])


# geração forma geral -----

alphas = c(.35,.5,.65,.80,.95)
#alphas = c(0.7258732 )
n = 144
N = 300
tic <- tictoc::tic()
for (al in 1:length(alphas)) {

  alpha <- alphas[al]
  cat("alpha=", alpha, "\n")
  at1 = a[1]
  bt1 = b[1]
  u1 = runif(N)
  y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha))) ^ (1 / at1)
  path.sample <-
    paste0(
      "scripts_tests/model_time_v2/Block_8/simulation_1/alpha",
      al,
      "_m2",
      ".txt"
    )
  y.aux <- data.frame(y.sim = y1, t = 1, alpha = alpha)
  write.table(
    y.aux,
    path.sample,
    sep = " ",
    append = T,
    quote = T,
    row.names = F,
    col.names = F
  )
  yt.y1 <- array(c(0), dim = c(n, N))
  yt.y1[1, ] <-  y1

  #seq.i = seq(3,5,2)
  for (i in 2:n) {
    #  i=2
    cat(i, "\n")
    at1 = a[i - 1]
    bt1 = b[i - 1]
    at = a[i]
    bt = b[i]
    #u2 <- runif(1,0,1)
    for (k in 1:N) {
      u2 <- runif(1, 0, 1)
      int.yt = try (uniroot(
        p.cond,
        interval = c(0, 1),
        u = 1 - u2,
        y.t.1 = median(yt.y1[i - 1, ], na.rm = T),
        at1 = at1,
        at = at,
        bt1 = bt1,
        bt = bt,
        alpha = alpha
      ),
      silent = T)
      if (class(int.yt) == "try-error") {
        next
      }
      else{
        yt.y1[i, k] <- int.yt$root
      }
    }

    path.sample <- paste0(
      "scripts_tests/model_time_v2/Block_8/simulation_1/alpha",
      al,
      "_m2",
      ".txt"
    )

    y.aux <- data.frame(y.sim = yt.y1[i, ],
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
toc <- tictoc::toc()

dd = 144
plot.ts(na.omit(yt.y1[dd, ]))
abline(h = median(yt.y1[dd, ], na.rm = T))
abline(h = mediana[dd])
abline(h = y[dd])
path.sample.int <-
  paste0("scripts_tests/model_time_v2/Block_7/sample_root_alpha_",
         al,
         "_m2_est1",
         ".txt")
ff = read.table(file = path.sample.int, sep = " ")
t.i = ff %>% filter(V2 == 144) %>% pull(V1) %>% na.omit()
plot.ts(t.i)
abline(h = median(t.i, na.rm = T))
abline(h = mediana[144])
abline(h = y[144])
ff[nrow(ff), ]
41100 / 300
plot.ts(c(ff[1, ]))
mean(c(ff[1, ]), na.rm = T)

ff = read.table(file = path.sample.int, sep = " ")
yy = ff %>% group_by(V2) %>% summarise(m = median(V1)) %>%
  pull(m)
plot.ts(yy)
#y.real = y
lines(yy, col = 2)
acf(abs(yy))
acf(y.real)
