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
formula <- RH ~ lles|lles-1
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
#emv$par = c(3.06857482, -0.06830916, 1.08626420, -4.89575000, 0.56244984) # m1
#emv$par = c(3.06857482, -0.06830916, -4.89575000, 0.56244984) # m1
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
# geração forma geral -----
alphas = c(.35,.5,.65,.80,.95)
alphas = c(.5,.65,.80,.95)
model=1
n = 144
N = 1
u2 = NULL
a.r = c(0,0.05,.10,.20,.30)
a.r = c(0)
B = 100
tic <- tictoc::tic()
for (bb in 1:B){
  for (aa.r in 1:length(a.r)) {
    #aa.r = 1
    arr = aa.r
    ar. = a.r[aa.r]
    for (al in 1:length(alphas)) {
      #al=1
      alpha <- alphas[al]
      cat("alpha=", alpha, "\n")
      at1 = a[1]
      bt1 = b[1]
      u1 = runif(N)
      # u2[1] = u1
      y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha))) ^ (1 / at1)
      path.sample <-
        paste0(
          "scripts_tests/model_time_v2/Block_8/simulation_4/model_1/",bb,"_alpha",
          al,
          "ar",
          arr,
          "_m",
          model,
          ".txt"
        )
      y.aux <- data.frame(
        y.sim = y1,
        t = 1,
        alpha = alpha,
        ar = ar.,
        u=u1
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
      yt.y1 <- array(c(0), dim = c(n, N))
      yt.y1[1,] <-  y1

      for (i in 2:n) {
        #i=64
        cat("i=", i, "\n")
        #cat(i, "\n")

        a[i] = exp(xbeta[i] +
                     ar. * (-log(-log(yt.y1[i - 1])) - (xbeta[i-1]-emv$par[1])  )   )

        at1 = a[i - 1]
        bt1 = b[i - 1]
        at = a[i]
        bt = b[i]


        for (k in 1:N) {
          #k=1
          #i=4
          # repeat{
          #   u2 = runif(1,.10,1)
          #   if(abs(u2-.10)>0.02 & abs(u2-1)>0.02)
          #     break
          # }
          u2 = runif(1,0,1)
          #repeat {
          int.yt = try (uniroot(
            p.cond,
            interval = c(0,1),
            u = 1 - u2,
            y.t.1 = median(yt.y1[i - 1,], na.rm = T),#sample(na.omit(yt.y1[i - 1,]),1,prob=NULL,replace=T),
            at1 = at1,
            at = at,
            bt1 = bt1,
            bt = bt,
            alpha = alpha
          ),
          silent = T)
          # if (class(int.yt) != "try-error") {
          #   yt.y1[i, k] <- int.yt$root
          #   break
          # }
          if (class(int.yt) == "try-error") {
            #yt.y1[i, k] <- int.yt$root
            next
          }

          # if (format(int.yt$root, scientific = FALSE) != 1) {
          #   break
          # }
          # else{
          #   yt.y1[i, k] <- int.yt$root
          #}
          yt.y1[i, k] <- int.yt$root
        }
        path.sample <-
          paste0(
            "scripts_tests/model_time_v2/Block_8/simulation_4/model_1/",bb,"_alpha",
            al,
            "ar",
            arr,
            "_m",
            model,
            ".txt"
          )

        y.aux <- data.frame(
          y.sim = yt.y1[i,],
          t = i,
          alpha = alpha,
          ar = ar.,
          u=u2

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
}
toc <- tictoc::toc()

path.sample = "scripts_tests/model_time_v2/Block_8/simulation_4/alpha3ar2_m1.txt"
ff = read.table(file = path.sample, sep = " ")
t.i = ff %>% pull(V1)
plot.ts(y)
abline(h = median(t.i, na.rm = T))
lines(t.i,col=3)
