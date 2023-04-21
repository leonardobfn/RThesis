snowfall.simulation = function(steps,
                               model,
                               alpha.value,
                               formula,
                               wd_sample,
                               wd,
                               wd_covariates,
                               n)
{
  # model=1
  #  steps = 650
  #  alpha.value = "alpha95"
  #print(MC)

  # setwd(wd_sample)
  # path.sample = paste0("Model_", model,"/simulations/",alpha.value,"/")
  # #file.sample = paste0( path.sample, data.labels[steps])
  # file.sample = paste0( path.sample,"data",steps,".txt")
  # data.sample = read.table(file.sample)
  # wd_covariates=wd.covariates
  #wd_sample = wd.sample
  # formula = FORMULA
  #n = 216
  setwd(wd_covariates)
  covi = read.table("covariates.txt")[1:n,]

  covi$semester <- as.factor(covi$semester)

  data = data.frame(covi, RH = rnorm(nrow(covi)))
  mf <- model.frame(Formula::Formula(formula), data = data)
  #y <- model.response(mf)
  cov_a <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 1)
  cov_delta <-
    model.matrix(Formula::Formula(formula),
                 data = data,
                 rhs = 2)

  ncx <- ncol(cov_a)

  ncv <- ncol(cov_delta)
  par_names <- c(paste0(colnames(cov_a), "_a"),
                 paste0(colnames(cov_delta), "_delta"),
                 "alpha")

  par_real = read.table(paste0("real_par_model_", model, ".txt"))
  beta.real <- par_real[c(1:ncx), 1]
  lambda.real <-
    par_real[c((1 + ncx):(ncx + ncv)), 1] %>% as.matrix()
  alpha_ger <-
    paste0("0.", stringr::str_sub(alpha.value, 6, 7)) %>% as.numeric()

  a = exp(cov_a %*% beta.real)
  delta = exp(cov_delta %*% lambda.real)
  b = 1 / delta

  # geração dos dados inicio----------
  nn = n
  yt.y1 = matrix(0, nn, 1)

  at1 = a[1]
  bt1 = b[1]
  #u1 = runif(1)
  repeat {
    u1 = runif(1)
    #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
    y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
    yt.y1[1] <-  y1
    if (format(y1) != 1 & format(y1) != 0)
      break
  }

  for (i in 2:nn) {
    #i = 84
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
      alpha = alpha_ger
    ),
    silent = T)

    if (class(int.yt) == "try-error" & i == 2) {
      repeat {
        repeat {
          u1 = runif(1)
          #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
          y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
          yt.y1[1] <-  y1
          if (format(y1) != 1 & format(y1) != 0)
            break
        }

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
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }
      }
    }
    #int.yt = list()
    #int.yt$root = 0.000000e+00
    if (class(int.yt) == "try-error") {
      int.yt = list()
      int.yt$root = 0
      test2 = "try-error"
    } else{
      test2 = "no-error"
    }

    if (class(int.yt) != "try-error" &
        i == 2 &
        (format(int.yt$root) == 0 | format(int.yt$root) == 1)) {
      test = "try-error"
      while (test == "try-error") {
        repeat {
          u1 = runif(1)
          #G = (1 - .8 ^ (1 / b[1])) ^ a[1]
          y1 = (1 - exp(-(1 / bt1) * (-log(1 - u1)) ^ (1 / alpha_ger))) ^ (1 / at1)
          yt.y1[1] <-  y1
          if (format(y1) != 1 & format(y1) != 0)
            break
        }

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
          alpha = alpha_ger
        ),
        silent = T)
        test = class(int.yt)
        if (test != "try-error") {
          if (format(int.yt$root) != 0 & format(int.yt$root) != 1)
            break
          else{
            test = "try-error"
          }
        }
        # if (class(int.yt) != "try-error") {
        #   break
        # }
      }
    }

    if ((class(int.yt) == "try-error" |
         test2 == "try-error") & i > 2) {
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
          alpha = alpha_ger
        ),
        silent = T)
        if (class(int.yt) == "try-error") {
          int.yt = list()
          int.yt$root <- NA

        }

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
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) != "try-error") {
          break
        }
      }
    }

    if (class(int.yt) != "try-error" &
        i > 2 &
        (format(int.yt$root) == 0 | format(int.yt$root) == 1)) {
      test = "try-error"
      while (test == "try-error") {
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
          alpha = alpha_ger
        ),
        silent = T)

        if (class(int.yt) == "try-error") {
          int.yt = list()
          int.yt$root <- NA
        }

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
          alpha = alpha_ger
        ),
        silent = T)

        test = class(int.yt)

        if (test != "try-error") {
          if (format(int.yt$root) != 0 & format(int.yt$root) != 1)
            break
          else{
            test = "try-error"
          }
        }

        # if (class(int.yt) != "try-error") {
        #   break
        # }

      }
    }

    yt.y1[i] <- int.yt$root
  }

  y = yt.y1
  setwd(wd)
  setwd(wd_sample)
  file.sample  = paste0(alpha.value, "/data", steps, ".txt")

  y.data = data.frame(
    y = y,
    alpha = rep(alpha_ger, length(y)),
    Model = rep(paste0("Model ", model), length(y))
  )
  #saveRDS(y.data, file.sample)
  write.table(y.data, file.sample)
  setwd(wd)
}
