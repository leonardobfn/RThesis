boot_kumar_no_par_3 <- function(idx, n, B, mtrim) {
  # rm(list=ls())
  # require(extraDistr)
  # require(tidyverse)
  # n <- 1000
  # B <- 1
  # mtrim <- .20
  # devtools::load_all()
  path_estimation <- "estimations_predicts/estimations.rds"
  results <- readRDS(path_estimation)

  p <- results$p
  q <- results$q
  link_kl <- results$link.kl
  link_delta <- results$link.delta

  row.deleted <- (1:(max(q, p)))
  data1 <- results$data %>%
    arrange(Group) %>%
    group_by(City) # %>% slice(-(1:(max(q,p))))

  quantil <- results$quantil
  data1$id <- seq(1, nrow(data1))
  row.deleted <- data1 %>%
    group_by(Group, City) %>%
    slice((1:(max(q, p)))) %>%
    data.frame() %>%
    select(id) %>%
    unlist()

  erro.aux = link_kl$linkfun(results$DADOS$kl_t) -  link_kl$linkfun(results$DADOS$y)
  eta.kl1 <- results$cov_kl[row.deleted,] %*% results$kl[,1]
  erro1.aux <-  eta.kl1 - link_kl$linkfun(data1[row.deleted, "RH"]%>%pull(RH))
  erro <-  klt <- matrix(0,nrow(data1),1)
  klt[-row.deleted,] <- link_kl$linkfun(results$DADOS$kl_t)
  klt[row.deleted,] <- eta.kl1
  erro[-row.deleted,] <-  erro.aux
  erro[row.deleted,] <-  erro1.aux

  blocos <- seq(1:9) %>% rep(each=12) %>% rep(14)
  data2 <- data1 %>% data.frame %>% mutate(erro=erro,blocos=blocos,eta.klt = klt)



  G = paste0("Group ",sample(c(1,2,3),3,replace=T))
  data3 <- data2 %>% filter(Group %in% G)%>%arrange(Group,City,Year,Month)
  data = NULL
  for( i in 1:9){
  blocos.sample <- sample(1:9,1,replace = T)
  aux = data3 %>% group_by(Group,City) %>% filter(blocos %in% blocos.sample)
  data <- rbind(data,aux)
  }
  data = data %>% mutate(RH=link_kl$linkinv(eta.klt+erro))
  data <- data %>% arrange(Group,City)


  N <- results$N
  erro <- results$erro
  formula <- results$formula
  w <- data$weights



  mod <- try(gev_kuma(
    formula = formula, link.kl = link_kl$name, link.delta = link_delta$name, w = w, data = data, N = N,
    grupos = data$Group, quantil = quantil, erro = erro, idx = data$City, p = p, q = q
  ), TRUE)

  if (class(mod) == "try-error") {
    ar_est <- rep(NA, nrow(results$ar))
    betas_est <- rep(NA, nrow(results$kl))
    lambdas_est <- rep(NA, nrow(results$delta))
    alfas_est <- rep(NA, nrow(results$alpha))
  }
  if (class(mod) != "try-error") {
    ar_est <- mod$ar
    betas_est <- mod$kl
    lambdas_est <- mod$delta
    alfas_est <- mod$alpha
  }
  names_alpha <- str_extract(levels(factor(data$Group)),"[0-9]")
  names <- c(
    results$ar %>% rownames(),
    paste0(rownames(mod$kl), "_kl"), paste0(rownames(mod$delta), "_delta"), paste0("alpha_",names_alpha)
  )


  estimation <- data.frame(Estimation = unlist(c(ar_est, betas_est, lambdas_est, alfas_est)), Par = names)
  path_monte_carlo <- "estimations_predicts/estimation_mc.txt"
  write.table(estimation, file = path_monte_carlo, append = T, quote = T, col.names = F)

  log.like <- data.frame(log_like=mod$objective)
  path_log.like <- "estimations_predicts/log_like.txt"
  write.table(log.like, file = path_log.like, append = T, quote = T, col.names = F)

  # return(estimation)
}
