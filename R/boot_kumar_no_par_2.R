boot_kumar_no_par_2 <- function(idx, n, B, mtrim) {
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
  row.deleted <- (1:(max(q, p)))
  data1 <- results$data %>%
    arrange(Group) %>%
    group_by(City) # %>% slice(-(1:(max(q,p))))
  quantil <- results$quantil



  G = paste0("Group ",sample(c(1,2,3),3,replace=T))
  data <- data1 %>% filter(Group %in% G)%>%arrange(Group,City,Year,Month)

  N <- results$N
  erro <- results$erro
  formula <- results$formula
  w <- data$weights
  link_kl <- results$link.kl
  link_delta <- results$link.delta


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
