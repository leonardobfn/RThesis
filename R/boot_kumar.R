boot_kumar <- function(idx, n, B, mtrim) {
  rm(list=ls())
  require(extraDistr)
  require(tidyverse)
  n <- 1000

  B <- 500
  mtrim <- 0.20
  devtools::load_all()
  path_estimation <- "estimations_predicts/estimations_m6.rds"

  B <- 1
  mtrim <- .20
  devtools::load_all()
  path_estimation <- "estimations_predicts/estimations.rds"

  results <- readRDS(path_estimation)
  link.kl <-results$link.kl$linkfun
  link.eta <-results$link.kl$linkinv
  p <- results$p
  q <- results$q
  row.deleted <- (1:(max(q, p)))
  data <- results$data %>%
    arrange(Group) %>%
    group_by(City) # %>% slice(-(1:(max(q,p))))
  quantil <- results$quantil

  cov_kl <- results$cov_kl
  cov_delta <- results$cov_delta
  ncx <- ncol(cov_kl)
  ncv <- ncol(cov_delta)
  alfa.est <- results$alpha$Estimate

  est <- results$DADOS %>%
    data.frame() %>%
    mutate(
      bq_t = 1 / delta,
      phi_ts = log(kl_t) / log(1 - (1 - quantil)^(delta)),
      aq_t = 1 / phi_ts
    ) %>%
    arrange(grupos)
  Lm <- est %>% group_by(grupos) %>% summarise(l=length(grupos)) %>% pull(l)

  yest.aux <- matrix(0, nrow = nrow(est), ncol = B)
  #qest.aux <- matrix(0, nrow = nrow(est), ncol = k)


#for (l in 1:k){
  for (b in 1:B) {
    estt <- NULL
    for (j in 1:length(alfa.est)) {
      yest.full = NULL

      # qest=NULL
      G <- paste("Group", j)
      alfa.aux <- alfa.est[j]
      zm = rstable_pos(1000,alfa.aux) %>% mean(trim=mtrim)
      #set.seed(10)
      ee <- est %>% dplyr::filter(grupos==G) %>% mutate(
        H=zm*weights,
        #y.eta = rgum(N = Lm[j],mi = log(-log(1-quantil^delta))+link.kl(kl_t),sigma = 1,A = H/delta),
        #yest.res = link.eta(y.eta)
         #eps = rgum(N = Lm[j],mi = log(-log(1-quantil^delta)),sigma = 1,A = H/delta),
         #yest.res = link.eta(link.kl(kl_t) + eps)
        #
        #Espy_z= bq_t*H*beta(1+1/aq_t,bq_t*H),
        #Vary_z= bq_t*H*beta(1+2/aq_t,bq_t*H)-(bq_t*H*beta(1+1/aq_t,bq_t*H))^2,
        betaa=rbeta(1,H,1),
        yest.res = (1-betaa^(delta))^(phi_ts)

      )


      # H <- unique(ee$H)
      # city.names <- ee$idx %>% unique
      # for(k in 1:length(H)){
      #   betaa = rbeta(1,H[k],1)
      #   delta <- ee %>% filter(idx==city.names[k]) %>% pull(delta)
      #   phi_ts <- ee %>% filter(idx==city.names[k]) %>% pull(phi_ts)
      #   klt <- ee %>% filter(idx==city.names[k]) %>% pull(kl_t)
      #   eps = rgum(107,log(-log(1-quantil^delta)),1,H[k]/delta)
      #   yest = link.eta(link.kl(klt) + eps)
      #   #aq_t <- ee %>% filter(idx==city.names[k]) %>% pull(aq_t)
      #   #bq_t <- ee %>% filter(idx==city.names[k]) %>% pull(bq_t)
      #   #yest = (1-betaa^(delta))^(phi_ts)
      #   #yest = rkumar(107,aq_t,bq_t*H[k]) %>% t() %>% t()
      #   yest.full <- rbind(yest.full,yest)
      # }
      #
      estt.aux <- ee %>% mutate(yest = yest.res)

      H <- unique(ee$H)
      city.names <- ee$idx %>% unique
      for(k in 1:length(H)){
        # betaa = rbeta(1,H[k],1)
        # delta <- ee %>% filter(idx==city.names[k]) %>% pull(delta)
        # phi_ts <- ee %>% filter(idx==city.names[k]) %>% pull(phi_ts)
        aq_t <- ee %>% filter(idx==city.names[k]) %>% pull(aq_t)
        bq_t <- ee %>% filter(idx==city.names[k]) %>% pull(bq_t)
        #yest = (1-betaa^(delta))^(phi_ts)
        yest = rkumar(107,aq_t,bq_t*H[k]) %>% t() %>% t()
        yest.full <- rbind(yest.full,yest)
      }

      estt.aux <- ee %>% mutate(yest = yest.full)

      estt <- rbind(estt,estt.aux)

      #
      # e = rkumar(nrow(ee),ee$aq_t,ee$bq_t*ee$H)
      # est_ <- ee %>% mutate(yk=e)
      # estt=rbind(estt,est_)
    }
    yest.aux[, b] <- estt$yest
  }
  #qest.aux[,l] <- apply(yest.aux, 1, quantile, prob = c(0.5))
  #qest.aux[,l] <- apply(yest.aux, 1, mean,trim=mtrim)
  #qest <- apply(qest.aux, 1, mean,trim=mtrim)
  #qest <- apply(qest.aux, 1, quantile, prob = c(0.5))
  # min.sample <- function(x){
  #   return(mean(sample(x,size = sqrt(B),prob=NULL,replace = T),trim=mtrim))
  # }
  #qest <- apply(yest.aux, 1,)
  #qest <- apply(yest.aux, 1, sample,1,replace=T,prob=NULL)
  qest <- apply(yest.aux, 1, median)
  data$id <- seq(1, nrow(data))
  row.deleted <- data %>%
    group_by(Group, City) %>%
    slice((1:(max(q, p)))) %>%
    data.frame() %>%
    select(id) %>%
    unlist()
  data[-row.deleted, "RH"] <- qest

  rht <- results$data%>%group_by(Group,City) %>% summarise(rht=RH) %>%pull(rht)
  dates <- seq.Date(as.Date("2000/01/01"),as.Date("2008/12/31"),by = "month",sep="/") %>% rep(14)
  fig1 = data %>% data.frame()%>%mutate(dates=dates,rht = rht) %>%
    select(dates,RH,Group,City,rht)%>%pivot_longer(cols = c(rht,RH)) %>%dplyr::group_split(Group) %>%
    purrr::map(
      ~ggplot(.,aes(colour=name)) + geom_line(aes(dates,value),size=.8)+
        #geom_point(aes(t,value))+
        ylab("Relative humidity")+ #scale_y_continuous(labels=scales::percent)+
        scale_x_date(date_labels = "%b-%y",breaks = c(seq(as.Date("2000/1/1"), as.Date("2008/12/31"), by = "years")),limits = c(as.Date("2000/01/01"),as.Date("2008/12/31")))+
        #scale_x_discrete(labels = c(month.abb))+
        facet_grid(~City,scale="free")+
        theme_bw()+#geom_hline(data = d,aes(yintercept = kp))+
        #geom_hline(yintercept = q)+
        theme(legend.position = "none",axis.text.x=element_text(angle=60,size=8,hjust = 1),
              axis.title.y = element_text(size=8)))%>%
    cowplot::plot_grid(plotlist = .,nrow=3)
  x11()
  fig1

  #saveRDS(data,"data_sintetic_analysis/data_sintetic_1.rds")

  #z1 = results$data %>% filter(City=="Manaus") %>% select(RH)
  # z2 = data %>%filter(City=="Manaus") %>% select(RH)
  # plot(1:120,z1$RH[1:120],type="o")
  # lines(1:120,z2$RH[1:120],type="o",col=2)
  N <- results$N
  erro <- results$erro
  formula <- results$formula
  w <- data$weights
  link_kl <- results$link.kl
  link_delta <- results$link.delta
  # tetastart <- c(ar,betas,lambdas,alfas)

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

  names <- c(
    results$ar %>% rownames(),
    paste0(colnames(cov_kl), "_kl"), paste0(colnames(cov_delta), "_delta"), paste0("alpha_", 1:length(alfa.est))
  )

  estimation <- data.frame(Estimation = unlist(c(ar_est, betas_est, lambdas_est, alfas_est)), Par = names)
  path_monte_carlo <- "estimations_predicts/estimation_mc.txt"
  write.table(estimation, file = path_monte_carlo, append = T, quote = T, col.names = F)

  log.like <- data.frame(log_like=mod$objective)
  path_log.like <- "estimations_predicts/log_like.txt"
  write.table(log.like, file = path_log.like, append = T, quote = T, col.names = F)

  # return(estimation)
}



