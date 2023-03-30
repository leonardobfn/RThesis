rm(list = ls())

# Packages -------
pack <- c("tidyverse", "extraDistr", "devtools", "Formula", "tictoc", "betareg", "cowplot","rgdal")

package.check <- lapply(
  pack,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

devtools::load_all() # loading my functions

# Database------

## Loading database------

data("data_1")
head(data_1)
## Parse database------
month_names <- factor(month.abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
database <- data_1 %>%
  filter(Station == "1") %>%
  mutate(
    semester= rep(c(1,2),each=6) %>% rep(21) %>% rep(14) %>% as.factor(),
    month_names = rep(month_names, 21) %>% rep(14),
    date = paste0(Year, "/", month_names),
    t = seq(1, 252) %>% rep(14),
    Group = rep(paste0("Group ", c(1, 2, 3, 3, 2, 1, 1, 3, 2, 3, 3, 3, 1, 1)), each = 252),
    cost = cos(2 * pi * as.numeric(Month) / 12),
    sent = sin(2 * pi * as.numeric(Month) / 12)
  )

head(database, 13)

y.mat <- matrix(database$RH,nrow = 252,ncol = 14)
colnames(y.mat) <- database$City %>% unique()
cor(y.mat)
## Calculating the weights ---------

u1 <- "Nova Olinda do Norte"
u2 <- "Maraã"
u3 <- "Itamarati"
u <- c(u1, u2, u3)

w <- list()
for (j in 1:3) {
  g <- paste0("Group ", j)
  dados_aux <- database %>% dplyr::filter(Group == g)
  u_lat_log <- dplyr::filter(data_1, City == u[j])[c("City", "Latitude", "Longitude")]
  lat_log <- unique(dados_aux[c("City", "Latitude", "Longitude")])
  aux <- rbind(lat_log, u_lat_log)
  w[j] <- weight(mun = aux$City, u_m = u[j], lat = as.numeric(aux$Latitude), long = as.numeric(aux$Longitude))
}
names(w) <- c(paste0("Group ", 1:length(u))) # u
w
sapply(w, sum)
# names(w)=u
w.data.frame <- tibble::enframe(w) %>%
  tidyr::unnest(cols = value) %>%
  select(value) %>%
  unlist() %>%
  rep(each = 252)
database <- database %>%
  arrange(Group) %>%
  mutate(weights = w.data.frame) %>%
  arrange(Group)
head(database)
# Estimation ---------
TT <- 108 #
formula <- RH ~  Altitude + sent + cost  |  semester
data <- database %>%
  group_by(City) %>%
  slice(1:TT) #%>% filter(City == "Barcelos")
# data <- estimations$data
grupos <- data$Group
w <- data$weights
quantil <- .5
N <- 500
erro <- 10^(-4)
link.kl <- "loglog"
link.delta <- "log"
p <- 1 # AR(p)
q <- 0 # MA(q)
idx <- data$City
# Arrumandos do dados -parte 2 -------

mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_kl <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)

ncx <- ncol(cov_kl)
ncv <- ncol(cov_delta)
grupos <- as.factor(grupos)
gr <- levels(grupos)
DADOS <- data.frame(idx, y, cov_kl, cov_delta, grupos, w) %>% arrange(grupos)
colnames(DADOS) <- c(
  "idx", "y",
  colnames(cov_kl),
  paste0(colnames(cov_delta), "_delta"),
  "grupos", "weights"
)

Lm <- DADOS %>%
  group_by(grupos) %>%
  summarise(Lm = length(unique(idx))) %>%
  select(Lm) %>%
  unlist() # number of station per group
# vmm <- logrr <- vector(l = length(gr))
 Esp <- matrix(0, length(gr), 2)


# link function----

if (link.kl != "loglog") {
  link_kl <- make.link(link.kl)
}
if (link.kl == "loglog") {
  link_kl <- structure(list(
    linkfun = function(mu) -log(-log(mu)),
    linkinv = function(eta) {
      pmax(pmin(
        exp(-exp(-eta)),
        1 - .Machine$double.eps
      ), .Machine$double.eps)
    },
    mu.eta = function(eta) {
      eta <- pmin(eta, 700)
      pmax(exp(-eta - exp(-eta)), .Machine$double.eps)
    }, dmu.deta = function(eta) {
      pmax(exp(-exp(-eta) -
                 eta) * expm1(-eta), .Machine$double.eps)
    }, valideta = function(eta) TRUE,
    name = "loglog"
  ), class = "link-glm")
}

if (link.delta == "identity" || link.delta == "log" || link.delta == "sqrt") {
  link_delta <- make.link(link.delta)
}

# Chute inicial ------------

if (p != 0 & q != 0) {
  ar.ma.start <- coef(arima(y, c(p, 0, q), include.mean = T))[1:(p + q)]
  # ar_start <- ar.ma.start[1:p]
  # ma_start <- ar.ma.start[-(1:p)]
}

if (p != 0 & q == 0) {
  ar.ma.start <- coef(arima(y, c(p, 0, q), include.mean = T))[1:(p)]
  # ar_start <- ar.ma.start[1:p]
}

if (p == 0 & q != 0) {
  ar.ma.start <- coef(arima(y, c(p, 0, q), include.mean = T))[1:(q)]
  # ar_start <- ar.ma.start[1:p]
}

beta_start <- optim(c(rep(0, ncx)), qbeta_fit, quantil = quantil, control = list(maxit = 600), link = link_kl, cov_kl = cov_kl, y = y)$par
delta_start <- -betareg::betareg(
  formula = formula, data = data,
  link.phi = "log", link = "loglog",
)$coefficients$precision

alpha_start <- c(rep(0.90,length(Lm)))

tetastart <- c(ar.ma.start, beta_start, delta_start, alpha_start)
tetastart_p <- tetastart
names(tetastart) <- c(
  names(ar.ma.start),
  colnames(cov_kl),
  colnames(cov_delta),
  paste0("alpha_", 1:length(alpha_start))
)


lower <- c(rep(-20,length(ar.ma.start)),
           rep(-20,length(beta_start)),
           rep(-20,length(delta_start)),
           rep(0.1,length(alpha_start)))

upper <- c(rep(20,length(ar.ma.start)),
           rep(20,length(beta_start)),
           rep(20,length(delta_start)),
           rep(0.95,length(alpha_start)))


EMVteta <- EMVteta.aux <- EMVteta.aux.h <- list()

## função para maximizar-------

GEV_kuma_fit_r1 <- function(teta)
{
  #teta <- tetastart
  ar.ma <- teta[1:(p + q)]
  BETA <- teta[-c(1:(p + q))][1:ncx]
  LAMBDA <- teta[-c(1:(p + q))][(ncx + 1):(ncx + ncv)]
  #alpha <- teta[-c(1:(p + q))][-c(1:(ncx + ncv))]
  vm.1 <- logr.1 <- vector(l = length(alpha_up))

  for (j in 1:length(Lm)) {
    data_group <- DADOS %>% dplyr::filter(grupos == gr[j])
    T_star <- (length(data_group$y) / Lm[j]) - max(p,q)

    cov_delta_grupo <- as.matrix((data_group[c(paste0(colnames(cov_delta), "_delta"))]))
    eta_delta <- cov_delta_grupo %*% LAMBDA
    data_group <- data_group %>% mutate(delta = link_delta$linkinv(eta_delta))

    cov_kl_grupo <- as.matrix((data_group[c(colnames(cov_kl))]))

    eta_1 <- cov_kl_grupo %*% BETA
    eta_2 <- link_kl$linkfun(data_group$y) - eta_1

    # eta_1.aux <- matrix(eta_1,T_star+1,Lm[j])
    # eta_2.aux <- matrix(eta_2,T_star+1,Lm[j])

    # i <- max(p, q)
    # ar <- seq(1:p)
    # ma <- seq(1:q)
    # res <- eta_klt <- matrix(0,nrow(eta_2.aux),Lm[j])
    # y.g = matrix(data_group$y,nrow(eta_2.aux),Lm[j])
    #
    # if (p != 0 & q != 0) { # Gkarma(p,q) model
    #   phi_ar <- ar.ma[1:p]
    #   phi_ma <- ar.ma[-(1:p)]
    #   for (t in (1 + i):nrow(eta_2.aux)) {
    #     for(l in 1:Lm[j]){
    #       eta_klt[t,l] <- eta_1.aux[t,l] + phi_ar %*% eta_2.aux[t - ar,l] +  phi_ma %*% res[t - ma,l]
    #       res[t,l] <- link_kl$linkfun(y.g[t,l]) - eta_klt[t,l]
    #     }
    #   }
    # }
    #
    #
    # if (p != 0 & q == 0) { # Gkarma(p,q) model
    #   phi_ar <- ar.ma[1:p]
    #   for (t in (1 + i):nrow(eta_2.aux)) {
    #     for(l in 1:Lm[j]){
    #       eta_klt[t,l] <- eta_1.aux[t,l] + phi_ar %*% eta_2.aux[t - ar,l]
    #     }
    #   }
    # }
    #
    # kl_t <- eta_klt[eta_klt != 0]
    # kl_t <- link_kl$linkinv(kl_t)
    #

    kl_t <- data_group %>%
      mutate(eta_1 = eta_1, eta_2 = eta_2) %>%
      group_by(idx) %>%
      summarise(kl_t = klt(eta_1, eta_2, BETA = BETA, ar.ma, q = q, p = p, link_kl = link_kl, y.d = data_group$y)) %>%
      data.frame() %>%
      pull(kl_t)

    data_group <- data_group %>%
      group_by(idx) %>%
      slice(-(1:max(p, q)))

    dados_g <- data_group
    # <- data.frame(data_group, kl_t, delta=data_group$delta, omega=data_group$weights)

    bq_t <- 1 / dados_g$delta[,1] # bq_t[bq_t==0] <- .Machine$double.eps

    phi_ts <- log(kl_t) / log(1 - (1 - quantil)^dados_g$delta[,1])

    # phi_ts[phi_ts<=0]=.Machine$double.eps

    aq_t <- 1 / phi_ts
    aq_t[aq_t <= 0] <- .Machine$double.eps

    p_kumar <- extraDistr::pkumar(dados_g$y, a = aq_t, b = bq_t, lower.tail = F)
    p_kumar[p_kumar == 0] <- 1
    d_kumar <- extraDistr::dkumar(dados_g$y, a = aq_t, b = bq_t)


    r <- d_kumar / p_kumar
    r[r == Inf | r == "NaN" | r == 0] <- 1
    vmteta <- -dados_g$weights * log(p_kumar)
    vm.1[j] <- sum(vmteta)
    logr.1[j] <- sum(log(r))

  }

  wlog <- DADOS %>%
    select(grupos, idx, weights) %>%
    unique.data.frame() %>%
    group_by(grupos) %>%
    summarise(wlog = sum(log(weights))) %>%
    select(wlog) %>%
    unlist()

  Q1 <- sum(T_star * wlog - vm.1 * Esp[, 1] + logr.1)
  Q2 <- sum(-log(pi)-(alpha_up+1)*Esp[,2]+log(sin(alpha_up*pi))+T_star * Lm * Esp[, 2])

  L <- Q1 + Q2

  return(L)
}



repeat{

  ar.ma <- tetastart[1:(p + q)]
  BETA <- tetastart[-c(1:(p + q))][1:ncx]
  LAMBDA <- tetastart[-c(1:(p + q))][(ncx + 1):(ncx + ncv)]
  alpha <- tetastart[-c(1:(p + q))][-c(1:(ncx + ncv))]


  for (j in 1:length(Lm)) {

    data_group <- DADOS %>% dplyr::filter(grupos == gr[j])
    T_star <- (length(data_group$y) / Lm[j]) - max(p,q)

    cov_delta_grupo <- as.matrix((data_group[c(paste0(colnames(cov_delta), "_delta"))]))
    eta_delta <- cov_delta_grupo %*% LAMBDA
    data_group <- data_group %>% mutate(delta = link_delta$linkinv(eta_delta))

    cov_kl_grupo <- as.matrix((data_group[c(colnames(cov_kl))]))

    eta_1 <- cov_kl_grupo %*% BETA
    eta_2 <- link_kl$linkfun(data_group$y) - eta_1

    # eta_1.aux <- matrix(eta_1,T_star+1,Lm[j])
    # eta_2.aux <- matrix(eta_2,T_star+1,Lm[j])
    #
    # i <- max(p, q)
    # ar <- seq(1:p)
    # ma <- seq(1:q)
    # res <- eta_klt <- matrix(0,nrow(eta_2.aux),Lm[j])
    # y.g = matrix(data_group$y,nrow(eta_2.aux),Lm[j])
    #
    # if (p != 0 & q != 0) { # Gkarma(p,q) model
    #   phi_ar <- ar.ma[1:p]
    #   phi_ma <- ar.ma[-(1:p)]
    #   for (t in (1 + i):nrow(eta_2.aux)) {
    #     for(l in 1:Lm[j]){
    #       eta_klt[t,l] <- eta_1.aux[t,l] + phi_ar %*% eta_2.aux[t - ar,l] +  phi_ma %*% res[t - ma,l]
    #       res[t,l] <- link_kl$linkfun(y.g[t,l]) - eta_klt[t,l]
    #     }
    #   }
    # }
    #
    #
    # if (p != 0 & q == 0) { # Gkarma(p,q) model
    #   phi_ar <- ar.ma[1:p]
    #   for (t in (1 + i):nrow(eta_2.aux)) {
    #     for(l in 1:Lm[j]){
    #       eta_klt[t,l] <- eta_1.aux[t,l] + phi_ar %*% eta_2.aux[t - ar,l]
    #     }
    #   }
    # }

    # kl_t <- eta_klt[eta_klt != 0]
    # kl_t <- link_kl$linkinv(kl_t)


    kl_t <- data_group %>%
      mutate(eta_1 = eta_1, eta_2 = eta_2) %>%
      group_by(idx) %>%
      summarise(kl_t = klt(eta_1, eta_2, BETA = BETA, ar.ma, q = q, p = p, link_kl = link_kl, y.d = data_group$y)) %>%
      data.frame() %>%
      pull(kl_t)

    data_group <- data_group %>%
      group_by(idx) %>%
      slice(-(1:max(p, q)))

    dados_g <- data_group
    # <- data.frame(data_group, kl_t, delta=data_group$delta, omega=data_group$weights)

    bq_t <- 1 / dados_g$delta[,1] # bq_t[bq_t==0] <- .Machine$double.eps

    phi_ts <- log(kl_t) / log(1 - (1 - quantil)^dados_g$delta[,1])

    # phi_ts[phi_ts<=0]=.Machine$double.eps

    aq_t <- 1 / phi_ts
    aq_t[aq_t <= 0] <- .Machine$double.eps

    p_kumar <- extraDistr::pkumar(dados_g$y, a = aq_t, b = bq_t, lower.tail = F)
    p_kumar[p_kumar == 0] <- 1
    d_kumar <- extraDistr::dkumar(dados_g$y, a = aq_t, b = bq_t)


    r <- d_kumar / p_kumar
    r[r == Inf | r == "NaN" | r == 0] <- 1
    vmteta <- -dados_g$weights * log(p_kumar)
    vmm <- sum(vmteta)
    #logrr[j] <- sum(log(r))

    Esp[j,1] <- (Lm[j]*T_star - alpha[j])/vmm #E(zm|y)
    Esp[j,2] <- digamma(Lm[j]*T_star-alpha[j])-log(vmm) #E(log.zm|y)
  }

  ## STEP M--------------
  x = pi/Esp[,2]
  angle <- 1/x
  alpha_up <- ifelse(x>0,atan(x)*pi^(-1),(pi+atan(x))*pi^(-1)) # verificar

  EMV <- try(optim(par = tetastart[-(length(tetastart):(length(tetastart)-2))],
                   fn = GEV_kuma_fit_r1 ,
                   method = "BFGS",
                   control = list(fnscale=-1)),
             TRUE)


  EMVteta$par <- c(EMV$par,alpha_up)

  names(EMVteta$par) <- c(
    names(ar.ma.start),
    colnames(cov_kl),
    colnames(cov_delta),
    paste0("alpha_", 1:length(alpha_start))
  )

  crit <- sum(((EMVteta$par - tetastart) / tetastart)^2)
  # if(crit<erro){
  #   cont = cont + 1
  #   EMVteta.aux <- rbind(EMVteta.aux,EMVteta$par)
  # }

  cat("Estimates:", sep = "\n")
  print(EMVteta$par)
  cat("Crit:", sep = "\n")
  print(crit)
  if (crit<erro) {
    break
  } else {
    tetastart <- EMVteta$par
  }


}

# Análise
ar.ma <- EMVteta$par[1:(p + q)]
BETA <- EMVteta$par[-c(1:(p + q))][1:ncx]
LAMBDA <- EMVteta$par[-c(1:(p + q))][(ncx + 1):(ncx + ncv)]
alpha <- EMVteta$par[-c(1:(p + q))][-c(1:(ncx + ncv))]
vm.1 <- logr.1 <- vector(l = length(alpha))
for (j in 1:length(Lm)) {
  data_group <- DADOS %>% dplyr::filter(grupos == gr[j])
  T_star <- (length(data_group$y) / Lm[j]) - max(p,q)

  cov_delta_grupo <- as.matrix((data_group[c(paste0(colnames(cov_delta), "_delta"))]))
  eta_delta <- cov_delta_grupo %*% LAMBDA
  data_group <- data_group %>% mutate(delta = link_delta$linkinv(eta_delta))

  cov_kl_grupo <- as.matrix((data_group[c(colnames(cov_kl))]))
  eta_1 <- cov_kl_grupo %*% BETA
  eta_2 <- link_kl$linkfun(data_group$y) - eta_1

  kl_t <- data_group %>%
    mutate(eta_1 = eta_1, eta_2 = eta_2) %>%
    group_by(idx) %>%
    summarise(kl_t = klt(eta_1, eta_2, BETA = BETA, ar.ma, q = q, p = p, link_kl = link_kl, y.d = data_group$y)) %>%
    data.frame() %>%
    select(kl_t)

  data_group <- data_group %>%
    group_by(idx) %>%
    slice(-(1:max(p, q)))

  dados_g <- data_group
  # <- data.frame(data_group, kl_t, delta=data_group$delta, omega=data_group$weights)

  bq_t <- 1 / dados_g$delta[,1] # bq_t[bq_t==0] <- .Machine$double.eps

  phi_ts <- log(kl_t[,1]) / log(1 - (1 - quantil)^dados_g$delta[,1])

  # phi_ts[phi_ts<=0]=.Machine$double.eps

  aq_t <- 1 / phi_ts
  aq_t[aq_t <= 0] <- .Machine$double.eps

  p_kumar <- extraDistr::pkumar(dados_g$y, a = aq_t, b = bq_t, lower.tail = F)
  p_kumar[p_kumar == 0] <- 1
  d_kumar <- extraDistr::dkumar(dados_g$y, a = aq_t, b = bq_t)


  r <- d_kumar / p_kumar
  r[r == Inf | r == "NaN" | r == 0] <- 1
  vmteta <- -dados_g$weights * log(p_kumar)
  vm.1[j] <- sum(vmteta)
  logr.1[j] <- sum(log(r))
}



