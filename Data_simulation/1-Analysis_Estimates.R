rm(list = ls())
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
compiler::enableJIT(3)

# Model 1------------
estimates_model_1 = NULL
methods <- c("Method_1", "Method_2", "Method_3")
#methods <- c("Method_1")
for (meth in methods) {
  path. = paste0("Data_simulation/Model_1/estimates/", meth, "/estimates/")
  LIST.FILES = list.files(path.)
  for (FILE in LIST.FILES) {
    path = paste0("Data_simulation/Model_1/estimates/",
                  meth,
                  "/estimates/",
                  FILE)
    estimates.aux = read.table(path)
    alpha. = paste0("0.", stringr::str_sub(FILE, start = 15, end = 16))
    estimates.aux = estimates.aux %>% mutate(alpha = alpha., Method = meth)
    estimates_model_1 = rbind(estimates_model_1, estimates.aux)
  }
}

estimates_model_1$alpha  = factor(estimates_model_1$alpha,
                                  levels = c("0.35", "0.50", "0.65", "0.80", "0.95"))

emv.1 = read.table("Data_simulation/real_par_model_1.txt")
estimates_model_1$V3  = factor(estimates_model_1$V3,
                               levels = c(emv.1$par_names, "alpha"))
estimates_model_1 = estimates_model_1 %>% rename(Parameter = V3)
seq.par <- seq(1,nrow(estimates_model_1))

alphas.par = c(0.35, 0.50, 0.65, 0.80, 0.95) %>% rep(each=1000) %>% rep(3)

seq.par.alpha <- seq(6,nrow(estimates_model_1),6)

seq.par[seq.par.alpha] <- alphas.par
seq.par[-seq.par.alpha] <- emv.1$real_values
conf = 0.95
estimates_model_1 %>% mutate(par_real = seq.par)
par.estimates_model_1 = estimates_model_1 %>%
  mutate(par_real = seq.par) %>%
  group_by(Method, alpha, Parameter) %>%
  summarise(
    par_real =  median(par_real),
    median = median(V2, na.rm = T),
    #sd_ = sd(V2, na.rm = T),
    MPE = mean((V2-par_real)/par_real,na.rm=T),
    relative_bias = ((median-par_real)/par_real)
    #inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    #inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
    #                         T)
    # par_real_in_interval = ifelse(par_real > inter.lower &
    #                                 par_real < inter.upper, 1, 0),
  ) %>%
  arrange(Method)

par.estimates_model_1 %>% filter(Parameter == "alpha")
write.csv(par.estimates_model_1,
                    "Data_simulation/estimates_values_model_1.csv")
read.csv("Data_simulation/estimates_values_model_1.csv")
# writexl::write_xlsx(par.estimates_model_1,
#           "Data_simulation/estimates_values_model_1.xlsx")
# Model 2------------


estimates_model_2 = NULL
methods <- c("Method_1", "Method_2", "Method_3")
#methods <- c("Method_2")
for (meth in methods) {
  path. = paste0("Data_simulation/model_2/estimates/", meth, "/estimates/")
  LIST.FILES = list.files(path.)
  for (FILE in LIST.FILES) {
    path = paste0("Data_simulation/model_2/estimates/",
                  meth,
                  "/estimates/",
                  FILE)
    estimates.aux = read.table(path)
    alpha. = paste0("0.", stringr::str_sub(FILE, start = 15, end = 16))
    estimates.aux = estimates.aux %>% mutate(alpha = alpha., Method = meth)
    estimates_model_2 = rbind(estimates_model_2, estimates.aux)
  }
}

estimates_model_2$alpha  = factor(estimates_model_2$alpha,
                                  levels = c("0.35", "0.50", "0.65", "0.80", "0.95"))

emv.2 = read.table("Data_simulation/real_par_model_2.txt")
estimates_model_2$V3  = factor(estimates_model_2$V3,
                               levels = c(emv.2$par_names, "alpha"))

estimates_model_2 = estimates_model_2 %>% rename(Parameter = V3)
seq.par <- seq(1,nrow(estimates_model_2))

alphas.par = c(0.35, 0.50, 0.65, 0.80, 0.95) %>% rep(each=1000) %>% rep(3)

seq.par.alpha <- seq(6,nrow(estimates_model_2),6)

seq.par[seq.par.alpha] <- alphas.par
seq.par[-seq.par.alpha] <- emv.2$real_values


conf = 0.95
par.estimates_model_2 = estimates_model_2 %>%
  mutate(par_real = seq.par) %>%
  group_by(Method, alpha, Parameter) %>%
  summarise(
    par_real =  median(par_real),
    median = median(V2, na.rm = T),
    #sd_ = sd(V2, na.rm = T),
    MPE = mean((V2-par_real)/par_real,na.rm=T),
    relative_bias = ((median-par_real)/par_real)
    #inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    #inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
    #                         T)
    # par_real_in_interval = ifelse(par_real > inter.lower &
    #                                 par_real < inter.upper, 1, 0),
  ) %>%
  arrange(Method)

par.estimates_model_2 %>% filter(Parameter == "alpha")

write.csv(par.estimates_model_1,
          "Data_simulation/estimates_values_model_2.csv")
# writexl::write_xlsx(par.estimates_model_2,
#                     "Data_simulation/estimates_values_model_2.xlsx")

# Model 3 ------------
estimates_model_3 = NULL
methods <- c("Method_1", "Method_2", "Method_3")
#methods <- c("Method_2")
for (meth in methods) {
  path. = paste0("Data_simulation/Model_3/estimates/", meth, "/estimates/")
  LIST.FILES = list.files(path.)
  for (FILE in LIST.FILES) {
    path = paste0("Data_simulation/Model_3/estimates/",
                  meth,
                  "/estimates/",
                  FILE)
    estimates.aux = read.table(path)
    alpha. = paste0("0.", stringr::str_sub(FILE, start = 15, end = 16))
    estimates.aux = estimates.aux %>% mutate(alpha = alpha., Method = meth)
    estimates_model_3 = rbind(estimates_model_3, estimates.aux)
  }
}

estimates_model_3$alpha  = factor(estimates_model_3$alpha,
                                  levels = c("0.35", "0.50", "0.65", "0.80", "0.95"))

emv.3 = read.table("Data_simulation/real_par_model_3.txt")
estimates_model_3$V3  = factor(estimates_model_3$V3,
                               levels = c(emv.3$par_names, "alpha"))
estimates_model_3 = estimates_model_3 %>% rename(Parameter = V3)
seq.par <- seq(1,nrow(estimates_model_3))

alphas.par = c(0.35, 0.50, 0.65, 0.80, 0.95) %>% rep(each=1000) %>% rep(3)

seq.par.alpha <- seq(5,nrow(estimates_model_3),5)

seq.par[seq.par.alpha] <- alphas.par
seq.par[-seq.par.alpha] <- emv.3$real_values


conf = 0.95
par.estimates_model_3 = estimates_model_3 %>%
  mutate(par_real = seq.par) %>%
  group_by(Method, alpha, Parameter) %>%
  summarise(
    par_real =  median(par_real),
    median = median(V2, na.rm = T),
    #sd_ = sd(V2, na.rm = T),
    MPE = mean((V2-par_real)/par_real,na.rm=T),
    relative_bias = ((median-par_real)/par_real)
    #inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    #inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
    #                         T)
    # par_real_in_interval = ifelse(par_real > inter.lower &
    #                                 par_real < inter.upper, 1, 0),
  ) %>%
  arrange(Method)

par.estimates_model_3 %>% filter(Parameter == "alpha")
write.csv(par.estimates_model_1,
          "Data_simulation/estimates_values_model_3.csv")
# writexl::write_xlsx(par.estimates_model_2,
#                     "Data_simulation/estimates_values_model_3.xlsx")
