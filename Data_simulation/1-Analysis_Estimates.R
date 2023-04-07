rm(list = ls())

#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
compiler::enableJIT(3)

# Model 1------------
estimates_model_1 = NULL
methods <- c("Method_1", "Method_2", "Method_3")
#methods <- c("Method_2")
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

emv = read.table("Data_simulation/real_par_model_1.txt")
estimates_model_1$V3  = factor(estimates_model_1$V3,
                               levels = c(emv$par_names, "alpha"))
conf = 0.95
par.estimates_model_1 = estimates_model_1 %>% group_by(Method, alpha, V3) %>%
  summarise(
    median = median(V2, na.rm = T),
    sd_ = sd(V2, na.rm = T),
    inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
                             T)
  ) %>%
  arrange(Method)

alphas = c(0.35, 0.50, 0.65, 0.80, 0.95)
seq_alpha = seq(6, 90, 6)
par_real = NULL
par_real[seq_alpha] = alphas
par_real[-seq_alpha] = c(emv$real_values)

par.estimates.complet_model_1 = par.estimates_model_1 %>%
  data.frame %>%
  mutate(par_real = par_real) %>%
  group_by(Method, alpha, V3) %>%
  summarise(
    par_real = par_real,
    median = median,
    inter.lower = inter.lower,
    inter.upper = inter.upper,
    par_real_in_interval = ifelse(par_real > inter.lower &
                                    par_real < inter.upper, 1, 0),
    sd_ = sd_,
    bias = median - par_real,
    relative_bias = bias / par_real
  ) %>% data.frame()
par.estimates.complet_model_1 %>% filter(V3 == "alpha")
writexl::write_xlsx(par.estimates.complet_model_1,
          "Data_simulation/estimates_values.xlsx")
# Model 2------------
estimates_model_2 = NULL
methods <- c("Method_1", "Method_2", "Method_3")
#methods <- c("Method_2")
for (meth in methods) {
  path. = paste0("Data_simulation/Model_2/estimates/", meth, "/estimates/")
  LIST.FILES = list.files(path.)
  for (FILE in LIST.FILES) {
    path = paste0("Data_simulation/Model_2/estimates/",
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

emv = read.table("Data_simulation/real_par_model_2.txt")
estimates_model_2$V3  = factor(estimates_model_2$V3,
                               levels = c(emv$par_names, "alpha"))
conf = 0.95
par.estimates_model_2 = estimates_model_2 %>% group_by(Method, alpha, V3) %>%
  summarise(
    median = median(V2, na.rm = T),
    sd_ = sd(V2, na.rm = T),
    inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
                             T)
  ) %>%
  arrange(Method)

alphas = c(0.35, 0.50, 0.65, 0.80, 0.95)
seq_alpha = seq(6, 90, 6)
par_real = NULL
par_real[seq_alpha] = alphas
par_real[-seq_alpha] = c(emv$real_values)

par.estimates.complet_model_2 = par.estimates_model_2 %>%
  data.frame %>%
  mutate(par_real = par_real) %>%
  group_by(Method, alpha, V3) %>%
  summarise(
    par_real = par_real,
    median = median,
    inter.lower = inter.lower,
    inter.upper = inter.upper,
    par_real_in_interval = ifelse(par_real > inter.lower &
                                    par_real < inter.upper, 1, 0),
    sd_ = sd_,
    bias = median - par_real,
    relative_bias = bias / par_real
  ) %>% data.frame()



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

emv = read.table("Data_simulation/real_par_model_3.txt")
estimates_model_3$V3  = factor(estimates_model_3$V3,
                               levels = c(emv$par_names, "alpha"))
conf = 0.95
par.estimates_model_3 = estimates_model_3 %>% group_by(Method, alpha, V3) %>%
  summarise(
    median = median(V2, na.rm = T),
    sd_ = sd(V2, na.rm = T),
    inter.lower = quantile(V2, c((1 - conf) / 2), na.rm = T),
    inter.upper = quantile(V2, c(conf + (1 - conf) / 2), na.rm =
                             T)
  ) %>%
  arrange(Method)

alphas = c(0.35, 0.50, 0.65, 0.80, 0.95)
seq_alpha = seq(6, 90, 6)
par_real = NULL
par_real[seq_alpha] = alphas
par_real[-seq_alpha] = c(emv$real_values)

par.estimates.complet_model_3 = par.estimates_model_3 %>%
  data.frame %>%
  mutate(par_real = par_real) %>%
  group_by(Method, alpha, V3) %>%
  summarise(
    par_real = par_real,
    median = median,
    inter.lower = inter.lower,
    inter.upper = inter.upper,
    par_real_in_interval = ifelse(par_real > inter.lower &
                                    par_real < inter.upper, 1, 0),
    sd_ = sd_,
    bias = median - par_real,
    relative_bias = bias / par_real
  ) %>% data.frame()
