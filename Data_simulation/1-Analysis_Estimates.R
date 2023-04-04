

rm(list = ls())
estimates = NULL

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
    estimates = rbind(estimates, estimates.aux)
  }
}

estimates$alpha  = factor(estimates$alpha,
                          levels = c("0.35", "0.50", "0.65", "0.80", "0.95"))

emv = read.table("Data_simulation/real_par_model_1.txt")
estimates$V3  = factor(estimates$V3,
                          levels = c(emv$par_names,"alpha"))
conf = 0.95
par.estimates = estimates %>% group_by(Method, alpha, V3) %>%
  summarise(median = median(V2, na.rm = T),
            sd_ = sd(V2,na.rm = T),
            inter.lower = quantile(V2,c((1-conf)/2),na.rm=T),
            inter.upper = quantile(V2,c(conf+(1-conf)/2),na.rm=T)
            ) %>%
  arrange(Method)

alphas = c(0.35,0.50,0.65,0.80,0.95)
seq_alpha = seq(6,90,6)
par_real = NULL
par_real[seq_alpha] = alphas
par_real[-seq_alpha] = c(emv$real_values)

par.estimates.complet = par.estimates %>%
  data.frame %>%
  mutate(par_real = par_real) %>%
  group_by(Method, alpha, V3) %>%
  summarise(par_real=par_real,
            median = median,
            inter.lower = inter.lower,
            inter.upper = inter.upper,
            par_real_in_interval = ifelse(par_real>inter.lower & par_real<inter.upper,1,0),
            sd_ = sd_,
            bias = median-par_real,
            relative_bias=bias/par_real ) %>% data.frame()

subb = par.estimates.complet %>% filter(Method == "Method_2")

#write.csv(par.estimates.complet,"Data_simulation/estimates_values.csv")
