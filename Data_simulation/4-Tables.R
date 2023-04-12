rm(list = ls())
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(xtable)
require(kableExtra)
compiler::enableJIT(3)
estimates_model_1 = read.csv("Data_simulation/estimates_values_model_1.csv")
emv.1 = read.table("Data_simulation/real_par_model_1.txt")
estimates_model_1$Parameter = factor(
  estimates_model_1$Parameter,
  levels = c(emv.1$par_names, "alpha"),
  labels = c(
    paste0("beta_0"),
    paste0("beta_1"),
    paste0("beta_2"),
    paste0("gamma_0"),
    paste0("gamma_1"),
    paste0("alpha")
  )
)


teste = estimates_model_1 %>% filter(alpha == "0.35") %>% arrange(Parameter)
rownames(teste) = paste0(teste$X, teste$Parameter)
teste0 = t(teste) %>% data.frame()
m1 = teste0[-c(1, 2, 3, 4), seq(1, ncol(teste0), 3)]
m2 = teste0[-c(1, 2, 3, 4), seq(2, ncol(teste0), 3)]
m3 = teste0[-c(1, 2, 3, 4), seq(3, ncol(teste0), 3)]
colnames(m1) = colnames(m2) = colnames(m3) = c(c(1:6) )
teste1 = rbind(m1, m2, m3)
betas.names = paste0("$\\beta_", 0:2, "$")
gamma.names = paste0("$\\gamma_", 0:1, "$")
alpha.names = paste0("$\\alpha$")
names.par = c("Methods", "Statistics", betas.names, gamma.names, alpha.names)
meth = rep(c("Method 1", "Method 2", "Method 3"), each = 5)
s = c("Real", "Mean", "Median", "MSE", "RB") %>% rep(3)
teste = teste1 %>% mutate(Method = meth, Statistics = s) %>%
  select(Method, Statistics, everything())
table = kbl(
  teste,
  align = ,
  format = "latex"  ,
  booktabs = T,
  col.names = names.par,
  escape = F,
  row.names = F,
  longtable = F,caption = "teste",label = "sadf"
  ) %>%
  kable_paper(full_width = F) %>%
  column_spec(1:2, bold = T) %>%
  collapse_rows(columns = 1:2, valign = "middle")

save_kable(table,"test2.tex")
