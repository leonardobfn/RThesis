rm(list = ls())
#packages--------------
require(tidyr)
require(dplyr)
require(extraDistr)
require(xtable)
require(kableExtra)
compiler::enableJIT(3)
al = c(0.35, 0.50, 0.65,0.80,0.95)
estimates_model_2 = read.csv("Data_simulation/estimates_values_model_2.csv")
emv.1 = read.table("Data_simulation/real_par_model_2.txt")
estimates_model_2$Parameter = factor(
  estimates_model_2$Parameter,
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

for (AL in al) {
  model_2 = estimates_model_2 %>% filter(alpha == AL) %>% arrange(Parameter)
  rownames(model_2) = paste0(model_2$X, model_2$Parameter)
  model_2 = t(model_2) %>% data.frame()
  m1 = model_2[-c(1, 2, 3, 4), seq(1, ncol(model_2), 3)]
  m2 = model_2[-c(1, 2, 3, 4), seq(2, ncol(model_2), 3)]
  m3 = model_2[-c(1, 2, 3, 4), seq(3, ncol(model_2), 3)]
  colnames(m1) = colnames(m2) = colnames(m3) = c(c(1:6))
  model_2 = rbind(m1, m2, m3)
  betas.names = paste0("$\\beta_", 0:2, "$")
  gamma.names = paste0("$\\gamma_", 0:1, "$")
  alpha.names = paste0("$\\alpha$")
  names.par = c("Methods", "Statistics", betas.names, gamma.names, alpha.names)
  meth = rep(c("Method 1", "Method 2", "Method 3"), each = 5)
  stats = c("Real", "Mean", "Median", "MSE", "RB") %>% rep(3)
  model_2 = model_2 %>% mutate(Method = meth, Statistics = stats) %>%
    select(Method, Statistics, everything())

  caption_model_2 = paste0("Estimates of parameters to Model 2 and ", "$\\alpha=", AL, "$")

  label_model_2 = paste0("Est_model_2_Alpha", AL)
  table_model_2 = kbl(
    model_2,
    align = ,
    format = "latex"  ,
    booktabs = T,
    col.names = names.par,
    escape = F,
    row.names = F,
    longtable = F,
    caption = caption_model_2,
    label = label_model_2
  ) %>%
    kable_paper(full_width = F) %>%
    column_spec(1:2, bold = T) %>%
    #row_spec(c(1:5, 11:15)-1, extra_latex_after = "\\rowcolor{gray}") %>%
    collapse_rows(columns = 1:2, valign = "middle")
  al.names = stringr::str_sub(AL, 3, 4)
  path_model_2 = paste0("Data_simulation/Tables/Table_model_2_Alpha",
                        al.names,
                        ".tex")
  save_kable(table_model_2, path_model_2)

}
