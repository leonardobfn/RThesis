rm(list = ls())
#packages--------------
require(ggplot2)
require(ggridges)
require(tidyr)
require(dplyr)
require(extraDistr)
compiler::enableJIT(3)
alphas.plot = c("0.35","0.65","0.95")
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


estimates_model_1$V3 = factor(estimates_model_1$V3,
                              levels = c(emv.1$par_names, "alpha"),
                              labels = c(expression(beta[0]),
                                         expression(beta[1]),
                                         expression(beta[2]),
                                         expression(gamma[0]),
                                         expression(gamma[1]),
                                         expression(alpha)))


estimates_model_1 = estimates_model_1 %>% rename(Parameter = V3)

for(i in 1:length(alphas.plot) ){

  fig_box_plot = estimates_model_1 %>% filter(alpha == alphas.plot[i]) %>%
    #ggplot(aes(Method,V2)) +
    ggplot()+
    #geom_boxplot(aes(Method,V2)) +
    #geom_ridgeline(scale = 2,alpha=0.75) +
    geom_density_ridges(aes(V2,Method),scale = 1)  +
    facet_wrap(Parameter~.,scales = "free",
               labeller = label_parsed) +
    xlab("Method")+
    ylab("Estimates")+
    theme_bw()
  path1 = paste0("Data_simulation/Figures/Model_1_",alphas.plot[i],".pdf")
  #path2 = paste0()
  pdf(path1)
  #postscript(path)
  print(fig_box_plot)
  dev.off()
}

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
estimates_model_2$V3 = factor(estimates_model_2$V3,
                              levels = c(emv.2$par_names, "alpha"),
                              labels = c(expression(beta[0]),
                                         expression(beta[1]),
                                         expression(beta[2]),
                                         expression(gamma[0]),
                                         expression(gamma[1]),
                                         expression(alpha)))


estimates_model_2 = estimates_model_2 %>% rename(Parameter = V3)

for(i in 1:length(alphas.plot) ){
  fig_box_plot = estimates_model_2 %>% filter(alpha == alphas.plot[i]) %>%
    #ggplot(aes(Method,V2)) +
    ggplot()+
    #geom_boxplot(aes(Method,V2)) +
    #geom_ridgeline(scale = 2,alpha=0.75) +
    geom_density_ridges(aes(V2,Method),scale = 1) +
    facet_wrap(Parameter~.,shrink = T,scales = "free",
               labeller = label_parsed) +
    xlab("Method")+
    ylab("Estimates")+
    theme_bw()
  path1 = paste0("Data_simulation/Figures/Model_2_",alphas.plot[i],".pdf")
  #path2 = paste0()
  pdf(path1)
  #postscript(path)
  print(fig_box_plot)
  dev.off()
}

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

estimates_model_3 %>% group_by(Method,alpha,V3) %>%
  summarise(iqr = IQR(V2,na.rm = T),
            q25 = quantile(V2,0.25,na.rm = T),
            q75 = quantile(V2,0.75,na.rm = T),
            l = q25-1.5*iqr,
            up =q75+1.5*iqr ,
            est = V2[V2<up & V2>l])


emv.3 = read.table("Data_simulation/real_par_model_3.txt")
estimates_model_3$V3 = factor(estimates_model_3$V3,
                              levels = c(emv.3$par_names, "alpha"),
                              labels = c(#expression(beta[0]),
                                expression(beta[1]),
                                expression(beta[2]),
                                expression(gamma[0]),
                                expression(gamma[1]),
                                expression(alpha)))

# removendo outilnes ---------


estimates_model_3 = estimates_model_3 %>% rename(Parameter = V3)

for(i in 1:length(alphas.plot) ){
  fig_box_plot = estimates_model_3 %>%
    group_by(Method,alpha,Parameter) %>%
    summarise(iqr = IQR(V2,na.rm = T),
              q25 = quantile(V2,0.25,na.rm = T),
              q75 = quantile(V2,0.75,na.rm = T),
              l = q25-1.5*iqr,
              up =q75+1.5*iqr ,
              V2 = V2[V2<up & V2>l]) %>%
    filter(alpha == alphas.plot[i]) %>%
    #ggplot(aes(Method,V2)) +
    ggplot()+
    #geom_boxplot(aes(Method,V2)) +
    #geom_ridgeline(aes(Method,V2),scale = 2,alpha=0.75) +
    geom_density_ridges(aes(V2,Method),scale = 1)  +
    facet_wrap(Parameter~.,shrink = T,scales = "free",
               labeller = label_parsed) +
    xlab("Method")+
    ylab("Estimates")+
    theme_bw()
  path1 = paste0("Data_simulation/Figures/Model_3_",alphas.plot[i],".pdf")
  #path2 = paste0()
  pdf(path1)
  #postscript(path)
  print(fig_box_plot)
  dev.off()
}

estimates_model_3 %>% filter(Method=="Method_1")
