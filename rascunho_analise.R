rm(list=ls())
require(tidyr)
require(dplyr)

dados1 = read.table("Data_simulation/Model_1/estimates/Method_2/estimates/estimates_alpha35_n144.txt")
dados2 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessian_Beta_alpha95_n144.txt")
#dados3 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessianCovLambda_alpha95.txt")
#dados4 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessian_alpha95.txt")

nrow(dados1)
nrow(dados2)
#nrow(dados3)
#nrow(dados3)
#nrow(dados4)
conf = .05

id = dados1 %>% slice(1:6000) %>% pull(V1) %>%sort()
length(table(id))
id.f = seq(1,1200) %in% id
which(id.f==F)
length(which(id.f==F) )
length(id)



est1 = dados1 %>% group_by(V2,V3,V6,V7,V8) %>% summarise(real_par = median(V4,na.rm=T)
                                                         ,m = median(V5,na.rm=T),
                                           l = quantile(V5,conf/2,na.rm=T),
                                           up = quantile(V5,1-(conf/2),na.rm=T))

est1
(dados1 %>% slice(which(is.na(V5)==T)) %>% nrow())/6
tt = dados1 %>% slice(which(V2>20)) %>% pull(V1)
est1 = dados1 %>% slice(tt) %>% group_by(V3) %>% summarise(m = median(V2,na.rm=T),
                                             l = quantile(V2,conf/2,na.rm=T),
                                             up = quantile(V2,1-(conf/2),na.rm=T))

est1
read.table("Data_simulation/real_par_model_1.txt")

