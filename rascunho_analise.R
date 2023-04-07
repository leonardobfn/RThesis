rm(list=ls())
require(tidyr)
require(dplyr)
# dados1 = read.table("Data_simulation/Model_1/estimates/Method_2/estimates/teste/estimatesalpha35.txt")
# dados2 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/teste/hessianCov_alpha35.txt")
# dados3 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/teste/hessian_alpha35.txt")

dados1 = read.table("Data_simulation/Model_1/estimates/Method_1/estimates/estimatesalpha35.txt")
dados2 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessianCov_alpha35.txt")
dados3 = read.table("Data_simulation/Model_1/estimates/Method_1/hessian/hessian_alpha35.txt")
#65,95,80

# dados1 = read.table("Data_simulation/Model_1/estimates/Method_2/estimates/estimatesalpha95.txt")
# dados2 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/hessianCov_alpha95.txt")

# dados2 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/hessianCovbeta_alpha95.txt")
# dados3 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/hessianCovLambda_alpha95.txt")
#dados4 = read.table("Data_simulation/Model_1/estimates/Method_2/hessian/hessian_alpha95.txt")

nrow(dados1)/6
nrow(dados2)
nrow(dados3)
#nrow(dados3)
#nrow(dados4)
conf = .05

id = dados1$V1 %>% unique()%>%sort()
id.f = seq(1,1000) %in% id
which(id.f==F)
length(which(id.f==F) )
length(id)

dados1 %>% filter(V1 %in% c(709 ,733))


est = dados1 %>% group_by(V3) %>% summarise(m = median(V2,na.rm=T),
                                           l = quantile(V2,conf/2,na.rm=T),
                                           up = quantile(V2,1-(conf/2),na.rm=T))

est
dados1 %>% slice(which(is.na(V2)==T)) %>% pull(V1) %>% unique() %>% length()

