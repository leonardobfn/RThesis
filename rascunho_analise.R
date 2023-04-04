require(tidyr)
require(dplyr)
dados = read.table("Data_simulation/Model_1/estimates/Method_1/estimates/estimatesalpha35")
head(dados,10)
conf = .05
est = dados %>% group_by(V3) %>% summarise(m = median(V2,na.rm=T),
                                     l = quantile(V2,conf/2,na.rm=T),
                                     up = quantile(V2,1-(conf/2),na.rm=T))
emv = read.table("Data_simulation/real_par_model_1.txt")



id = dados$V1 %>% unique()%>%sort()
id.f = seq(1,1000) %in% id
length(which(id.f==F) )
length(id)

est = dados %>% group_by(V3) %>% summarise(m = median(V2,na.rm=F),
                                           l = quantile(V2,conf/2,na.rm=F),
                                           up = quantile(V2,1-(conf/2),na.rm=F))

dados %>% slice(which(is.na(V2)==T)) %>% pull(V1) %>% unique() %>% length()
