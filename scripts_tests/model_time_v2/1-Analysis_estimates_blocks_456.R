

block.3 <- read.table("scripts_tests/model_time_v2/Block_3/EMV.txt")
block.4 <- read.table("scripts_tests/model_time_v2/Block_4/EMV.txt")
block.5 <- read.table("scripts_tests/model_time_v2/Block_5/EMV.txt")
block.6 <- read.table("scripts_tests/model_time_v2/Block_6/EMV.txt")

# estimates <- data.frame(id = block.3$V2,
#                         b3.emv.marg=block.3$V3,
#                         b3.emv.cond.z=block.3$V4,
#                         b4.emv.gbase=block.4$V3,
#                         b4.emv.cond.z=block.4$V4,
#                         b5.emv.cond.z=block.5$V3,
#                         par.names=block.5$V4)

estimates <- data.frame(id = block.4$V2,
                        b4.emv.cond.z=block.4$V4,
                        b5.emv.cond.z=block.5$V3,
                        b6.emv.marg=block.6$V3,
                        par.names=block.5$V4)

estimates.mun <- estimates %>%  filter(id=="Manaus") %>%
  pivot_longer(cols = colnames(estimates)[-c(1,5)])

estimates.mun %>%
  ggplot() +
  geom_text(aes(x=name,y=value ,label=round(value,3))) +
  facet_wrap(.~par.names,scales = "free")

