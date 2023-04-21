require(ggplot2)
require(tidyverse)
require(dplyr)
alpha.values = c("alpha35", "alpha65", "alpha95")
models <- c("Model_1", "Model_2", "Model_3")
sc = .80
width=7*sc
height=7*sc
y = NULL
for (i in 1:length(models)) {
  for (j in 1:length(alpha.values)) {
    paths <-
      paste0("Data_simulation/",
             models[i],
             "/simulations/",
             alpha.values[j],
             "/data10.txt")
    y.aux = read.table(paths)
    alpha = stringr::str_sub(alpha.values[j], 6, 7)
    y.aux. = data.frame(y.aux, Model = models[i], alpha = alpha)
    y = rbind(y, y.aux.)
  }
}


#------------ séries------

fig_series = y %>%
  rename(t = V2, y = V1) %>% ggplot() +
  geom_line(aes(t, y)) +
  facet_grid(V3~Model  ,labeller = label_bquote(alpha == .(V3)  )) +
  theme_bw()
#fig_series

path1 = paste0("Data_simulation/Figures/series",".pdf")
#path1 = paste0("Data_simulation/Figures/series",".eps")
#path2 = paste0()
pdf(path1,width = width,height = height,title = "y")
#cairo_ps(path1)
print(fig_series)
dev.off()
#------------ ACF series------

ic_alpha= function(alpha, acf_res){
  return(qnorm((1 + (1 - alpha))/2)/sqrt(acf_res$n.used))
}

acff = y %>% group_by(Model,V3) %>%
    reframe(acff=acf(V1,plot = F)$acf[,,1],
            lag = acf(V1,plot=F)$lag,
            lim1=ic_alpha(0.05,acf(V1)))

fig_acf = acff %>%
  ggplot() +
  geom_segment(mapping = aes(x = lag,y=0, xend = lag, yend = acff))+
  geom_hline(aes(yintercept = 0)) +
  geom_hline(aes(yintercept = lim1), linetype = 2, color = 'blue') +
  geom_hline(aes(yintercept = -lim1), linetype = 2, color = 'blue')+
  facet_grid(V3~Model  ,labeller = label_bquote(alpha == .(V3)  )) +
  theme_bw()
#path1 = paste0("Data_simulation/Figures/series",".eps")
path1 = paste0("Data_simulation/Figures/acf",".pdf")
#path2 = paste0()
pdf(path1,width = width,height = height,title = "acf")
#cairo_ps(path1)
print(fig_acf)
dev.off()
