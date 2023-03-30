block.4 <- read.table("scripts_tests/model_time_v2/Block_4/EMV.txt")
block.5 <- read.table("scripts_tests/model_time_v2/Block_5/EMV.txt")
block.6 <- read.table("scripts_tests/model_time_v2/Block_6/EMV.txt")

estimates <- data.frame(id = block.4$V2,
                        b4.emv.cond.z=block.4$V4,
                        b5.emv.cond.z=block.5$V3,
                        b6.emv.marg=block.6$V3,
                        par.names=block.5$V4)
data("data_1")
head(data_1)
## Parse database------
month_names <- factor(month.abb, levels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))
database <- data_1 %>%
  filter(Station == "1") %>%
  mutate(
    semester= rep(c(1,2),each=6) %>% rep(21) %>% rep(14) %>% as.factor(),
    month_names = rep(month_names, 21) %>% rep(14),
    date = paste0(Year, "/", month_names),
    t = seq(1, 252) %>% rep(14),
    Group = rep(paste0("Group ", c(1, 2, 3, 3, 2, 1, 1, 3, 2, 3, 3, 3, 1, 1)), each = 252),
    cost = cos(2 * pi * as.numeric(Month) / 12),
    sent = sin(2 * pi * as.numeric(Month) / 12),
    lles = (  log(10^((7.5*TBS)/(237.3+TBS)))  )
  )

head(database, 13)
citys <- database$City %>% unique()
TT=143
p=10
data <- database %>% filter(City==citys[p]) %>% slice((252-TT):252)
data$precp[data$precp==0] <- 1
#p=10
formula <- RH ~  lles + sent+cost|semester+log(precp)
mf <- model.frame(Formula::Formula(formula), data = data)
y <- model.response(mf)
cov_a <- model.matrix(Formula::Formula(formula), data = data, rhs = 1)
cov_delta <- model.matrix(Formula::Formula(formula), data = data, rhs = 2)
# write.table(cov_a,"scripts_tests/Block_1/cov_a_manaus_m4s1.txt")
# write.table(cov_delta,"scripts_tests/Block_1/cov_delta_manaus_m4s1.txt")
ncx <- ncol(cov_a)
ncv <- ncol(cov_delta)


theta  <- (estimates %>% filter(id==citys[p]))[,4]


Beta = theta[1:ncx] # real
lambda = theta[(ncx+1):c(ncx+ncv)] # real
alpha = theta[-c(1:c(ncx+ncv))]
xbeta = cov_a %*% Beta
wlambda = cov_delta %*% lambda
delta = exp(wlambda)
b = 1 / delta
a = NULL

a[1] = exp(xbeta[1])
for (t in 2:(TT+1)) {
  a[t] = exp(xbeta[t])
}
mediana = ( 1-exp( - delta*(-log(0.5))^(1/alpha) ) )^(1/a)
plot.ts(y)
lines(mediana,col=2)
gt <- extraDistr::dkumar(y,a,b)
Gt <- extraDistr::pkumar(y,a,b,lower.tail = F)
LAMBDA <- Gt#-log(Gt)
r <- gt/Gt # derivada de -log(Gt)

LAMBDA.t.esp <- LAMBDA[143]^(1-alpha)/alpha

qx = .70
gt <- extraDistr::dkumar(y,a,b)
Gt <- extraDistr::pkumar(qx,a,b,lower.tail = F)
Gt[Gt==0] <- .Machine$double.eps
LAMBDA <- -log(Gt)
r <- gt/Gt # derivada de -log(Gt)

prob=NULL
for(t in 2:length(gt)){
  prob[t] <- (LAMBDA[t-1])^(1-alpha)*(LAMBDA[t]+LAMBDA[t-1])^(alpha-1)*
    exp( (LAMBDA[t-1])^alpha - (LAMBDA[t]+LAMBDA[t-1])^(alpha)    )
}
# prob = prob[-1]
# plot.ts(prob,type="o")
# lines(Gt,type="o",col=2)


data.full <- readRDS("data/data_3.rds")


prob.emp = data.full %>% filter(Year>=2009,City==citys[p]) %>%
  group_by(Year,Month) %>%
  summarise(prob.emp = length(RH[RH>qx])/length(RH)  ) %>%
  pull(prob.emp)
#lines(prob.emp,col=3)


dd = data.frame(id=(1:(length(y)-1)),y=y[-1],p1=prob,p2=Gt[-1],p3=prob.emp[-1])  %>%
  pivot_longer(cols = c(p1,p2,p3,y))
f1 =  dd %>% data.frame() %>%  ggplot() +
  geom_line(aes(id,value),stat = "identity")+
  geom_point(aes(id,value)) +
  facet_wrap(.~name,ncol = 1,scales = "free_y")

f1



# extemogram
par(mfrow = c(2, 2))
extremogram::extremogram1(y,.70,70,1)
acf(y)
extremogram::extremogram1(rnorm(144),.50,70,1)
