library(tidyverse)
library(grf)

set.seed(1)

org.data = read.csv("./acic_data.csv")
n = nrow(org.data)
X = org.data %>% select(-Z,-Y) %>% as.matrix

cat(" - Generating synthetic population...")

re.ind = sample(n,n)
# fit E[Y(0)|X] on 1/5 data
fit.data = org.data[re.ind[1:floor(n/5)],]
fit.X = X[re.ind[1:floor(n/5)],]
fit.model = regression_forest(fit.X[fit.data$Z==0,], fit.data$Y[fit.data$Z==0], num.threads=1)
ps.model = regression_forest(fit.X, fit.data$Z, num.threads=1)
qtl.model = quantile_forest(fit.X[fit.data$Z==1,], fit.data$Y[fit.data$Z==1], num.threads=1)

# construct \hat{m}_0 on the remaining data
data = org.data[re.ind[(floor(n/5)+1):n],]
X = data %>% select(-Z,-Y) %>% as.matrix
hat.m0 = predict(fit.model, X, num.threads=1)$predictions
gam = rnorm(76) * 0.105
tau.x = 0.228 + 0.05 * (X[,7] < 0.07) - 0.05 * (X[,8] <- 0.69) - 0.08 * (X[,3] %in% c(1,13,14)) + gam[X[,1]]
p.score = predict(ps.model, X, num.threads=1)$predictions


qtls = predict(qtl.model, quantile=c(0.25,0.75), newdata = X)
hat.rx = qtls[,2] - qtls[,1]

cat("Done.\n")

pop.data = list("X" = X, "mu0" = hat.m0, "tau.x" = tau.x, "p.score" = p.score, "r.x" = hat.rx)
save(pop.data, file="synthetic_population.RData")

