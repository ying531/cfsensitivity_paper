#!/usr/bin/env Rscript
########################################
## input configurations
########################################
args <- commandArgs(trailingOnly = TRUE)
alpha_ind <- as.integer(args[1])
Gamma_ind <- as.integer(args[2])
seed <- as.integer(args[3])

alphas = seq(0.1,0.9,by=0.1)
gammas = c(1.5,2,2.5,3,5)
# coverage target 1-alpha
alpha = alphas[alpha_ind] 
# confounding level Gamma
Gamma = gammas[Gamma_ind] 

########################################
## load libraries
########################################
suppressPackageStartupMessages(library(grf))
options(warn=-1)

########################################
## load util functions
########################################
source("../utils/util_synthetic.R")
source("../utils/util_ate.R")
load("synthetic_population.RData")
cat(paste(" - Running counterfactual prediction on synthetic data with PAC algorithm, alpha", alpha, ", Gamma",Gamma,
          ", seed", seed, "\n"), sep = '')

########################################
## Output direcroty
########################################
out_dir <- "../results/realdata/"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}


set.seed(seed)
delta = 0.05
########################################
## generate synthetic data
########################################
n = 10000
data = syn.gen(n, pop.data, Gamma, 1)
# data splitting
re.ind = sample(n,n)
train.ind = re.ind[1:floor(3*n/10)]
calib.ind = re.ind[floor(3*n/10+1):floor(9*n/10)]
test.ind = re.ind[floor(9*n/10+1):n]

########################################
## fit on the training fold
########################################
train.X = (data$X[train.ind,])[data$T[train.ind]==1,]
train.Y = (data$Y[train.ind])[data$T[train.ind]==1]
# train the nonconformity score function
train.score = conform.score(train.X, train.Y, "cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model
# estimate hat{e}(x)
hat.p = mean(data$T)
e.model = regression_forest(data$X[train.ind,], data$T[train.ind], num.threads = 1)

########################################
## calibration 
########################################
calib.X = (data$X[calib.ind,])[data$T[calib.ind]==1,] 
calib.Y = (data$Y[calib.ind])[data$T[calib.ind]==1]
calib.T = data$T[calib.ind]
calib.ex = predict(e.model, newdata=calib.X)$predictions 
n_calib = length(calib.Y)
# lower and upper bounds of weight function
calib.lx = hat.p * (1 + (1-calib.ex) / (calib.ex*Gamma))
calib.ux = hat.p * (1 + Gamma * (1-calib.ex) / (calib.ex))
calib.nc.w = hat.p / calib.ex
# non-conformity score on calibration data
calib.score = conform.score(calib.X, calib.Y, "cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.all = data.frame("score" = calib.score, "lx" = calib.lx, "ux" = calib.ux, "wx" = calib.nc.w, "ex" = calib.ex)
calib.all = calib.all[order(calib.all$score),]
rownames(calib.all) = 1:dim(calib.all)[1]

########################################
## test fold
########################################
test.X = data$X[test.ind,] 
test.Y1 = data$Y1[test.ind]
test.ex = predict(e.model, newdata=test.X, num.threads=1)$predictions
test.lx = hat.p * (1 + 1/Gamma * (1-test.ex)/test.ex)
test.ux = hat.p * (1 + Gamma* (1-test.ex)/(test.ex))
test.pred = predict(t.mdl, test.X, quantile=c(alpha/2, 1-alpha/2)) 
n_test = length(test.Y1)


###################################
# the PAC-type algorithm  
###################################

cat(" - Computing the PAC-type algorithm...")

calib.ex = calib.all$ex
calib.score = calib.all$score
# compute the 1-alpha quantile of hat{G}_n(t)
rand_ind = sample(n_calib)
v.kstar = ate.single.pac(Gamma, hat.p, calib.ex, calib.score, delta, alpha, rand_ind)

# conformal prediction interval
c.test.lo = test.pred[,1] - v.kstar
c.test.hi = test.pred[,2] + v.kstar
c.cover = (c.test.lo <= test.Y1) * (c.test.hi >= test.Y1)

########################################
# output summary of test   
########################################
res = data.frame("c.cov" = mean(c.cover), "c.len" = mean(c.test.hi-c.test.lo),
                 "n_train" = length(train.ind), "n_calib"=n_calib,
                 "gamma" = Gamma, "alpha" = alpha, "seed" = seed, "method" = "PAC") 

save.path = paste(out_dir, "syn_pred_pac_alpha_",alpha_ind,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep='')
write.csv(res, save.path)
