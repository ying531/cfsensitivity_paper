#!/usr/bin/env Rscript
########################################
## input configurations
########################################
args <- commandArgs(trailingOnly = TRUE)
p <- as.integer(args[1])
n <- as.integer(args[2])
alpha_ind <- as.integer(args[3])
Gamma_ind <- as.integer(args[4])
seed <- as.integer(args[5])

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
source("../utils/util_ate.R")
cat(paste(" - Running the script with marginally-valid algorithm and ground truth, alpha ", alpha, ", Gamma ",Gamma,
          ", n ", n, ", p ", p, ", seed ", seed, "\n"), sep = '')

########################################
## Output direcroty
########################################
if(!dir.exists("../results")){
  dir.create("../results")
}
out_dir <- "../results/simulation/"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

########################################
## Parameter
########################################
alpha0 = 0
n_test = 500
beta = matrix(c(-0.531,0.126,-0.312,0.018,rep(0,p-4)), nrow=p)
# generate true probability of treatment
pp = mean(data.gen.ate(n*1000,p,Gamma,beta,alpha0,obs=FALSE)$T)


set.seed(seed)
########################################
## fit on the training fold
########################################
train.data = data.gen.ate(n,p,Gamma,beta,alpha0,obs=TRUE)
train.X = (train.data$X[train.data$T==1,])[1:n,]
train.Y = (train.data$Y1[train.data$T==1])[1:n]
# train the nonconformity score function
train.score = conform.score(train.X, train.Y, "cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model

########################################
## calibration 
########################################
calib.data = data.gen.ate(n,p,Gamma,beta,alpha0,obs=TRUE)
calib.X = (calib.data$X[calib.data$T==1,])[1:n,]
calib.Y = (calib.data$Y1[calib.data$T==1])[1:n]
calib.ex = (calib.data$ex[calib.data$T==1])[1:n]
n_calib = length(calib.Y)

# lower and upper bounds of weight function
calib.lx = pp * (1 + (1-calib.ex) / (calib.ex*Gamma))
calib.ux = pp * (1 + Gamma * (1-calib.ex) / (calib.ex))
calib.nc.w = pp / (calib.data$ex[calib.data$T==1])[1:n]
# non-conformity score on calibration data
calib.score = conform.score(calib.X, calib.Y, "cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.all = data.frame("score"=calib.score, "lx"=calib.lx, "ux"=calib.ux, "wx"=calib.nc.w, "ex" = calib.ex)
calib.all = calib.all[order(calib.all$score),]
rownames(calib.all) = 1:dim(calib.all)[1]

########################################
## generate test fold
########################################
test.data = data.gen.ate(n_test,p,Gamma,beta,alpha0,obs=FALSE)
test.X = test.data$X
test.Y1 = test.data$Y1
test.ex = test.data$ex
test.lx = pp*(1+ 1/Gamma * (1-test.ex)/test.ex)
test.ux = pp*(1+ Gamma* (1-test.ex)/(test.ex))
test.pred = predict(t.mdl, test.X, quantile=c(alpha/2, 1-alpha/2)) 

########################################
## the confounding-aware algorithm  
########################################

cat(" - Computing the robust weighted conformal inference...")

# partial sums for confounding-aware
sum.num = rep(0,n_calib) # for numerator
sum.den = rep(0,n_calib) # for denominator
sum.num[1] = calib.all$lx[1]
sum.den[1] = calib.all$lx[1] + sum(calib.all$ux[2:n_calib])
for (k in 2:n_calib){
  sum.num[k] = sum.num[k-1] + calib.all$lx[k]
  sum.den[k] = sum.den[k-1] - calib.all$ux[k] + calib.all$lx[k] 
}

# confounding-aware prediction 
c.test.lo = rep(0,n_test)
c.test.hi = rep(0,n_test)
for (ii in 1:n_test){
  ratios = sum.num / (sum.den + test.ux[ii])
  kstar = min(which(ratios>1-alpha))
  v.kstar = calib.all$score[kstar]
  c.test.lo[ii] = test.pred[ii,1]-v.kstar
  c.test.hi[ii] = test.pred[ii,2]+v.kstar
}
# evaluate coverage on test data
c.cover = (c.test.lo <= test.Y1) * (c.test.hi >= test.Y1)

cat("Done.\n")

########################################
## the confounding-unaware algorithm  
########################################

cat(" - Computing the vanilla weighted conformal inference...")

nc.sum = rep(0,n_calib)
nc.sum[1] = calib.all$wx[1]
for (k in 2:n_calib){
  nc.sum[k] = nc.sum[k-1] + calib.all$wx[k]
}

nc.test.lo = rep(0,n_test)
nc.test.hi = rep(0,n_test)
nc.test.weight = pp / test.ex

for (ii in 1:n_test){
  nc.ratios = nc.sum / (nc.sum[n_calib] + nc.test.weight[ii])
  nc.kstar = min(which(nc.ratios>1-alpha))
  nc.v.kstar = calib.all$score[nc.kstar]
  nc.test.lo[ii] = test.pred[ii,1] - nc.v.kstar
  nc.test.hi[ii] = test.pred[ii,2] + nc.v.kstar
}
nc.cover = (nc.test.lo <= test.Y1) * (nc.test.hi >= test.Y1)

cat("Done.\n")

########################################
# output summary of test   
########################################
res = data.frame("c.cov" = mean(c.cover), "c.len" = mean(c.test.hi-c.test.lo),
                 "nc.cov" = mean(nc.cover), "nc.len" = mean(nc.test.hi-nc.test.lo), 
                 "n" = n, "p" = p, "n_calib" = n_calib,
                 "gamma" = Gamma, "alpha" = alpha, "seed" = seed, "method" = "marginal")

save.path = paste(out_dir, "pred_marginal_p_",p,"_n_",n,"_alpha_",alpha_ind,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep='')
write.csv(res, save.path)