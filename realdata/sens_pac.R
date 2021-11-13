#!/usr/bin/env Rscript
########################################
## input configurations
########################################
args <- commandArgs(trailingOnly = TRUE)
pos <- as.integer(args[1]) # whether test for positive ITE
seed <- as.integer(args[2]) 

########################################
## load libraries
########################################
suppressPackageStartupMessages(library(grf))
options(warn=-1)

########################################
## load util functions
########################################
source("../utils/util_att.R")
cat(paste(" - Running sensitivity analysis with PAC algorithm, seed", seed, "\n"), sep = '')

########################################
## Output direcroty
########################################
if(!dir.exists("../results")){
  dir.create("../results")
}
out_dir <- "../results/realdata/"
if(!dir.exists(out_dir)){
  dir.create(out_dir)
}

########################################
## Parameter
########################################
alpha = 0.1
delta = 0.05
set.seed(seed)

########################################
## Data processing
########################################
org.data = read.csv("./acic_data.csv")
n = nrow(org.data)
X = org.data %>% select(-Z,-Y) %>% as.matrix
if (pos == 2){ # test for negative ITE
  org.data$Y = -org.data$Y
}
re.ind = sample(n,n)
train.idx = re.ind[1:floor(n/3)]
# Z=0 as calibration set, Z=1 as test set
ct.idx = re.ind[(floor(n/3)+1):n]
calib.data = org.data[ct.idx,] %>% filter(Z==0)
test.data = org.data[ct.idx,] %>% filter(Z==1)
pp = mean(org.data$Z)

########################################
## fit on the training fold
########################################
train.X = org.data[train.idx,] %>% filter(Z==0) %>% select(-Z,-Y)
train.Y = (org.data[train.idx,] %>% filter(Z==0))$Y
train.score = conform.score(train.X, train.Y, "one-side-cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model
ps.model = regression_forest((org.data %>% select(-Z,-Y) %>% as.matrix)[train.idx,], org.data$Z[train.idx], num.threads=1)

########################################
## calibration 
########################################
calib.X = calib.data %>% select(-Z,-Y) %>% as.matrix
calib.Y = calib.data$Y
calib.T = calib.data$Z
n_calib = nrow(calib.X)
# non-conformity score on calibration data
calib.score = conform.score(calib.X, calib.Y, "one-side-cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.ex = predict(ps.model, newdata = calib.X)$predictions
calib.all = data.frame("score"=calib.score, "ex" = calib.ex)
calib.all = calib.all[order(calib.all$score),]
calib.score = calib.all$score
calib.ex = calib.all$ex

########################################
## test fold
########################################
test.X = test.data %>% select(-Z,-Y) %>% as.matrix
test.Y1 = test.data$Y
test.T = test.data$Z
test.ex = predict(ps.model, newdata=test.X)$predictions
test.pred = predict(t.mdl, newdata = test.X, quantile=1-alpha)
n_test = nrow(test.X)

############################################# 
### sensitivity analysis invert prediction 
############################################# 

cat(" - Running the sensitivity analysis...")

rand_ind = sample(n_calib)
# grid search of gamma values
g = 1
dd = att.single.pac(g, pp, calib.ex, calib.score, delta, alpha, rand_ind)
while (dd < Inf ){
  g = g * 2
  dd = att.single.pac(g, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  print(c(g,dd))
}
# grid search of gamma values
ggs = c(seq(1, 15, length.out = 200), seq(15, 25, length.out = 100))
rr = c()
for (g in ggs){
  rr = c(rr, att.single.pac(g, pp, calib.ex, calib.score, delta, alpha, rand_ind))
}
# find the smallest gamma s.t. Y(1) < \hat{Y}(gamma)
hat.gamma = rep(1, n_test)
for (ii in 1:n_test){
  if (test.pred[ii]+rr[1] > test.Y1[ii]){
    hat.gamma[ii] = 0
  }else{
    add = test.pred[ii] + rr
    idx = min(length(ggs), min(which(add > test.Y1[ii])))
    hat.gamma[ii] = ggs[idx]
  }
}  

cat("Done.\n")

########################################
## output summary of test
## along with original test fold
########################################
if (pos == 1){
  all.test = cbind(test.data, hat.gamma)
  write.csv(all.test, paste(out_dir, "sens_positive_pac_hgamma_seed_",seed,".csv",sep=''))
}else{
  org.data$Y = -org.data$Y
  all.test = cbind(test.data, hat.gamma)
  write.csv(all.test, paste(out_dir, "sens_negative_pac_hgamma_seed_",seed,".csv",sep=''))
}

