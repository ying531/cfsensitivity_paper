library(grf)

args <- commandArgs(trailingOnly = TRUE)
Gamma_ind <- as.integer(args[1]) # 1:6
diff_ind <- as.integer(args[2]) # 1:4
diff_cf <- as.integer(args[3]) # 1:random 2:fixed
seed <- as.integer(args[4]) # 1:2

alpha = 0.1 # coverage target 1-alpha
gammas = seq(1,2, by=0.2) 
Gamma = gammas[Gamma_ind]
diffs = seq(-1,1, by=0.5) 
diff.ite = diffs[diff_ind] # effect size
diff.cf = c(TRUE, FALSE)
dcf = diff.cf[diff_cf] # whether random ITE

n = 2000 # fixed sample size
p = 4 # fixed covariate dimension
beta = matrix(c(-0.531,0.126,-0.312,0.018), nrow=p)

# load util functions
source("../utils/util_att.R")
if (dcf){
  cat(paste(" - Running the script with marginally-valid algorithm and ground truth, Gamma", Gamma, 
            ", effect size", diff.ite, ", random ITE, seed", seed, "\n"), sep = '')
}else{
  cat(paste(" - Running the script with marginally-valid algorithm and ground truth, Gamma", Gamma, 
            ", effect size", diff.ite, ", fixed ITE, seed", seed, "\n"), sep = '')
}

alpha0 = 1
n_test = 500
pp = mean(data.gen.att(n*1000, p, Gamma, beta, alpha0, obs=FALSE, a=diff.ite, diff.confound=dcf)$T)


set.seed(seed)
train.data = data.gen.att(n, p, Gamma, beta, alpha0, obs=TRUE, a=diff.ite, diff.confound=dcf)
train.X = (train.data$X[train.data$T==0,])[1:n,]
train.Y = (train.data$Y0[train.data$T==0])[1:n]
# train the nonconformity score function
train.score = conform.score(train.X, train.Y, "one-side-cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model

# generate calibration fold
calib.data = data.gen.att(n,p,Gamma,beta,alpha0,obs=TRUE, a=diff.ite, diff.confound=dcf)
calib.X = (calib.data$X[calib.data$T==0,])[1:n,]
calib.Y = (calib.data$Y0[calib.data$T==0])[1:n]
calib.ex = (calib.data$ex[calib.data$T==0])[1:n]
calib.nc.w = (1-pp)* calib.ex / (pp * (1-calib.ex))
n_calib = length(calib.Y)
  
# non-conformity score on the calibration set
calib.score = conform.score(calib.X, calib.Y, "one-side-cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.all = data.frame("score"=calib.score, "ex" = calib.ex, "wx" = calib.nc.w)
# re-order by the magnitude of scores
calib.all = calib.all[order(calib.all$score),]
calib.score = calib.all$score
calib.ex = calib.all$ex
  
# generate test fold
test.data = data.gen.att(n_test,p,Gamma,beta,alpha0,obs=FALSE, a=diff.ite, diff.confound=dcf)
test.X = test.data$X
test.Y0 = test.data$Y0
test.Y1 = test.data$Y1
test.T = test.data$T
test.ex = test.data$ex
test.pred = predict(t.mdl, test.X, quantile=1-alpha) 
  
true.ite = test.Y1 - test.Y0
  
################################################## 
### counterfactual prediction to check validity ##
################################################## 

cat(" - Computing the counterfactual prediction...")
  
# partial sums for confounding-aware
calib.lx = (1-pp)*calib.ex/(Gamma*pp*(1-calib.ex))
calib.ux = (1-pp)*calib.ex*Gamma/(pp*(1-calib.ex))
sum.num = rep(0,n_calib) # for numerator
sum.den = rep(0,n_calib) # for denominator
sum.num[1] = calib.lx[1]
sum.den[1] = calib.lx[1] + sum(calib.ux[2:n_calib])
for (k in 2:n_calib){
  sum.num[k] = sum.num[k-1] + calib.lx[k]
  sum.den[k] = sum.den[k-1] - calib.ux[k] + calib.lx[k]
}

test.lx = (1-pp) * test.ex / (Gamma*pp*(1-test.ex)) 
test.ux = (1-pp) * Gamma * test.ex / (pp*(1-test.ex)) 
test.pred = predict(t.mdl, test.X, quantile=1-alpha)

# confounding-aware prediction
c.test = rep(0,n_test)
for (ii in 1:n_test){
  ratios = sum.num / (sum.den + test.ux[ii])
  kstar = min(which(ratios>1-alpha))
  v.kstar = calib.score[kstar]
  c.test[ii] = test.pred[ii]+v.kstar
}
c.cover = test.Y0 <= c.test
  
# confounding-unaware prediction
nc.sum = rep(0,n_calib)
nc.sum[1] = calib.all$wx[1]
for (k in 2:n_calib){
  nc.sum[k] = nc.sum[k-1] + calib.all$wx[k]
}
nc.test = rep(0,n_test) 
nc.test.weight = (1-pp) * test.ex / (pp*(1-test.ex))
  
for (ii in 1:n_test){
  nc.ratios = nc.sum / (nc.sum[n_calib] + nc.test.weight[ii])
  nc.kstar = min(which(nc.ratios>1-alpha))
  nc.v.kstar = calib.score[nc.kstar]
  nc.test[ii] = test.pred[ii] + nc.v.kstar
}
nc.cover = test.Y0 <= nc.test

cat("Done.\n")
  
############################################# 
### sensitivity analysis invert prediction ##
############################################# 
  
cat(" - Running the sensitivity analysis...")

hat.gamma = rep(1, n_test)
for (ii in 1:n_test){
  hat.gamma[ii] = att.find.gamma(test.pred[ii], test.Y1[ii], test.ex[ii], pp, calib.ex, calib.score, alpha)
} 

cat("Done.\n")

res = data.frame("c.cov" = mean(c.cover[test.T==1]), "nc.cov" = mean(nc.cover[test.T==1]), 
                 "gamma"=Gamma, "ite" = diff.ite, "seed"=seed, "method" = "marginal")

write.csv(hat.gamma, paste("./results/sens_mgn_gvalue_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))
write.csv(true.ite, paste("./results/sens_mgn_gvalue_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))

write.csv(res, paste(".results/sens_mgn_res_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))