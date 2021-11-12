library(grf)

args <- commandArgs(trailingOnly = TRUE)
Gamma_ind <- as.integer(args[1]) 
diff_ind <- as.integer(args[2]) 
diff_cf <- as.integer(args[3]) # 1:random 2:fixed
seed <- as.integer(args[4]) 

alpha = 0.1 # coverage target 1-alpha
delta = 0.05 # confidence 
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
  cat(paste(" - Running the script with PAC algorithm and ground truth, Gamma", Gamma, 
            ", effect size", diff.ite, ", random ITE, seed", seed, "\n"), sep = '')
}else{
  cat(paste(" - Running the script with PAC algorithm and ground truth, Gamma", Gamma, 
            ", effect size", diff.ite, ", fixed ITE, seed", seed, "\n"), sep = '')
}


alpha0 = 1
n_test = 1000

set.seed(seed)
# generate training fold
train.data = data.gen.att(n, p, Gamma, beta, alpha0, obs=TRUE, a=diff.ite, diff.confound=dcf)
train.X = (train.data$X[train.data$T==0,])[1:n,]
train.Y = (train.data$Y0[train.data$T==0])[1:n]
# train the nonconformity score function
train.score = conform.score(train.X, train.Y, "one-side-cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model
# fit hat{e}(x)
e.model = regression_forest(train.data$X, train.data$T, num.threads = 1)
pp = mean(train.data$T)

# generate calibration fold
calib.data = data.gen.att(n, p, Gamma, beta, alpha0, obs=TRUE, a=diff.ite, diff.confound=dcf)
calib.X = (calib.data$X[calib.data$T==0,])[1:n,]
calib.Y = (calib.data$Y0[calib.data$T==0])[1:n]
calib.ex = predict(e.model, newdata=calib.X)$predictions
n_calib = length(calib.Y)

# lower and upper bounds of weight function
calib.lx = (1-pp)*calib.ex/(Gamma*pp*(1-calib.ex))
calib.ux = (1-pp)*calib.ex*Gamma/(pp*(1-calib.ex))

# non-conformity score on the calibration set
calib.score = conform.score(calib.X, calib.Y, "one-side-cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.all = data.frame("score"=calib.score, "ex" = calib.ex)
calib.all = calib.all[order(calib.all$score),]
rownames(calib.all) = 1:dim(calib.all)[1]

# generate test fold
test.data = data.gen.att(n_test,p,Gamma,beta,alpha0,obs=FALSE, a=diff.ite, diff.confound=dcf)
test.X = test.data$X
test.T = test.data$T
test.Y0 = test.data$Y0
test.Y1 = test.data$Y1
test.pred = predict(t.mdl, test.X, quantile=1-alpha) 

true.ite = test.Y1 - test.Y0

calib.ex = calib.all$ex
calib.score = calib.all$score
  
################################################## 
### counterfactual prediction to check validity ##
################################################## 

cat(" - Computing the counterfactual prediction...")

rand_ind = sample(n_calib)
c.test = test.pred + att.single.pac(Gamma, pp, calib.ex, calib.score, delta, alpha, rand_ind)
c.cover = c.test >= test.Y0

cat("Done.\n")

############################################# 
### sensitivity analysis invert prediction ##
############################################# 

cat(" - Running the sensitivity analysis...")

g = 1
dd = att.single.pac(g, pp, calib.ex, calib.score, delta, alpha, rand_ind)
while (dd < Inf){
  g = g * 2
  dd = att.single.pac(g, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  print(c(g, dd))
}
# grid search of gamma values
ggs = c(seq(1, g/2, length.out = 200), seq(g/2, g, length.out = 100))
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
    idx = min(which(add > test.Y1[ii]))
    hat.gamma[ii] = ggs[idx]
  }
}  

cat("Done.\n")

res = data.frame("c.cov" = mean(c.cover[test.T==1]), "set" = diff_cf, "gamma"=Gamma, "ite" = diff.ite, "seed"=seed, "method" = "pac_est") 

write.csv(hat.gamma, paste("./results/sens_pac_est_hgamma_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))
write.csv(true.ites, paste("./results/sens_pac_est_ite_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))

write.csv(res, paste("./results/sens_whp_est_res_diff_",diff_ind,"_set_",diff_cf,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep=''))