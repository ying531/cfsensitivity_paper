library(grf)

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

# load util functions
source("../utils/util_ate.R")
cat(paste(" - Running the script with PAC-type algorithm and ground truth, alpha ", alpha, ", Gamma ",Gamma,
          ", n ", n, ", p ", p, ", seed ", seed, "\n"), sep = '')

# PAC, estimated e(x)
alpha0 = 0
delta = 0.05
n_test = 10000
beta = matrix(c(-0.531,0.126,-0.312,0.018,rep(0,p-4)), nrow=p)
# generate true probability of treatment
pp = mean(data.gen.ate(n*1000,p,Gamma,beta,alpha0,obs=FALSE)$T)


set.seed(seed)
train.data = data.gen.ate(n,p,Gamma,beta,alpha0,obs=TRUE)
train.X = (train.data$X[train.data$T==1,])[1:n,]
train.Y = (train.data$Y1[train.data$T==1])[1:n]
# train the nonconformity score function
train.score = conform.score(train.X, train.Y, "cqr", trained_model=NULL, quantile=1-alpha)
t.mdl = train.score$model
# estimate hat{e}(x)
hat.p = mean(train.data$T)
e.model = regression_forest(train.data$X, train.data$T, num.threads = 1)
train.ex = predict(e.model)$predictions
M = max(hat.p * (1 + Gamma * (1-train.ex)/(train.ex)))
  
# generate calibration fold
calib.data = data.gen.ate(n,p,Gamma,beta,alpha0,obs=TRUE)
calib.X = (calib.data$X[calib.data$T==1,])[1:n,]
calib.Y = (calib.data$Y1[calib.data$T==1])[1:n]
calib.ex = predict(e.model, newdata=calib.X)$predictions
n_calib = length(calib.Y)

# lower and upper bounds of weight function
calib.lx = hat.p * (1 + (1-calib.ex) / (calib.ex*Gamma))
calib.ux = hat.p * (1 + Gamma * (1-calib.ex) / (calib.ex))

# non-conformity score on calibration data
calib.score = conform.score(calib.X, calib.Y, "cqr", trained_model=t.mdl, quantile=1-alpha)$score
calib.all = data.frame("score"=calib.score, "lx"=calib.lx, "ux"=calib.ux, "wx"=calib.nc.w, "ex" = calib.ex)
calib.all = calib.all[order(calib.all$score),]
rownames(calib.all) = 1:dim(calib.all)[1]

# calculate estimation error and actual gap
calib.exu.true = calib.data$exu[calib.data$T==1]
calib.ex.true = calib.data$ex[calib.data$T==1]
calib.wx.true = pp / calib.exu.true  
gap = max( mean(pmax(0, calib.lx- calib.wx.true )), mean(pmax(calib.wx.true - calib.ux, 0)) )
l.err = mean(abs(calib.lx - pp * (1 + (1-calib.ex.true)/(calib.ex.true*Gamma))))
u.err = mean(abs(calib.ux - pp * (1 + Gamma*(1-calib.ex.true)/(calib.ex.true))))

# generate test fold 
test.data = data.gen.ate(n_test,p,Gamma,beta,alpha0,obs=FALSE)
test.X = test.data$X
test.Y1 = test.data$Y1
test.ex = predict(e.model, newdata=test.X)$predictions
test.pred = predict(t.mdl, test.X, quantile=c(alpha/2, 1-alpha/2)) 
test.T = test.data$T


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

res = data.frame("c.cov" = mean(c.cover), "c.len" = mean(c.test.hi-c.test.lo),
                 "gap" = gap, "l.err" = l.err, "u.err" = u.err,
                 "n"=n, "p"=p, "n_calib" = n_calib, 
                 "gamma" = Gamma, "alpha" = alpha, "seed" = seed, "method" = "PAC_estimated") 

save.path = paste("./results/pred_pac_est_p_",p,"_n_",n,"_alpha_",alpha_ind,"_gamma_",Gamma_ind,"_seed_",seed,".csv",sep='')
write.csv(res, save.path)