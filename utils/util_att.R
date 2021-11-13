############################################################
#####      data generating process of Sec6.2       #########
# n: output sample size 
# p: covariate dimension 
# Gamma: confounding level
# beta, alpha0: linear coefficients
# obs = TRUE generate observations (training data) of size n
# obs = FALSE generate all data of size n
# diff.confound = TRUE: random ITE
# diff.confound = FALSE: fixed ITE
# a: effect size
############################################################

data.gen.att <- function(n, p, Gamma, beta=NULL, alpha0=0, obs=TRUE, a, diff.confound=FALSE){
  X = matrix(runif(n*p),nrow=n,ncol=p)
  U = rnorm(n) * abs(1+0.5*sin(2.5*X[,1]))
  Y0 = X %*% beta + U
  if (diff.confound){
    Y1 = Y0 + a * U
  }else{
    Y1 = Y0 + a
  }
  prop.x = exp(alpha0 + X%*%beta)/(1 + exp(alpha0+X%*%beta))
  l.ex = prop.x / (prop.x + Gamma*(1-prop.x))
  u.ex = prop.x / (prop.x + (1-prop.x)/Gamma)
  if (sum(l.ex == u.ex) == length(l.ex)){
    p.x = 0
  }else{
    p.x = (u.ex - prop.x) / (u.ex - l.ex)
  }
  
  t.x = qnorm( p.x ) * abs(1+0.5*sin(2.5*X[,1]))
  prop.xu = (l.ex)*(U<=t.x) + (u.ex)*(U>t.x)
  TT = rbinom(n, size=1, prob=prop.xu)
  
  if (obs==FALSE){
    return(list("T"=TT, "X"=X, "U"=U, "Y0"=Y0, "Y1"=Y1, "ex"=prop.x, "exu"=prop.xu))
  }else{
    n_useful = sum(1-TT)
    while (n_useful < n){
      add.data = data.gen.att(n, p, Gamma, beta, alpha0, FALSE, a, diff.confound)
      X = rbind(X, add.data$X)
      U = c(U, add.data$U)
      Y0 = c(Y0, add.data$Y0)
      Y1 = c(Y1, add.data$Y1)
      prop.x = c(prop.x, add.data$ex)
      prop.xu = c(prop.xu, add.data$exu)
      TT = c(TT, add.data$T)
      n_useful = sum(1-TT)
    }
    return(list("T"=TT, "X"=X, "U"=U, "Y0"=Y0, "Y1"=Y1, "ex"=prop.x, "exu"=prop.xu))
  
  }
}

############################################################
#####   output nonconformity score and trained model  ######
# X: covariate matrix
# Y: response
# train and output trained model if trained_model = NULL
# otherwise, use trained_model to get nonconformity score
# quantile is the coverage (1-alpha)
############################################################

conform.score <- function(X,Y, method='cqr', trained_model = NULL, quantile=0.9){
  if (method == 'cqr'){
    if (is.null(trained_model)){
      trained_model = quantile_forest(X,Y,num.threads = 1)
    }
    # fit quantiles
    qs = predict(trained_model, X,quantile=c((1-quantile)/2, 1-(1-quantile)/2))
    q_lo = qs[,1]
    q_hi = qs[,2]
    score = pmax( Y-q_hi, q_lo-Y )
  }else if (method == 'one-side-cqr'){
    if (is.null(trained_model)){
      trained_model = quantile_forest(X,Y,num.threads = 1)
    }
    # fit quantiles
    qs = predict(trained_model, X,quantile=quantile)
    score = Y - qs
  }
  return(list("model"=trained_model, "score"=score))
}


#################################### 
# compute WS-R lower bound for cdf  
# at a single point V[i]  
####################################

wsr.cdf.single <- function(delta, lx, ux, M, i, rand_given, rand_ind){
  n = length(lx)
  if (i == n){
    lxx = lx/M
    uxx = rep(1, n)
  }else{
    lxx = c(lx[1:i], rep(0,n-i))/M
    uxx = 1 + (c(rep(0,i),-ux[(i+1):n]))/M
  }
  l.mean = mean(lxx) * M
  u.mean = mean(uxx) * M + 1 - M
  if (rand_given){
    wsr1 = wsr_lower(delta/2, lxx[rand_ind])*M
    wsr2 = wsr_lower(delta/2, uxx[rand_ind])*M +1-M
  }else{
    wsr1 = wsr_lower(delta/2, sample(lxx))*M
    wsr2 = wsr_lower(delta/2, sample(uxx))*M +1-M
  }
  return( min( max(wsr1,wsr2), max(l.mean, u.mean)))
}

#################################### 
# compute WS-R lower bound for cdf  
# at all points V[i]  
####################################

wsr.cdf <- function(delta, lx, ux, M, rand_given, rand_ind){
  n = length(lx)
  w_all <- sapply(1:n, FUN = wsr.cdf.single, delta=delta, lx=lx, ux=ux, M=M, rand_given = rand_given, rand_ind=rand_ind)
  return(w_all)
}

#################################### 
# find the hat{v} of pac procedure  
# for one single Gamma  
####################################

att.single.pac <- function(Gamma, pp, calib.ex, calib.score, delta, alpha, rand_ind){
  # lower and upper bounds
  calib.lx = (1-pp)*calib.ex/(Gamma*pp*(1-calib.ex))
  calib.ux = (1-pp)*calib.ex*Gamma/(pp*(1-calib.ex))
  # aggregate
  calib.all = data.frame("score"=calib.score, "lx"=calib.lx, "ux"=calib.ux, "ex" = calib.ex)
  # compute the 1-alpha quantile
  M = max(calib.ux)
  eta.gamma = wsr.qtl(delta, calib.all, M, alpha, rand_ind)
  return(eta.gamma)
}


############################################# 
# bi-search for quantile of wsr lower bound  
#############################################

wsr.qtl <- function(delta, calib.all, M, alpha, rand_ind){
  n_calib = nrow(calib.all)
  l.i = 1
  r.i = n_calib - 1
  m.i = floor((l.i+r.i)/2)
  left.cdf = wsr.cdf.single(delta, calib.all$lx, calib.all$ux, M, l.i, rand_given=TRUE, rand_ind)
  right.cdf = wsr.cdf.single(delta, calib.all$lx, calib.all$ux, M, r.i, rand_given=TRUE, rand_ind)
  right.cdf.check = 1 - calib.all$ux[n_calib]/n_calib
  if (right.cdf < 1-alpha){
    return(Inf)
  }
  mid.cdf = wsr.cdf.single(delta, calib.all$lx, calib.all$ux, M, m.i, rand_given=TRUE, rand_ind)
  gap = min(m.i - l.i, r.i - m.i)
  while( gap >0 ){
    if (mid.cdf < 1-alpha){
      l.i = m.i
      m.i = floor((l.i+r.i)/2)
      left.cdf = mid.cdf
      mid.cdf = wsr.cdf.single(delta, calib.all$lx, calib.all$ux, M, m.i, rand_given=TRUE, rand_ind)
    }else{
      r.i = m.i
      m.i = floor((l.i+r.i)/2)
      right.cdf = mid.cdf
      mid.cdf = wsr.cdf.single(delta, calib.all$lx, calib.all$ux, M, m.i, rand_given=TRUE, rand_ind)
    }
    gap = min(m.i - l.i, r.i - m.i)
    
  }
  if (mid.cdf >= 1-alpha){
    return(calib.all$score[m.i])
  }else{
    return(calib.all$score[r.i])
  }
}


##############
## invert for sensitivity analysis
mgn.inv <- function(ux, sum.num, sum.den, alpha){
  ratios = sum.num / (sum.den + ux)
  kstar = min(which(ratios>1-alpha))
  return(kstar)
}

att.Y <- function(Gamma, test.ex, pp, sum.num, sum.den, test.pred, calib.score, alpha){
  kstar.list = sapply((1-pp)*Gamma*test.ex/(pp*(1-test.ex)), mgn.inv, sum.num = sum.num, sum.den = sum.den, alpha = alpha)
  c.test = test.pred + calib.score[kstar.list]
  return(c.test)
}

################################################################
######### find the 1-alpha invertal (0, hat{Y}(gamma)] #########
# for one single test point 
# t.ex is the propensity score of test point
# pp is the population P(T=1)
# t.pred is the estimated (1-alpha) quantile for the test point
# calib.ex is the propensity scores of calibration set (sorted according to scores)
# c.score is the nonconformity scores of the calibration set
################################################################

att.single <- function(Gamma, t.ex, pp, t.pred, calib.ex, c.score, alpha){
  # lower and upper bounds
  # ATT type bound
  calib.lx = (1-pp)*calib.ex/(Gamma*pp*(1-calib.ex))
  calib.ux = (1-pp)*calib.ex*Gamma/(pp*(1-calib.ex))
  n_calib = length(calib.ex)
  sum.num = rep(0,n_calib) # for numerator
  sum.den = rep(0,n_calib) # for denominator
  sum.num[1] = calib.lx[1]
  sum.den[1] = calib.lx[1] + sum(calib.ux[2:n_calib])
  for (k in 2:n_calib){
    sum.num[k] = sum.num[k-1] + calib.lx[k]
    sum.den[k] = sum.den[k-1] - calib.ux[k] + calib.lx[k] 
  }
  ux = (1-pp)*Gamma* t.ex/(pp*(1- t.ex))
  ratios = sum.num / (sum.den + ux)
  kstar = min(which(ratios>1-alpha))
  if (kstar == Inf){
    return(Inf)
  }
  att = t.pred + c.score[kstar]
  return(att)
}


#################################################################
#### bisearch for the largest gamma s.t. Y1 >= hat{Y}(gamma) ####
################# with marginally valid procedure ###############
# t.pred: predicted quantile of one-sided CQR
# t.ex: propensity score of test point
# pp: population P(T=1)
# t.pred: estimated (1-alpha) quantile for the test point
# calib.ex: the propensity scores of calibration set (sorted according to scores)
# c.score: the nonconformity scores of the calibration set
#################################################################
att.find.gamma <- function(t.pred, Y1, t.ex, pp, calib.ex, c.score, alpha){
  if (Y1 - t.pred > max(c.score)){
    return(Inf)
  }
  lg = 1
  rg = 1
  l.y = att.single(rg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
  if (Y1 < l.y){
    return(0)
  }
  r.y = att.single(rg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
  while(Y1 >= r.y){
    rg = rg *2 
    r.y = att.single(rg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
  }
  mg = (lg + rg) / 2
  m.y = att.single(mg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
  gap = min(mg-lg, rg-mg)
  while(gap > 0.01){
    if (Y1 > m.y){
      lg = mg
      l.y = m.y
      mg = (lg + rg) / 2
      m.y = att.single(mg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
    }else{
      rg = mg
      r.y = m.y
      mg = (lg + rg) / 2
      m.y = att.single(mg, t.ex, pp, t.pred, calib.ex, c.score, alpha)
    }
    gap = min(mg-lg, rg-mg)
  }
  return(mg)
}

#################################################################
#### bisearch for the largest gamma s.t. Y1 >= hat{Y}(gamma) ####
#################     with PAC-type procedure     ###############
# t.pred: predicted quantile of one-sided CQR
# t.ex: propensity score of test point
# pp: population P(T=1)
# t.pred: estimated (1-alpha) quantile for the test point
# calib.ex: the propensity scores of calibration set (sorted according to scores)
# c.score: the nonconformity scores of the calibration set
#################################################################

att.find.gamma.pac <- function(t.pred, Y1, pp, calib.ex, calib.score, delta, alpha, rand_ind){
  lg = 1
  rg = 1
  l.Y = t.pred + att.single.pac(lg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  r.Y = t.pred + att.single.pac(rg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  while((r.Y < Y1)){
    rg = rg * 2
    r.Y = t.pred + att.single.pac(rg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  }
  # bi-search for the smallest gamma s.t. Y(1) < \hat{Y}(gamma)
  mg = (lg + rg) / 2
  m.y = t.pred + att.single.pac(mg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
  gap = min(mg-lg, rg-mg)
  while(gap > 0.01){
    if (Y1 > m.y){
      lg = mg
      l.y = m.y
      mg = (lg + rg) / 2
      m.y = t.pred + att.single.pac(mg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
    }else{
      rg = mg
      r.y = m.y
      mg = (lg + rg) / 2
      m.y = t.pred + att.single.pac(mg, pp, calib.ex, calib.score, delta, alpha, rand_ind)
    }
    gap = min(mg-lg, rg-mg)
  }
  return(mg)
}



######################################
# basic functions for WSR inequality  
######################################

compute_k_lower <- function(x, mu, nu){
  kterms = 1 + nu * (x - mu)
  ks = rep(kterms[1],length(kterms))
  for (ii in 2:length(ks)){
    ks[ii] = ks[ii-1]*kterms[ii]
  }
  return(max(ks))
}

compute_k_upper <- function(x, mu, nu){
  kterms = 1 - nu * (x - mu)
  ks = rep(kterms[1],length(kterms))
  for (ii in 2:length(ks)){
    ks[ii] = ks[ii-1]*kterms[ii]
  }
  return(max(ks))
}

wsr_lower <- function(delta, x){
  n = length(x)
  mu_hat <- (1/2 + cumsum(x)) / (1 : n + 1)
  sig_hat = (1/4 + cumsum((x-mu_hat)^2))/(1:n +1)
  # nu[i] = xxxx/simga_{i-1}
  nu <- pmin(1, sqrt(2 * log(1 / delta) / (n * sig_hat^2)))
  nu[2: length(nu)] = nu[1:(length(nu)-1)]
  nu[1] = pmin(1, sqrt(2 * log(1 / delta) / (n /4)))
  # u_list <- (1 : 1000) / 1000 
  u_list <- (0 : 1000) / 1000 
  k_all <- sapply(u_list, FUN = compute_k_lower, x = x, nu = nu)
  u_ind <- min(which(k_all <= 1 / delta))
  if (u_ind == Inf){
    u_ind = 1001
  }
  if (u_ind == 2){
    u_ind = 1 # if positive values in u_list satisfies, then set to zero. 
  }
  bnd <- u_list[u_ind]
  return(bnd)
}

wsr_upper <- function(delta, x){
  n = length(x)
  mu_hat <- (1/2 + cumsum(x)) / (1 : n + 1)
  sig_hat = (1/4 + cumsum((x-mu_hat)^2))/(1:n +1)
  nu <- pmin(1, sqrt(2 * log(1 / delta) / (n * sig_hat^2)))
  u_list <- (1 : 1000) / 1000 
  k_all <- sapply(u_list, FUN = compute_k_upper, x = x, nu = nu)
  u_ind <- min(which(k_all > 1 / delta))
  if (u_ind == Inf){
    u_ind = 1000
  }
  bnd <- u_list[u_ind]
  return(bnd)
}