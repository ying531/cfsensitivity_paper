############################################################
#####      data generating process of Sec6.1       #########
# n: output sample size 
# p: covariate dimension 
# Gamma: confounding level
# beta, alpha0: linear coefficients
# obs = TRUE generate observations (training data) of size n
# obs = FALSE generate all data of size n
############################################################

data.gen.ate <- function(n, p, Gamma, beta, alpha0=0, obs=TRUE){
  X = matrix(runif(n*p),nrow=n,ncol=p)
  U = rnorm(n) * abs(1+0.5*sin(2.5*X[,1]))
  Y1 = X %*% beta + U
  prop.x = exp(alpha0 + X%*%beta)/(1+exp(alpha0+X%*%beta))
  p.x = (1/(prop.x + (1-prop.x)/Gamma ) -1)/( 1/(prop.x + (1-prop.x)/Gamma) - 1/(prop.x + Gamma*(1-prop.x)))
  t.x = qnorm(1-p.x/2) * abs(1+0.5*sin(2.5*X[,1]))
  prop.xu = (prop.x/(prop.x+Gamma*(1-prop.x)))*(abs(U)>t.x) + (prop.x/(prop.x+ (1-prop.x)/Gamma))*(abs(U)<=t.x)
  TT = rbinom(n, size=1, prob=prop.xu)
  
  if (obs==FALSE){
    return(list("T"=TT, "X"=X, "U"=U, "Y1"=Y1, "ex"=prop.x, "exu"=prop.xu))
  }else{
    n_useful = sum(TT)
    while (n_useful < n){
      add.data = data.gen.ate(n,p,Gamma,beta,alpha0,FALSE)
      X = rbind(X, add.data$X)
      U = c(U, add.data$U)
      Y1 = c(Y1, add.data$Y1)
      prop.x = c(prop.x, add.data$ex)
      prop.xu = c(prop.xu, add.data$exu)
      TT = c(TT, add.data$T)
      n_useful = sum(TT)
    }
    return(list("T"=TT, "X"=X, "U"=U, "Y1"=Y1, "ex"=prop.x, "exu"=prop.xu))
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
conform.score <- function(X, Y, method='cqr', trained_model = NULL, quantile=0.9){
  if (method == 'cqr'){
    if (is.null(trained_model)){
      trained_model = quantile_forest(X,Y,num.threads = 1)
    }
    # fit quantiles
    qs = predict(trained_model, X,quantile=c((1-quantile)/2, 1-(1-quantile)/2))
    q_lo = qs[,1]
    q_hi = qs[,2]
    score = pmax( Y-q_hi, q_lo-Y )
  }
  return(list("model"=trained_model, "score"=score))
}


#################################### 
# compute WS-R lower bound for cdf  
# at a single point V[i]  
####################################
wsr.cdf.single <- function(delta, lx, ux, M, i, rand_given, rand_ind){
  n = length(lx)
  lxx = c(lx[1:i], rep(0,n-i))/M
  uxx = 1 + (c(rep(0,i),-ux[(i+1):n]) )/M
  if (rand_given){
    wsr1 = wsr_lower(delta/2, lxx[rand_ind])*M
    wsr2 = wsr_lower(delta/2, uxx[rand_ind])*M +1-M
  }else{
    wsr1 = wsr_lower(delta/2, sample(lxx))*M
    wsr2 = wsr_lower(delta/2, sample(uxx))*M +1-M
  }
  return( max(wsr1,wsr2))
}

#################################### 
# compute WS-R lower bound for cdf  
# at all points V[i]  
####################################

wsr.cdf <- function(delta, lx, ux, M, rand_given, rand_ind){
  n = length(lx)
  w_all <- sapply(1:n, FUN = wsr.cdf.single, delta=delta, lx=lx, ux=ux, 
                  M=M, rand_given = rand_given, rand_ind=rand_ind)
  return(w_all)
}

#################################### 
# find the hat{v} of pac procedure  
# for one single Gamma  
####################################

ate.single.pac <- function(Gamma, pp, calib.ex, calib.score, delta, alpha, rand_ind){
  # lower and upper bounds for likelihood ratio
  calib.lx = pp*(1+ (1-calib.ex)/(calib.ex*Gamma))
  calib.ux = pp*(1+ Gamma* (1-calib.ex)/(calib.ex))
  # aggregate information
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
  u_list <- (1 : 1000) / 1000 
  k_all <- sapply(u_list, FUN = compute_k_lower, x = x, nu = nu)
  u_ind <- min(which(k_all <= 1 / delta))
  if (u_ind == Inf){
    u_ind = 1000
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
