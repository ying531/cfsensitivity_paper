


syn.gen <- function(n.syn, pop.data, Gamma, setting=1){
  id.syn = sample(length(pop.data$mu0), n.syn, TRUE)
  mu.0 = pop.data$mu0[id.syn]
  ex.syn = pop.data$p.score[id.syn]
  X.syn = pop.data$X[id.syn,]
  rx.syn = pop.data$r.x[id.syn]
  # add confounder
  sig.U = rx.syn 
  U = rnorm(n.syn) * sig.U
  flip = rbinom(n.syn, size=1, prob=0.5) *2-1
  Y1 = mu.0 + pop.data$tau.x[id.syn] + abs(U) * flip
  Y0 = mu.0 + sig.U *( 3- abs(U)/sig.U) * (1-flip) 
  
  # confound to different levels
  # confounded propensity score
  p.thd.syn = (1/(ex.syn + (1-ex.syn)/Gamma ) -1) / ( 1/(ex.syn + (1-ex.syn)/Gamma) - 1/(ex.syn + Gamma*(1-ex.syn)) )
  t.syn = qnorm(1 - p.thd.syn/2) * sig.U
  exu.syn = (ex.syn/(ex.syn + Gamma*(1-ex.syn))) * (abs(U)>t.syn) + (ex.syn/(ex.syn+ (1-ex.syn)/Gamma)) * (abs(U)<=t.syn)
  T = rbinom(n.syn, size=1, prob=exu.syn)
  # observations
  Y = Y1 * T + Y0 * (1-T)
  syn.data = list("X" = X.syn, "ps.x" = ex.syn, "ps.xu" = exu.syn, "Y1" = Y1, "Y0" = Y0, "Y" = Y, "T" = T, 
                  "Gamma" = Gamma, "U" = U,
                  "mu0" = mu.0, "mu1" = mu.0 + pop.data$tau.x[id.syn], "taux" = pop.data$tau.x[id.syn], "rx" = rx.syn)
  
  return(syn.data)
}