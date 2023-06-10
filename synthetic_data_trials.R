#################### EM for estimating the MLE #########################

##################### Trial on synthetic data ##########################

trials = 10

n=500
p=9

#load initial values

alpha = 1
beta = matrix(runif(p,0,1))
gamma = matrix(runif(p,0,1))
delta = 1
sigma_sq = 0.5
tau = 0.2

U = matrix(rbern(n,0.5))
X = matrix(rnorm(n*p, 0, 1),ncol=9)
W = matrix(rbern(n,inv.logit(X %*% gamma + alpha*U)))
Y = matrix(tau*W + X %*% beta + delta*U + rnorm(n,0,sqrt(sigma_sq)))

lkhds = c()

beta = matrix(runif(p,0,1))
gamma = matrix(runif(p,0,1))
sigma_sq = runif(1,0,1)
tau = runif(1,0,1)

#the for-loop can be replaced by one call of the optimizer function

for (i in 1:trials){
  
  lkhds = c(lkhds,log_lkhd(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau))
  
  probs = u_probs(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau)
  
  # beta update
  
  beta = beta_update(X,Y,W,alpha,gamma,delta,tau,probs)
  
  lkhds = c(lkhds,log_lkhd(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau))
  
  probs = u_probs(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau)
  
  # gamma update
  
  gamma = gamma_update(X,Y,W,alpha,beta,delta,tau,probs)
  
  lkhds = c(lkhds,log_lkhd(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau))
  
  probs = u_probs(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau)
  
  # sigma_sq update
  
  sigma_sq = sigma_update(X,Y,W,alpha,beta,gamma,delta,tau,probs)
  
  lkhds = c(lkhds,log_lkhd(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau))
  
  probs = u_probs(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau)
  
  # tau update
  
  tau = tau_update(X,Y,W,alpha,beta,gamma,delta,probs)
}

