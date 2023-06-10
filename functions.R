library(glmnet)
library(tidyverse)
library(boot)
library(Rlab)


# Given a table with the covariates X, the treatments W and the outcomes Y,
# the code uses EM with CAVI updates to estimate the MLE.


# expected data dimensions for all functions
# X - nxp matrix
# Y - nx1 matrix
# W - nx1 binary matrix
# alpha - real number
# beta - nx1 matrix
# gamma - nx1 matrix
# delta - real number
# tau - real number
# sigma_sq - real, non-negative number

# ratio_to_prob is only there to make the posterior probabilities computationally
# stable.

ratio_to_prob<- function(r){
  
  if(r>1e10 || is.na(r)){return(1)}
  else{return(r/(1+r))}
}


# u_probs calculates the posterior probabilities of the latent variables being
# equal to 1.

u_probs <- function(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau){
  
  n = nrow(X)
  p = ncol(X)
  
  ratios = exp((2*(Y-tau*W - X %*% beta)-delta)*delta/(2*sigma_sq))
  ratios = ratios * (1+exp(X %*% gamma))
  ratios = ratios * exp(alpha*W)
  ratios = ratios /  (1+ exp(alpha + X %*% gamma))
  
  new_probs = c()
  
  for(i in 1:n){
    new_probs = c(new_probs,ratio_to_prob(ratios[i]))
  }
  
  return(new_probs)
}

# beta_update is the CAVI update for beta holding all other parameters fixed.

beta_update<- function(X,Y,W,alpha,gamma,delta,tau,probs){
  
  outcomes = Y - tau*W - delta*probs
  
  new_df = data.frame(X,outcomes)
  
  new_beta = coef(lm(outcomes~.-1,data=new_df))
  
  return(unname(new_beta))
}

# gamma_update is the CAVI update for gamma holding all other parameters fixed.

gamma_update<- function(X,Y,W,alpha,beta,delta,tau,probs){
  
  n = nrow(X)
  
  wts = c(probs,1-probs)
  offsets = c(rep(alpha,n),rep(0,n))
  
  new_X = rbind(X,X)
  
  new_W = matrix(W,ncol=1)
  new_W = rbind(W,W)
  
  new_df = data.frame(cbind(new_W,new_X))
  colnames(new_df)[1] = "W"
  
  new_gamma = coef(glm(W~.-1, family=binomial,data = new_df,
                       weights = wts, offset = offsets ))
  
  return(unname(new_gamma))
}


# sigma_update is the CAVI update for the variance \sigma^2.


sigma_update <- function(X,Y,W,alpha,beta,gamma,delta,tau,probs){
  
  n = nrow(X)
  
  a0 = (Y - tau*W - X %*% beta)^2
  a1 = (Y - tau*W - X %*% beta - delta)^2
  
  new_sigma = (sum(probs*a1) + sum((1-probs)*a0))/n
  
  return(new_sigma)
}

# tau_update is the CAVI update for the tau

tau_update<- function(X,Y,W,alpha,beta,gamma,delta,probs){
  
  n_treat = sum(W)
  
  Y_treat = Y[which(W==1)]
  X_treat = X[which(W==1),]
  probs_treat = probs[which(W==1)]
  
  new_tau = sum(Y_treat - X_treat %*% beta - delta*probs_treat)/n_treat
  
  return(new_tau)
}


# log_lkhd calculates the log likelihood of the observed data for given values
# of the parameters.


log_lkhd<-function(X,Y,W,alpha,beta,gamma,delta,sigma_sq,tau){
  
  a0 = Y - tau*W - X %*% beta
  a1 = a0 - delta
  
  a0 = dnorm(a0,mean=0,sd=sqrt(sigma_sq))
  a1 = dnorm(a1,mean=0,sd=sqrt(sigma_sq))
  
  b0 = inv.logit(X %*% gamma)
  b1 = inv.logit(X %*% gamma + alpha)
  
  b0 = dbern(W,b0)
  b1 = dbern(W,b1)
  
  lkhd = sum(log((b0*a0 + b1*a1)/2))
  
  return(lkhd)
}


# optimizer returns the optimized values of the parameters and the sequence of 
# log-likelihoods to check that the likelihood increases with each step.
# Calling this function usually gives some warnings. That is because I am running
# a weighted logistic regression and the glm function is fussy about getting
# non-integer weights. It doesn't cause any issues.


optimizer <- function(X,Y,W,alpha,delta,trials){
  
  n = nrow(X)
  p = ncol(X)
  
  beta = matrix(runif(p,0,1))
  gamma = matrix(runif(p,0,1))
  sigma_sq = runif(1,5,10)
  tau = runif(1,0,1)
  
  lkhds = c()
  
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
  
  return(list("beta" = beta,
              "gamma" = gamma,
              "sigma_sq" = sigma_sq,
              "tau" = tau,
              "likelihoods" = lkhds))
}


# r_sq_ratios returns the two R-squared values that we eventually plot in the
# Imbens procedure.


r_sq_ratios <-function(X,Y,W,alpha,delta,trials){
  
  n = nrow(X)
  
  optimized_0 = optimizer(X,Y,W,0,0,trials)
  optimized = optimizer(X,Y,W,alpha,delta,trials)
  
  Sigma = mean((Y-mean(Y))^2)
  
  r_sq_outcomes_0 = 1 - optimized_0$sigma_sq/Sigma
  r_sq_outcomes = 1 - optimized$sigma_sq/Sigma
  
  outcomes_ratio = (r_sq_outcomes-r_sq_outcomes_0)/(1-r_sq_outcomes_0)
  
  gamma_0 = matrix(optimized_0$gamma, ncol=1)
  gamma = matrix(optimized$gamma, ncol=1)
  
  x_cov = t(X) %*% X
  x_cov = x_cov/n
  
  r_sq_treatments_0 = 1 - (pi^2/3)/(t(gamma_0) %*% x_cov %*% gamma_0+pi^2/3)
  r_sq_treatments = 1 - (pi^2/3)/(t(gamma) %*% x_cov %*% gamma + alpha^2/4 + pi^2/3)
  
  treatments_ratio = (r_sq_treatments - r_sq_treatments_0)/(1-r_sq_treatments_0)
  
  return(list("outcomes" = outcomes_ratio,
              "treatments" = treatments_ratio))
  
}
