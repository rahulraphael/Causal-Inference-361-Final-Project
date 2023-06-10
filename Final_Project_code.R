################## Trials on actual data ##################

# Remove the comments below and install the packages. You need them to get the
# Lalonde Dataset

#install.packages("devtools")
#devtools::install_github("jjchern/lalonde")

# Importing the datasets

lalonde::nsw-> nsw_data
lalonde::nsw_dw -> nsw_dw 
lalonde::psid_controls -> psid_controls
lalonde::psid_controls2 - psid_controls2
lalonde::psid_controls3-> psid_controls3
lalonde::cps_controls-> cps_controls
lalonde::cps_controls2-> cps_controls2
lalonde::cps_controls3-> cps_controls3

# Removes the unncessary first column.

nsw_dw = nsw_dw[,-c(1)]

# Separating covariates, outcomes and treatment assignments.

X = as.matrix(unname(nsw_dw[,2:8]))
W = matrix(unname(nsw_dw$treat))
Y = matrix(unname(nsw_dw$re78))

# trials is the number of iterations in the optimizer function.

trials = 300

# These steps find the estimates for alpha=delta=0.

optimized_0 = optimizer(X,Y,W,0,0,trials)

sigma_sq_0 = optimized_0$sigma_sq
gamma_0 = optimized_0$gamma

# Sigma is the sample variance of the outcomes.

Sigma = mean((Y-mean(Y))^2)

# x_cov is the sample covariance matrix of the covariates.

n = nrow(X)

x_cov = t(X-colMeans(X)) %*% (X-colMeans(X))
x_cov = x_cov/n


# The values of alpha and delta that reduce the estimated tau by roughly 1000
# follow the relationship alpha*delta = c.
# The vector prods contains a range of values of c that all reduce the estimated
# tau by around 1000.
# The vector alphas contains the range of values of alpha that we will try. The
# range can be made larger but it's unnecessary to increase it beyond 5. Also,
# making it too large makes the glm calls in the optimizer function unstable.

prods = seq(3500,4500,by=100)
q = length(prods)

alphas = seq(0.3,4.1,by=0.2)
a = length(alphas)

# outputs_df is the data frame that will store the alpha,delta values and the
# corresponding estimated tau and R-squared values.

outputs_df = matrix(0,ncol=5,nrow=a*q)



for (i in 1:a){
  for (j in 1:q){
    
    alpha = alphas[i]
    delta = min(prods[j]/alpha,1e4)
    
    optimized = optimizer(X,Y,W,alpha,delta,trials)
    
    gamma_1 = optimized$gamma
    sigma_sq_1 = optimized$sigma_sq
    
    outputs_df[(i-1)*q+j,1] = alpha
    outputs_df[(i-1)*q+j,2] = delta
    outputs_df[(i-1)*q+j,3] = optimized$tau
    
    r_sq_outcomes_0 = 1 - sigma_sq_0/Sigma
    r_sq_outcomes_1 = 1 - sigma_sq_1/Sigma
    
    outputs_df[(i-1)*q+j,4] = (r_sq_outcomes_1 - r_sq_outcomes_0)/(1-r_sq_outcomes_0)
    
    gamma_0 = matrix(optimized_0$gamma, ncol=1)
    gamma_1 = matrix(optimized$gamma, ncol=1)
    
    r_sq_treatments_0 = 1 - (pi^2/3)/(t(gamma_0) %*% x_cov %*% gamma_0+pi^2/3)
    r_sq_treatments_1 = 1 - (pi^2/3)/(t(gamma_1) %*% x_cov %*% gamma_1 + alpha^2/4 + pi^2/3)
    
    outputs_df[(i-1)*q+j,5] = (r_sq_treatments_1 - r_sq_treatments_0)/(1-r_sq_treatments_0)
  }
}



