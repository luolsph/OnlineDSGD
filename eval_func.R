# offline debiased lasso with hdi package

delasso.all <- function(X, Y, family){
  fit <- lasso.proj(x = X, y = Y, family = family)
  CI <- confint(fit, level = 0.95)
  CI.lower <- CI[, 1]
  CI.upper <- CI[, 2]

  betahat <- (CI.lower + CI.upper) / 2
  sd <- (CI.upper - betahat) / 1.96
  phi <- fit$sigmahat
  cbind(betahat = betahat, sd = sd, phi = phi)
}


##########################################
## Evaluation functions
##########################################
eval.func.offline <- function(alpha0, alphahat, sd){
  
  alpha_b <- abs(alphahat - alpha0) # p * B matrix

  pvalue <- 2 * pnorm(-abs((alphahat - alpha0) / sd))  # point-wise division
  
  covprob <- ifelse(pvalue >= 0.05, 1, 0)

  CI_length <- sd * 1.96 * 2

  cbind(a.bias = mean(alpha_b), a.se = mean(sd), a.cp = mean(covprob), a.CI = mean(CI_length))
}

eval.func <- function(alpha0, alphahat, sd){
  
  alpha_b <- abs(alphahat - alpha0) # p * B matrix
 # a.bias <- apply(alpha_b, 2, mean) # average over p covariates
 # a.se <- apply(sd, 2, mean)
  CI_length <- sd * 1.96 * 2 
 # a.CI <- apply(CI_length, 2, mean)

  pvalue <- 2 * pnorm(-abs((alphahat - alpha0) / sd))  # point-wise division
  
  covprob <- ifelse(pvalue >= 0.05, 1, 0)
#  a.cp <- apply(covprob, 2, mean)

  cbind(a.bias = alpha_b, a.se = sd, a.cp = covprob, a.CI = CI_length)
}

eval.func.small <- function(alpha0, alphahat, sd){
  
  alpha_b <- abs(alphahat - alpha0) # p * B matrix
  a.bias <- apply(alpha_b, 2, mean) # average over p covariates
  a.se <- apply(sd, 2, mean)
  CI_length <- sd * 1.96 * 2 
  a.CI <- apply(CI_length, 2, mean)
  
  pvalue <- 2 * pnorm(-abs((alphahat - alpha0) / sd))  # point-wise division
  
  covprob <- ifelse(pvalue >= 0.05, 1, 0)
  a.cp <- apply(covprob, 2, mean)
  
  cbind(a.bias, a.se, a.cp, a.CI)
}
