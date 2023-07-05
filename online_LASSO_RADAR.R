online_LASSO_full_RADAR <- function (B, subset_index, tempdatadir, p, original_eta, lambda_seq, intercept = FALSE, T_list, R_1, out_index){  
  
  maxit <- 1
  tol <- 1e-6
  k <- 1

  sub_length <- length(subset_index)

  beta_de_mat <- sd_de_mat <-  matrix(rep(NA, sub_length * length(out_index)), nrow = sub_length)
  sigma_mat <- matrix(rep(NA, 1 * length(out_index)), nrow = 1)
  
  beta_lambda_new <- matrix(rep(0, p * length(lambda_seq)), nrow = p)
  
  initial_theta <- matrix(rep(0, p * length(lambda_seq)), nrow = p)
  
  beta_tilde <- matrix(rep(0, p * length(1)), nrow = p)
    
  tilde_theta <- matrix(rep(0, p * length(1)), nrow = p)
  
  sigma_new <- 0

  gamma_new <- matrix(rep(0, (p - 1) * sub_length), nrow = p - 1)
  
  gamma_tilde <- matrix(rep(0, (p - 1) * sub_length), nrow = p - 1)

  tilde_theta_gamma <- gamma_tilde
  
  N_new <- 0
  
  zz_r <- ztx_r <- zty_r <- rep(0, p) # p * 1 vec
  ztX_r <- matrix(rep(0, p * sub_length), nrow = p) 


  lambda_s_seq <- c()
  sigma_vec <- c()

  time_load <- 0
  ptm <- proc.time()
  

  for(b in 1 : B){
    

    load <- proc.time()
    load(paste(tempdatadir, "/", b, ".RData", sep = ""))
    time_load <- time_load + (proc.time() - load)[3]
    
    if(intercept == TRUE) X <- cbind(1, X)
    
    if (b == 1){
      set.seed(101)
      n <- nrow(X)
      if (n > 1){
        sample <- sample.int(n, size = floor(0.80 * n), replace = F)
        index1 <- sample - 1 # adjust to C++ index starting from 0
        index2 <- setdiff(seq(0, n - 1, 1), index1)}
      if (n == 1){ 
        index1 <- 0
        index2 <- 0}
    }
    
    beta_lambda <- beta_lambda_new
    N_new <- length(y)
    # Cpp code from here
    
    
    if (b == (T_list[k] + 1)){  # initialize
      
      mu <- matrix(rep(0, p * length(lambda_seq)), nrow = p)
      beta_lambda <- matrix(rep(tilde_theta, length(lambda_seq)), nrow = p)
      theta <- matrix(rep(0, p *length(lambda_seq)), nrow = p)
      tilde_theta <- matrix(rep(0, p * length(1)), nrow = p)             
      count_int_k <- 1
      R_k <- R_1 / (sqrt(2) ** (k - 1))
      
      mu_gamma <- matrix(rep(0, (p - 1) * length(subset_index)), nrow = p - 1)
      gamma_new <- tilde_theta_gamma
      theta_gamma <- matrix(rep(0, (p - 1) * length(subset_index)), nrow = p - 1)
      tilde_theta_gamma <- matrix(rep(0, (p - 1) * length(subset_index)), nrow = p - 1)     
      R_k <- R_1 / (sqrt(2) ** (k - 1))
      
      if(k == 2){sigma_new <- 0}
      
      k <- k + 1
    }
    
    eta <- (original_eta / sqrt(count_int_k)) * (R_k / R_1)   
    
    estimate <- online_Lasso_full(X, y, beta_lambda, subset_index, N_new,
      index1, index2, gamma_new, zz_r, ztx_r, zty_r, ztX_r, 
      lambda_seq, eta, b, maxit, tol, beta_tilde, gamma_tilde, R_k, mu, theta, tilde_theta, count_int_k,
      mu_gamma, theta_gamma, tilde_theta_gamma, k)
    
    
    count_int_k <- count_int_k + 1
    
    sigma_sub <- estimate$sigma_ols

    sigma_new_square <- (sigma_new^2) * (b - T_list[2] - 1) / (b - T_list[2]) + (sigma_sub^2) / (b - T_list[2])
    sigma_new <- sqrt(sigma_new_square)
    
    beta_de <- estimate$beta_de
    sd_de <- estimate$sd_de  * sigma_new
    lambda_s <- estimate$lambda_s
    
    beta_lambda_new <- estimate$beta_lambda_new
    
    zz_r <- estimate$zz_r
    ztx_r <- estimate$ztx_r
    zty_r <- estimate$zty_r
    ztX_r <- estimate$ztX_r
    gamma_new <- estimate$gamma_new

    lambda_s_seq <- c(lambda_s_seq, lambda_s)
    
    if(b %in% out_index){
        ind <- which(out_index == b)
        beta_de_mat[, ind] <- beta_de  # dim is p * B 
        sd_de_mat[, ind] <- sd_de
        sigma_mat[, ind] <- sigma_new
        
    }
    
    mu <- estimate$mu
    theta <- estimate$theta
    tilde_theta <- estimate$tilde_theta
    mu_gamma <- estimate$mu_gamma
    theta_gamma <- estimate$theta_gamma
    tilde_theta_gamma <- estimate$tilde_theta_gamma
    

  }
  
  out <- list()
  out$beta_lambda <- beta_lambda_new
  out$beta_de_mat <- beta_de_mat  # debiased estimator, each row corresponds to a step
  out$sd_de_mat <- sd_de_mat
  
  out$sigma_mat <- sigma_mat

  out$lambda_s_seq <- lambda_s_seq
  
  time_total <- (proc.time() - ptm)[3]
  time_run <- time_total - time_load 
  
  out$time <- time_total
  out$run <- time_run
  
  return(out)
}


