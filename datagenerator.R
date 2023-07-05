library(MASS)
##########################################
## data generator
##########################################
sigmagenerator <- function(corst, p, rho = 0){
  if(corst == "ind") { Sigma <- diag(rep(1,p)) }
  if(corst == "cs")  { Sigma <- matrix(rho, p, p); diag(Sigma) <- 1 }
  if(corst == "AR-1") { Sigma <- rho^(abs(outer(1:p, 1:p, "-"))) }
  return(Sigma)
}

datagenerator <- function(n, p, B, C, tempdatadir, beta0, phi, intercept = FALSE, 
  corst_x, rho_x, categorical = FALSE, seed = NULL){
  n          # sample size in each data batch
  p          # dimension of fixed effects
  tempdatadir # directory for simulated data
  B          # number of data batches
  C          # number of blocks in covariance matrix
  beta0      # true coefficients
  phi        # variance of the white noise in the observed process
  intercept  # if TRUE, the beta0[1] is the coefficient of intercept
  corst_x    # c("ind", "cs", "ar1")
  rho_x      # parameter of correlation
  seed       # random seed
  
  if(length(n)!=1 & length(n)!= B){ stop("n must be a number or a vector of length B.") }
  if(!is.null(seed)){ set.seed(seed) }
  seed.list <- sample(1:1e8, B, replace=FALSE)
  
  dir.create(tempdatadir)
  
  X_total <- c()
  for(c in 1 : C){
    d <- p / C # number of covariates per block
    Sigma_x <- sigmagenerator(corst = corst_x, p = d, rho = rho_x)
    set.seed(seed.list[B] + c)
    X_sub <- mvrnorm(n = n * B, mu = rep(0, d), Sigma = Sigma_x)
    X_total <- cbind(X_total, X_sub)
  }
  
  for(b in 1:B){

    start <- n * (b - 1) + 1
    end <- n * b
    X <- X_total[start:end,]
    if (n == 1){
      dim(X) <- c(1,length(X))
    }
    if(categorical == TRUE) X[, 2] <- rep(seq(0.1, 1.0, 0.1), n)
    if(intercept == TRUE) X[, 1] <- 1 
    
    set.seed(seed.list[b])
    epsilon <- rnorm(n, 0, sqrt(phi))
    epsilon <- c(t(epsilon))

    Xbeta <- drop(as.matrix(X) %*% beta0)
    
    y <- Xbeta + epsilon

    if(intercept == TRUE) X <- as.matrix(X[, -1])
    save(y, X, file = paste(tempdatadir, "/", b, ".RData", sep = ""))
    
  }
  
}



fulldata <- function(B, tempdatadir){
  X.full <- c(); y.full <- c(); 
  for(b in 1:B){
    load(paste(tempdatadir, "/", b, ".RData", sep=""))
    X.full <- rbind(X.full, X)
    y.full <- c(y.full, y) 
  }
  list(X = X.full, y = y.full)
}



