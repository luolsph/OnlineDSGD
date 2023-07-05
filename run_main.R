
library(MASS)
library(Rcpp)
library(Matrix)
library(hdi)


source("datagenerator.R")
source("eval_func.R")

source("online_LASSO_ASGD.R")
sourceCpp("online_Lasso_ASGD.cpp")

source("online_LASSO_RADAR.R")
sourceCpp("online_Lasso_RADAR.cpp")

source("offline_LASSO_RADAR.R")
sourceCpp("offline_Lasso_RADAR.cpp")



#number of simulations
nsim <- 200
NB <- 100
B <- 100
p <- 400


C <- 1 # partition big X into C columns (fast generating X)


out_index <- c(40, 70, 100)*B/NB
  

S <- length(out_index)
# For RADAR algorithm
T_list <- c(0, ceiling(2^c(0:10) * log(p) *(B/NB)))*3

set.seed(1)
k <- 6
R <- k / 2

# k <- 40#6
# R <- 10#k / 2

true <- sample(1 : p, size = k, replace = FALSE)
beta0 <- rep(0, p)

beta01 <- 1
beta02 <- 0.01

true_beta <- rep(c(1, 0.01), each = k / 2)

beta0[true] <- true_beta

phi <- 1

# corst_x <- "ind"
corst_x <- "AR-1"
rho_x <- 0.5

family <- "gaussian"

# three sets based on signal strength
A1 <- which(beta0 == 0)
A2 <- which(beta0 == 0.01)
A3 <- which(beta0 == 1)

subset_index <- c(1 : p) #c(A1, A2, A3)
# beta0_sub <- beta0[subset_index]

intercept <- FALSE

nb <- round(NB / B)

# tuning parameters
eta1 <- 0.0025 #0.0025 #0.005  # tuning parameter for ASGD
eta2 <- 0.0025 # tuning parameter for RADAR
lambda_seq <- seq(0.30, 0.45, 0.05)

dir.create(paste(p))
dir.create('data')

outputfilename <- paste(family, "_", "N", NB, "B", B, "p", p, "s0", k, corst_x, eta1, eta2, sep = "")

outputfilename_store <- paste(p, "/", family, "_", "N", NB, "B", B, "p", p, "s0", k, corst_x, eta1, eta2, sep = "")

RADAR_on <- TRUE
ASGD_on <- TRUE

RADAR_off <-  TRUE
hdi <- TRUE
lambda_c <- 3

for(s in c(1 : nsim)){
    print(s)

    tempdatadir <- paste('data', "/", "Temp_", outputfilename, "_", s, sep = "")

    datagenerator(nb, p, B, C, tempdatadir, beta0 = beta0, phi = phi, intercept = intercept, 
      corst_x = corst_x, rho_x = rho_x, categorical = FALSE, seed = s)
    
    
    if(hdi == TRUE){
    time.all.read <- system.time(data.A.full <- fulldata(B = B, tempdatadir = tempdatadir))[3]
    
    time.run.delasso <- system.time(result.A.delasso <- delasso.all(X = data.A.full$X, Y = data.A.full$y,
                                                                    family = "gaussian"))[3]
    print("debiased lasso on full data: done!")}
    
    
    if(ASGD_on == TRUE){
    result.B.lasso.ASGD <- online_LASSO_full_ASGD(B = B, subset_index = subset_index, tempdatadir = tempdatadir, 
      p = p, original_eta = eta1, lambda_seq = lambda_seq, 
      intercept = intercept, out_index = out_index)
    print("online lasso ASGD: done!")}

    if(RADAR_on == TRUE){
    result.B.lasso.RADAR <- online_LASSO_full_RADAR(B = B, subset_index = subset_index, tempdatadir = tempdatadir, 
      p = p, original_eta = eta2, lambda_seq = lambda_seq, 
      intercept = intercept, T_list = T_list, R_1 = R, out_index = out_index)
    print("online lasso RADAR: done!")}

    if(RADAR_off == TRUE){
      result.B.lasso.RADAR_off <- offline_LASSO_full_RADAR(B = B, subset_index = subset_index, tempdatadir = tempdatadir, 
        p = p, original_eta = eta2, lambda_seq = lambda_seq, 
        intercept = intercept, T_list = T_list, R_1 = R, out_index = out_index, lambda_c = lambda_c, phi = phi)
      print("offline lasso RADAR: done!")}

    unlink(tempdatadir, recursive = TRUE)
     
    if(hdi == TRUE){
    out_delasso <- c(
      eval.func.offline(beta0[A1], result.A.delasso[A1, 1], result.A.delasso[A1, 2]), 
      eval.func.offline(beta0[A2], result.A.delasso[A2, 1], result.A.delasso[A2, 2]), 
      eval.func.offline(beta0[A3], result.A.delasso[A3, 1], result.A.delasso[A3, 2]),
      time.run.delasso + time.all.read, time.run.delasso) 
    out_ese_delasso <- result.A.delasso[, 1]
    write.table(as.matrix(t(out_delasso)), file = paste(outputfilename_store, "delasso.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE) 
    write.table(as.matrix(t(out_ese_delasso)), file = paste(outputfilename_store, "delasso.ese.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE) 
    }
    
    
    
    # evaluation at every time point
    if(ASGD_on == TRUE){
    out_online_ASGD <- c(
     eval.func.small(beta0[A1], result.B.lasso.ASGD$beta_de_mat[A1, ], result.B.lasso.ASGD$sd_de_mat[A1, ]), 
     eval.func.small(beta0[A2], result.B.lasso.ASGD$beta_de_mat[A2, ], result.B.lasso.ASGD$sd_de_mat[A2, ]), 
     eval.func.small(beta0[A3], result.B.lasso.ASGD$beta_de_mat[A3, ], result.B.lasso.ASGD$sd_de_mat[A3, ]),
     result.B.lasso.ASGD$time, result.B.lasso.ASGD$run)
    
    out_ese_ASGD <- c(result.B.lasso.ASGD$beta_de_mat)
    out_tuning <- result.B.lasso.ASGD$lambda_s_seq
    write.table(as.matrix(t(out_online_ASGD)), file = paste(outputfilename_store, "online.asgd.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE)  
    
    write.table(as.matrix(t(out_ese_ASGD)), file = paste(outputfilename_store, "online.ese.asgd.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE)
    write.table(as.matrix(t(out_tuning)), file = paste(outputfilename_store, "para.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE)  
    }

    
    if(RADAR_on == TRUE){
    out_online_RADAR <- c(
     eval.func.small(beta0[A1], result.B.lasso.RADAR$beta_de_mat[A1, ], result.B.lasso.RADAR$sd_de_mat[A1, ]), 
     eval.func.small(beta0[A2], result.B.lasso.RADAR$beta_de_mat[A2, ], result.B.lasso.RADAR$sd_de_mat[A2, ]), 
     eval.func.small(beta0[A3], result.B.lasso.RADAR$beta_de_mat[A3, ], result.B.lasso.RADAR$sd_de_mat[A3, ]),
     result.B.lasso.RADAR$time, result.B.lasso.RADAR$run)
    
    write.table(result.B.lasso.RADAR$sigma_mat, file = paste(outputfilename_store, "sigma.online.radar.csv", sep = ""), sep = ",",
                col.name = FALSE, row.names = s,append = TRUE)   
    
    out_ese_RADAR <- c(result.B.lasso.RADAR$beta_de_mat) # p * B # this can be used to produce histogram
    
    
    write.table(as.matrix(t(out_online_RADAR)), file = paste(outputfilename_store, "online.radar.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE) 
    write.table(as.matrix(t(out_ese_RADAR)), file = paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",",
                col.names = FALSE, row.names = s, append = TRUE)
    }
    
    
    if(RADAR_off == TRUE){
      out_offline_RADAR <- c(
        eval.func.small(beta0[A1], result.B.lasso.RADAR_off$beta_de_mat[A1, ], result.B.lasso.RADAR_off$sd_de_mat[A1, ]), 
        eval.func.small(beta0[A2], result.B.lasso.RADAR_off$beta_de_mat[A2, ], result.B.lasso.RADAR_off$sd_de_mat[A2, ]), 
        eval.func.small(beta0[A3], result.B.lasso.RADAR_off$beta_de_mat[A3, ], result.B.lasso.RADAR_off$sd_de_mat[A3, ]),
        result.B.lasso.RADAR_off$time, result.B.lasso.RADAR_off$run)
      
      write.table(result.B.lasso.RADAR_off$sigma_mat, file = paste(outputfilename_store, "sigma.offline.radar.csv", sep = ""), sep = ",",
                  col.name = FALSE, row.names = s,append = TRUE)   
      
      out_ese_RADAR <- c(result.B.lasso.RADAR_off$beta_de_mat) # p * B # this can be used to produce histogram
      
      
      write.table(as.matrix(t(out_offline_RADAR)), file = paste(outputfilename_store, "offline.radar.csv", sep = ""), sep = ",",
                  col.names = FALSE, row.names = s, append = TRUE) 
      write.table(as.matrix(t(out_ese_RADAR)), file = paste(outputfilename_store, "offline.ese.radar.csv", sep = ""), sep = ",",
                  col.names = FALSE, row.names = s, append = TRUE)
    }
    
  

}

if(hdi == TRUE){
### results of delasso ###
out_delasso <- read.table(paste(outputfilename_store, "delasso.csv", sep = ""), sep = ",")[, -1]
colnames(out_delasso) <- c("delasso.bias.A1", "delasso.sd.A1", "delasso.cp.A1", "delasso.length.A1", 
                           "delasso.bias.A2", "delasso.sd.A2", "delasso.cp.A2", "delasso.length.A2", 
                           "delasso.bias.A3", "delasso.sd.A3", "delasso.cp.A3", "delasso.length.A3", 
                           "c.time.delasso", "r.time.delasso")
result_delasso <- apply(out_delasso, 2, mean)
print(result_delasso)

# bias of delasso
ese_delasso <- read.table(paste(outputfilename_store, "delasso.ese.csv", sep = ""), sep = ",")[, -1]
delasso_est <- apply(ese_delasso, 2, mean)
c(mean(abs(delasso_est[A1] - 0)), mean(abs(delasso_est[A2] - 0.01)), mean(abs(delasso_est[A3] - 1)))


# ESE of delasso
ese_delasso <- read.table(paste(outputfilename_store, "delasso.ese.csv", sep = ""), sep = ",")[, -1]
delasso_sd <- apply(ese_delasso, 2, sd)
c(mean(delasso_sd[A1]), mean(delasso_sd[A2]), mean(delasso_sd[A3]))}



### results of sgd debias ###
if(ASGD_on == TRUE){
  out_online <- read.table(paste(outputfilename_store, "online.asgd.csv", sep = ""), sep = ",")[1:nsim, -1]
  out_online <- apply(out_online, 2, mean)
  for (c in 1:3) {
    #c <- 2
    online_colname <- c(paste("online.bias.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("online.sd.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("online.cp.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("online.length.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = "")
    )
    out_online_A <- out_online[(4 * S * (c - 1) + 1) : (4 * S * c)]
    names(out_online_A) <- online_colname
    print(out_online_A)
    
    start = length(out_index)*2 + 1
    end = length(out_index)*3 
    #print(out_online_A[start:end])
    result = c(out_online_A[start:end], NB/B, eta2, lambda_seq, T_list)
    dim(result) = c(1,length(result))
    
  }
  
  
  
  # estimate of radar
  ese_online <- read.table(paste(outputfilename_store, "online.ese.asgd.csv", sep = ""), sep = ",")[1:nsim, -1]
  online_est <- apply(ese_online, 2, mean)
  ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    ese_online_seq[b, ] <- online_est[((b - 1) * 3 + 1) : (b * 3)]
  }
  colnames(ese_online_seq) <- c("A1", "A2", "A3")
  print(ese_online_seq)
  
  
  
  # bias of radar
  ese_online <- read.table(paste(outputfilename_store, "online.ese.asgd.csv", sep = ""), sep = ",")[1:nsim, -1]
  online_est <- apply(ese_online, 2, mean)
  ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    online_est_b <- online_est[((b - 1) * p + 1) : (b * p)]
    ese_online_seq[b, ] <- c(mean(abs(online_est_b[A1] - 0)), mean(abs(online_est_b[A2] - 0.01)), 
                             mean(abs(online_est_b[A3] - 1)))
  }
  colnames(ese_online_seq) <- c("A1", "A2", "A3")
  print(ese_online_seq)
  
  
  # ESE of RADAR
  ese_online <- read.table(paste(outputfilename_store, "online.ese.asgd.csv", sep = ""), sep = ",")[1:nsim, -1]
  online_sd <- apply(ese_online, 2, sd)
  ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    online_sd_b <- online_sd[((b - 1) * p + 1) : (b * p)]
    ese_online_seq[b, ] <- c(mean(online_sd_b[A1]), mean(online_sd_b[A2]), mean(online_sd_b[A3]))
  }
  colnames(ese_online_seq) <- c("A1", "A2", "A3")
  print(ese_online_seq)
  
  time <- out_online[(length(out_online) -1) : length(out_online)]
  names(time) <- c("c.time.online","r.time.online")
  print(time)}



### results of RADAR debias ###
if(RADAR_on == TRUE){
out_online <- read.table(paste(outputfilename_store, "online.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
out_online <- apply(out_online, 2, mean)
for (c in 1:3) {
  #c <- 2
  online_colname <- c(paste("online.bias.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                      paste("online.sd.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                      paste("online.cp.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                      paste("online.length.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = "")
  )
  out_online_A <- out_online[(4 * S * (c - 1) + 1) : (4 * S * c)]
  names(out_online_A) <- online_colname
  print(out_online_A)
  
  start = length(out_index)*2 + 1
  end = length(out_index)*3 
  #print(out_online_A[start:end])
  result = c(out_online_A[start:end], NB/B, eta2, lambda_seq, T_list)
  dim(result) = c(1,length(result))
  
  write.table(as.matrix(result), file = paste("pool_summary_minibatch.csv", sep = ""), sep = ",",
              col.names = FALSE, append = TRUE)  
}



# estimate of radar
ese_online <- read.table(paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
#ese_online <- read.table(paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",")[, -1]
online_est <- apply(ese_online, 2, mean)
ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
for (b in 1 : S){
  ese_online_seq[b, ] <- online_est[((b - 1) * 3 + 1) : (b * 3)]
}
colnames(ese_online_seq) <- c("A1", "A2", "A3")
print(ese_online_seq)



# bias of radar
ese_online <- read.table(paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",")[, -1]
online_est <- apply(ese_online, 2, mean)
ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
for (b in 1 : S){
  online_est_b <- online_est[((b - 1) * p + 1) : (b * p)]
  ese_online_seq[b, ] <- c(mean(abs(online_est_b[A1] - 0)), mean(abs(online_est_b[A2] - 0.01)), 
                           mean(abs(online_est_b[A3] - 1)))
}
colnames(ese_online_seq) <- c("A1", "A2", "A3")
print(ese_online_seq)


# ESE of RADAR
ese_online <- read.table(paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",")[, -1]
online_sd <- apply(ese_online, 2, sd)
ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
for (b in 1 : S){
  online_sd_b <- online_sd[((b - 1) * p + 1) : (b * p)]
  ese_online_seq[b, ] <- c(mean(online_sd_b[A1]), mean(online_sd_b[A2]), mean(online_sd_b[A3]))
}
colnames(ese_online_seq) <- c("A1", "A2", "A3")
print(ese_online_seq)

time <- out_online[(length(out_online) -1) : length(out_online)]
names(time) <- c("c.time.online","r.time.online")
print(time)}



### results of RADAR offline debias ###
if(RADAR_off == TRUE){
  out_offline <- read.table(paste(outputfilename_store, "offline.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
  out_offline <- apply(out_offline, 2, mean)
  for (c in 1:3) {
    #c <- 2
    offline_colname <- c(paste("offline.bias.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("offline.sd.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("offline.cp.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = ""),
                        paste("offline.length.A", sprintf("%02d", c), ".", sprintf("%02d", 1:S), sep = "")
    )
    out_offline_A <- out_offline[(4 * S * (c - 1) + 1) : (4 * S * c)]
    names(out_offline_A) <- offline_colname
    print(out_offline_A)
    
    start = length(out_index)*2 + 1
    end = length(out_index)*3 
    #print(out_online_A[start:end])
    result = c(out_offline_A[start:end], NB/B, eta2, lambda_seq, T_list)
    dim(result) = c(1,length(result))
    
  }
  
  
  
  # estimate of radar
  ese_offline <- read.table(paste(outputfilename_store, "offline.ese.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
  #ese_online <- read.table(paste(outputfilename_store, "online.ese.radar.csv", sep = ""), sep = ",")[, -1]
  offline_est <- apply(ese_offline, 2, mean)
  ese_offline_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    ese_offline_seq[b, ] <- offline_est[((b - 1) * 3 + 1) : (b * 3)]
  }
  colnames(ese_offline_seq) <- c("A1", "A2", "A3")
  print(ese_offline_seq)
  
  
  
  # bias of radar
  ese_offline <- read.table(paste(outputfilename_store, "offline.ese.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
  offline_est <- apply(ese_offline, 2, mean)
  ese_online_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    offline_est_b <- offline_est[((b - 1) * p + 1) : (b * p)]
    ese_offline_seq[b, ] <- c(mean(abs(offline_est_b[A1] - 0)), mean(abs(offline_est_b[A2] - 0.01)), 
                             mean(abs(offline_est_b[A3] - 1)))
  }
  colnames(ese_offline_seq) <- c("A1", "A2", "A3")
  print(ese_offline_seq)
  
  
  # ESE of RADAR
  ese_offline <- read.table(paste(outputfilename_store, "offline.ese.radar.csv", sep = ""), sep = ",")[1:nsim, -1]
  offline_sd <- apply(ese_offline, 2, sd)
  ese_offline_seq <- matrix(rep(NA, S * 3), nrow = S)
  for (b in 1 : S){
    offline_sd_b <- offline_sd[((b - 1) * p + 1) : (b * p)]
    ese_offline_seq[b, ] <- c(mean(offline_sd_b[A1]), mean(offline_sd_b[A2]), mean(offline_sd_b[A3]))
  }
  colnames(ese_offline_seq) <- c("A1", "A2", "A3")
  print(ese_offline_seq)
  
  time <- out_offline[(length(out_online) -1) : length(out_offline)]
  names(time) <- c("c.time.offline","r.time.offline")
  print(time)}



if(ASGD_on == TRUE){
tuning <- read.table(paste(outputfilename_store, "para.csv", sep = ""), sep = ",")[, -1]
#table(tuning)
apply(tuning, 2, table)}




