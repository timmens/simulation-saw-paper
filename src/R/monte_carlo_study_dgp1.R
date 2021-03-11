
## devtools::install_github("https://github.com/timmens/sawr", force = TRUE)
library("sawr")
library("foreach")
library("doParallel")
## Reading in all dgp definitions
source("r_simulation_codes/dgp.R")

N     <- c(30, 60, 120, 300)
T     <- 2 ^ c(5, 6, 7) + 1
nSim  <- 1000
nN    <- length(N)
nT    <- length(T)

s_est1_mean <- numeric(nN * nT)
s_est1_sd   <- numeric(nN * nT)
s_est2_mean <- numeric(nN * nT)
s_est2_sd   <- numeric(nN * nT)
mise1_mean  <- numeric(nN * nT)
mise1_sd    <- numeric(nN * nT)
mise2_mean  <- numeric(nN * nT)
mise2_sd    <- numeric(nN * nT)
hd1_mean    <- numeric(nN * nT)
hd1_sd      <- numeric(nN * nT)
hd2_mean    <- numeric(nN * nT)
hd2_sd      <- numeric(nN * nT)
hd_mean     <- numeric(nN * nT)
hd_sd       <- numeric(nN * nT)
s_01        <- numeric(nN * nT)
s_02        <- numeric(nN * nT)

rng_number <- 123
set.seed(rng_number)

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

nIter         <- nT * nN
starting_time <- Sys.time()
##
for (t in 1:nT) {
  
  t_tmp <- T[t]
  for (n in 1:nN) {
    n_tmp <- N[n]
    cat(sprintf("dgp = %d; n = %d; t = %d\n", 1, n_tmp, t_tmp))
    
    tmp_result_matrix <- foreach::foreach(r = 1:nSim, .combine = cbind) %dopar% {
    #tmp_result_matrix <- c()
    #for (r in 1:nSim) {
      
      data  <- dgp1(t_tmp, n_tmp, sderr_e=sqrt(0.75))
      beta1 <- data$beta1 # true beta 1
      beta2 <- data$beta2 # true beta 2
      tau1  <- data$tau1  # true tau 1
      tau2  <- data$tau2  # true tau 2
      ## 
      results  <- sawr::saw_fun(data$Y ~ data$X1 + data$X2)
      
      # tmp_check <- cbind(c(tau1,tau2), c(results$tausList[[1]], results$tausList[[2]]))
      # colnames(tmp_check) <- c("True","Estim")
      # rownames(tmp_check) <- c("tau_11","tau_12","tau_21","tau_22","tau_23")
      # tmp_check
      
      tausList <- results$tausList
      
      s_est_tmp_1 <- sum(!is.na(tausList[[1]]))
      s_est_tmp_2 <- sum(!is.na(tausList[[2]]))
      
      mise_tmp_1 <- mean((beta1 - results$betaMat[, 1])^2)
      mise_tmp_2 <- mean((beta2 - results$betaMat[, 2])^2)
      
      ## Hausdorff distance normalized by T
      hd_tmp1    <- dist_hausdorff(tau1, tausList[[1]]) 
      hd_tmp2    <- dist_hausdorff(tau2, tausList[[2]]) 
      hd_tmp     <- dist_hausdorff(c(tau1, tau2), unlist(tausList)) 
      
      tmp_results <- c(s_est_tmp_1,
                       s_est_tmp_2,
                       mise_tmp_1,
                       mise_tmp_2,
                       ifelse(is.na(hd_tmp1), NA_real_, hd_tmp1),
                       ifelse(is.na(hd_tmp2), NA_real_, hd_tmp2),
                       ifelse(is.na(hd_tmp),  NA_real_, hd_tmp))
      tmp_results
      #tmp_result_matrix <- cbind(tmp_result_matrix, tmp_results)
    }
    
    index              <- (t - 1) * nN + n 
    s_est1_mean[index] <- sum(tmp_result_matrix[1, ]) / nSim
    s_est2_mean[index] <- sum(tmp_result_matrix[2, ]) / nSim
    s_est1_sd[index]   <- sd(tmp_result_matrix[1, ])
    s_est2_sd[index]   <- sd(tmp_result_matrix[2, ])
    mise1_mean[index]  <- mean(tmp_result_matrix[3, ], na.rm = TRUE)
    mise1_sd[index]    <- sd(tmp_result_matrix[3, ], na.rm = TRUE)
    mise2_mean[index]  <- mean(tmp_result_matrix[4, ], na.rm = TRUE)
    mise2_sd[index]    <- sd(tmp_result_matrix[4, ], na.rm = TRUE)
    hd1_tmp            <- tmp_result_matrix[5, ]
    hd2_tmp            <- tmp_result_matrix[6, ]
    hd_tmp             <- tmp_result_matrix[7, ]
    hd1_mean[index]    <- mean(hd1_tmp[!is.infinite(hd1_tmp)], na.rm = TRUE)
    hd2_mean[index]    <- mean(hd2_tmp[!is.infinite(hd2_tmp)], na.rm = TRUE)
    hd1_sd[index]      <- sd(hd1_tmp[!is.infinite(hd1_tmp)], na.rm = TRUE)
    hd2_sd[index]      <- sd(hd2_tmp[!is.infinite(hd2_tmp)], na.rm = TRUE)
    hd_mean[index]     <- mean(hd_tmp[!is.infinite(hd_tmp)], na.rm = TRUE)
    hd_sd[index]       <- sd(hd_tmp[!is.infinite(hd_tmp)], na.rm = TRUE)
    s_01[index]        <- sum(is.na(tmp_result_matrix[5, ])) / nSim
    s_02[index]        <- sum(is.na(tmp_result_matrix[6, ])) / nSim
    cat(sprintf("%2.2f percent done\n", index / nIter * 100))
  }
}
parallel::stopCluster(cl)

params                <- list(N = N, T = T)
result_df             <- expand.grid(params)[c(2, 1)]
additional_info       <- character(nrow(result_df))
additional_info[1:3]  <- c(paste0("nsim = ", nSim), 
                          paste0("seed = ", rng_number), 
                          paste0("ellapsed time = ", Sys.time() - starting_time))
## 
result_df$s_est1_mean <- s_est1_mean
result_df$s_est2_mean <- s_est2_mean
result_df$s_est1_sd   <- s_est1_sd
result_df$s_est2_sd   <- s_est2_sd
result_df$mise1_mean  <- mise1_mean
result_df$mise1_sd    <- mise1_sd
result_df$mise2_mean  <- mise2_mean
result_df$mise2_sd    <- mise2_sd
result_df$hd1_mean    <- hd1_mean
result_df$hd2_mean    <- hd2_mean
result_df$hd1_sd      <- hd1_sd
result_df$hd2_sd      <- hd2_sd
result_df$hd_mean     <- hd_mean
result_df$hd_sd       <- hd_sd
result_df$s_01        <- s_01
result_df$s_02        <- s_02
result_df$additional_info <- paste(additional_info, collapse = '; ')

filename <- paste0(paste0("r_simulation_results/rsimulaton-dgp1-", 
                gsub(" ", "-", gsub(":", "-", as.character(Sys.time())))), ".csv")

write.csv(result_df, file = filename)
