
## devtools::install_github("timmens/sawr")
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
nIter <- nT * nN

s_est_mean <- numeric(nIter)
s_est_sd   <- numeric(nIter)
mise_sd    <- numeric(nIter)
mise_mean  <- numeric(nIter)
s_0        <- numeric(nIter)

rng_number <- 123
set.seed(rng_number)

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

for (t in 1:nT) {
  
  t_tmp <- T[t]
  true_beta <- rep(2, t_tmp)
  for (n in 1:nN) {
    n_tmp <- N[n]
    cat(sprintf("dgp = %d; n = %d; t = %d\n", 7, n_tmp, t_tmp))
    
    tmp_result_matrix <- foreach::foreach(r = 1:nSim, .combine = cbind) %dopar% {
      
      data <- DGP(t_tmp, n_tmp, NULL, 7)
      Z    <- make_instruments(data$X) 
      
      results <- sawr::saw_fun(data$Y ~ Z)
      estimated_taus <- results$tausList[[1]]
      
      s_est_mean_tmp <- sum(!is.na(estimated_taus))
      
      mise_tmp    <- mean((true_beta - results$betaMat)^2)
      
      tmp_results <- c(s_est_mean_tmp, 
                       mise_tmp)
      tmp_results
    }
    
    index             <- (t - 1) * nN + n 
    s_est_mean[index] <- mean(tmp_result_matrix[1, ], na.rm = TRUE)
    s_est_sd[index]   <- sd(tmp_result_matrix[1, ], na.rm = TRUE)
    
    mise_vec          <- tmp_result_matrix[2, ]
    mise_mean[index]  <- mean(mise_vec[!is.infinite(mise_vec)], na.rm = TRUE)
    mise_sd[index]    <- sd(mise_vec[!is.infinite(mise_vec)],   na.rm = TRUE)
    
    s_0[index]        <- sum(is.na(tmp_result_matrix[1, ])) / nSim
    
    cat(sprintf("%2.2f percent done\n", index / nIter * 100))
    
  }
}
parallel::stopCluster(cl)

params    <- list(N = N, T = T)
result_df <- expand.grid(params)[c(2, 1)]

additional_info      <- character(nrow(result_df))
additional_info[1:2] <- c(paste0("nsim = ", nSim), 
                          paste0("seed = ", rng_number))

result_df$s_est_mean <- s_est_mean
result_df$s_est_sd   <- s_est_sd
result_df$mise_mean  <- mise_mean
result_df$mise_sd    <- mise_sd
result_df$s_0        <- s_0
result_df$additional_info <- additional_info

filename <- 
  paste0(paste0("r_simulation_results/rsimulaton-dgp7-", 
                gsub(" ", "-", gsub(":", "-", as.character(Sys.time())))), ".csv")

write.csv(result_df, file = filename)
