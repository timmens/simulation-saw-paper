
## devtools::install_github("https://github.com/timmens/sawr", force = TRUE)
library("sawr")
library("foreach")
library("doParallel")
## Reading in all dgp definitions
source("r_simulation_codes/dgp.R")

S             <- c(1, 2, 3)
N             <- c(30, 60, 120, 300)
T             <- 2 ^ c(5, 6, 7) + 1
nSim          <- 1000
dgp_seq       <- c(2:6)
nDGP          <- length(dgp_seq)
nN            <- length(N)
nT            <- length(T)
nS            <- length(S)
##
nIter         <- nDGP * nT * nS * nN
s_est_mean    <- numeric(nIter)
s_est_sd      <- numeric(nIter)
mdcj_mean     <- numeric(nIter)
mdcj_sd       <- numeric(nIter)
mise_mean     <- numeric(nIter)
mise_sd       <- numeric(nIter)
hd_mean       <- numeric(nIter)
hd_sd         <- numeric(nIter)
s_0           <- numeric(nIter)
## 
rng_number    <- 123
set.seed(rng_number)

cl            <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)
starting_time <- Sys.time()
##
for (dgp in dgp_seq) {
  ##
  for (t in 1:nT) {
    ##
    t_tmp <- T[t]
    for (s in 1:nS) { 
      ## 
      s_tmp          <-  S[s]
      true_beta_tau  <-  make_beta(t_tmp, s_tmp)
      ##
      for (n in 1:nN) { 
        ##
        n_tmp <- N[n]
        ##
        cat(sprintf("dgp = %d; n = %d; s = %d; t = %d\n", dgp, n_tmp, s_tmp, t_tmp))
        ##
        tmp_result_matrix <- foreach::foreach(r = 1:nSim, .combine = cbind) %dopar% {
        #tmp_result_matrix <- c()
        #for (r in 1:nSim) {
          
          data    <- DGP(t_tmp, n_tmp, true_beta_tau$beta, dgp)
          
          ## Estimation
          results <- sawr::saw_fun(data$Y ~ data$X, dot=FALSE)
          ##
          estimated_taus <- results$tausList[[1]]
          s_est_mean_tmp <- sum(!is.na(estimated_taus))
          mdcj_tmp       <- MDCJ(true_beta_tau$tau, estimated_taus, s_tmp)
          mise_tmp       <- mean((true_beta_tau$beta - results$betaMat)^2)
          hd_mean_tmp    <- dist_hausdorff(true_beta_tau$tau, estimated_taus)
          tmp_results    <- c(s_est_mean_tmp,
                              ifelse(is.na(mdcj_tmp), NA_real_, mdcj_tmp),
                              mise_tmp,
                              ifelse(is.na(hd_mean_tmp), NA_real_, hd_mean_tmp))

          tmp_results
          #tmp_result_matrix <- cbind(tmp_result_matrix, tmp_results)
        }
        index             <- (dgp - 2) *  nT   * nS   * nN + 
                                       (t - 1) * nS   * nN + 
                                              (s - 1) * nN + 
                                                         n 
        
        s_est_mean[index] <- mean(tmp_result_matrix[1, ], na.rm = TRUE)
        s_est_sd[index]   <- sd(tmp_result_matrix[1, ], na.rm = TRUE)
        mdcj_vec          <- tmp_result_matrix[2, ]
        mdcj_mean[index]  <- mean(mdcj_vec[!is.infinite(mdcj_vec)], na.rm = TRUE)
        mdcj_sd[index]    <- sd(mdcj_vec[!is.infinite(mdcj_vec)], na.rm = TRUE)
        mise_vec          <- tmp_result_matrix[3, ]
        mise_mean[index]  <- mean(mise_vec[!is.infinite(mise_vec)], na.rm = TRUE)
        mise_sd[index]    <- sd(mise_vec[!is.infinite(mise_vec)],   na.rm = TRUE)
        hd_vec            <- tmp_result_matrix[4, ]
        hd_mean[index]    <- mean(hd_vec[!is.infinite(hd_vec)], na.rm = TRUE)
        hd_sd[index]      <- sd(hd_vec[!is.infinite(hd_vec)], na.rm = TRUE)
        s_0[index]        <- sum(is.na(tmp_result_matrix[2, ])) / nSim
        ##
        cat(sprintf("%2.2f percent done\n", index / nIter * 100))
      }
    }
  }
}
##
parallel::stopCluster(cl)
##
params                    <- list(N = N, S = S, T = T, DGP = dgp_seq)
result_df                 <- expand.grid(params)[c(4, 3, 2, 1)]
additional_info           <- character(nrow(result_df))
additional_info[1:3]      <- c(paste0("nsim = ", nSim), 
                               paste0("seed = ", rng_number), 
                               paste0("ellapsed time = ", Sys.time() - starting_time))
##
result_df$s_est_mean      <- s_est_mean
result_df$s_est_sd        <- s_est_sd
result_df$mise_mean       <- mise_mean
result_df$mise_sd         <- mise_sd
result_df$mdcj_mean       <- mdcj_mean
result_df$mdcj_sd         <- mdcj_sd
result_df$hd_mean         <- hd_mean
result_df$hd_sd           <- hd_sd
result_df$s_0             <- s_0
result_df$additional_info <- paste(additional_info, collapse = '; ')
##
filename <- paste0(paste0("r_simulation_results/rsimulaton-dgp2-to-dgp6s-", 
                          gsub(" ", "-", gsub(":", "-", as.character(Sys.time())))), ".csv")

write.csv(result_df, file = filename)
