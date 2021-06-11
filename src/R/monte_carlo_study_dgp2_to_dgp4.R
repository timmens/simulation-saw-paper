# dgp 2 to 4

library("yaml")
library("sawr")
library("readr")
library("foreach")
library("doParallel")

source("src/R/dgp.R")  # exports function DGP


# read config parameters from file
config <- yaml.load_file(file.path("src", "config.yaml"))
test_run <- config$test
sample_sizes <- config$sample_sizes
time_periods <- config$time_periods
n_sims <- config$n_sims
if (n_sims != 500) warning("Check n_sims.")

jumps <- c(1, 2, 3)  # S
dgps <- c(2, 3, 4)


n_iter <- length(sample_sizes) * length(time_periods) * length(jumps) * length(dgps)

s_est_mean  <- numeric(n_iter)
s_est_sd    <- numeric(n_iter)
mdcj_mean   <- numeric(n_iter)
mdcj_sd     <- numeric(n_iter)
mise_mean   <- numeric(n_iter)
mise_sd     <- numeric(n_iter)
hd_mean     <- numeric(n_iter)
hd_sd       <- numeric(n_iter)
s_0         <- numeric(n_iter)
taed_mean   <- numeric(n_iter)  # time average euclidian distance (taed)
taed_sd     <- numeric(n_iter)

# simulation
seed <- 123
set.seed(seed)
starting_time <- Sys.time()

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

for (dgp in dgps) {
  
  for (t in time_periods) {
    
    for (s in jumps) { 
      
      for (n in sample_sizes) { 
        
        true_beta_tau <- make_beta(t, s, n)
        
        cat(sprintf("dgp = %d; n = %d; s = %d; t = %d\n", dgp, n, s, t))
        
        tmp_result_matrix <- foreach::foreach(r = 1:n_sims, .combine = cbind) %dopar% {
          
          # create data
          data  <- DGP(t, n, true_beta_tau$beta, dgp)
          
          # apply method
          if (dgp == 2) {
            # in endogeneous case use cross validation
            results <- sawr::fit_saw_iter(y=data$Y, X=data$X, Z=data$Z, max_iter=20)
          } else {
            results <- sawr::fit_saw(y=data$Y, X=data$X, Z=data$Z)
          }
          
          ## evaluate metrics
          estimated_taus <- results$jump_locations[[1]]
          s_est_mean_tmp <- sum(!is.na(estimated_taus))
          mdcj_tmp       <- MDCJ(true_beta_tau$tau, estimated_taus, s)
          mise_tmp       <- mean((true_beta_tau$beta - results$beta_matrix)^2)
          hd_mean_tmp    <- dist_hausdorff(true_beta_tau$tau, estimated_taus) / t
          
          # time-average of euclidian distance
          gamma_true <- beta_to_gamma(true_beta_tau$beta)
          taed_tmp <- dist_euclidean_time_average(results$gamma_hat, gamma_true)
          
          inner_loop_results    <- c(
            s_est_mean_tmp,
            ifelse(is.na(mdcj_tmp), NA_real_, mdcj_tmp),
            mise_tmp,
            ifelse(is.na(hd_mean_tmp), NA_real_, hd_mean_tmp),
            taed_tmp
          )
          inner_loop_results
        }
        
        .dgp <- which(dgp == dgps)
        .t <- which(t == time_periods)
        .n <- which(n == sample_sizes)
        .s <- which(s == jumps)
        index <- (.dgp - 1) * length(time_periods) * length(jumps) * length(sample_sizes) + 
                           (.t - 1) * length(jumps) * length(sample_sizes) + 
                                  (.s - 1) * length(sample_sizes) + 
                                             .n 
        
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
        
        s_0[index]        <- sum(is.na(tmp_result_matrix[2, ])) / n_sims
        
        taed_mean[index] <- mean(tmp_result_matrix[5, ], na.rm = TRUE)
        taed_sd[index] <- sd(tmp_result_matrix[5, ], na.rm = TRUE)
        
        cat(sprintf("%2.2f percent done\n", index / n_iter * 100))
        cat(s_est_mean[index])
      }
    }
  }
}
parallel::stopCluster(cl)

# save results to data frame then to csv file
params                    <- list(N = sample_sizes, S = jumps, T = time_periods, dgp = dgps)
result_df                 <- expand.grid(params)[c(4, 3, 2, 1)]

result_df$s_est_mean      <- s_est_mean
result_df$s_est_sd        <- s_est_sd
result_df$mise_mean       <- mise_mean
result_df$mise_sd         <- mise_sd
result_df$mdcj_mean       <- mdcj_mean
result_df$mdcj_sd         <- mdcj_sd
result_df$hd_mean         <- hd_mean
result_df$hd_sd           <- hd_sd
result_df$s_0             <- s_0
result_df$taed_mean <- taed_mean
result_df$taed_sd <- taed_sd

# write to file
date_time_str = substr(gsub(" ", "-", gsub(":", "-", as.character(Sys.time()))), 1, 16)
output_dir = file.path("bld", "R")
if (test_run) {
  file_name = file.path(output_dir, paste0("simulation_dgp2-to-dgp4_", date_time_str, ".csv"))
} else{
  file_name = file.path(output_dir, "simulation_dgp2-to-dgp4.csv")
}
write_csv(result_df, file_name)
