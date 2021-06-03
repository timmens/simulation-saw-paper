# dgp 6 (no-jumps)

library("yaml")
library("sawr")
library("readr")
library("foreach")
library("doParallel")

source("src/R/dgp.R")  # exports function dgp1

# read config parameters from file
config <- yaml.load_file(file.path("src", "config.yaml"))
test_run <- config$test
sample_sizes <- config$sample_sizes
time_periods <- config$time_periods
n_sims <- config$n_sims
if (n_sims != 500) warning("Check n_sims.")

n_iter <- length(sample_sizes) * length(time_periods)

s_est_mean <- numeric(n_iter)
s_est_sd   <- numeric(n_iter)
mise_sd    <- numeric(n_iter)
mise_mean  <- numeric(n_iter)
s_0        <- numeric(n_iter)
taed_mean  <- numeric(n_iter)  # time average euclidian distance (taed)
taed_sd    <- numeric(n_iter)

seed <- 123
set.seed(seed)
starting_time <- Sys.time()

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

for (t in time_periods) {
  
  true_beta <- rep(1, t)
  
  for (n in sample_sizes) {
    
    cat(sprintf("dgp = %d; n = %d; t = %d\n", 6, n, t))
    
    foreach_result_matrix <- foreach::foreach(r = 1:n_sims, .combine = cbind) %dopar% {
      
      data <- dgp6(t, n, beta = NULL)
      
      results <- sawr::fit_saw(y=data$Y, X=data$X)
      estimated_taus <- results$jump_locations[[1]]
      
      s_est_tmp <- sum(!is.na(estimated_taus))
      
      mise_tmp    <- mean((true_beta - results$beta_matrix)^2)
      
      # time-average of euclidian distance
      gamma_true <- beta_to_gamma(true_beta)
      taed_tmp <- dist_euclidian_time_average(results$gamma_hat, gamma_true)
      
      inner_loop_results <- c(s_est_tmp, mise_tmp, taed_tmp)
      inner_loop_results
    }
    
    .t <- which(t == time_periods)
    .n <- which(n == sample_sizes)
    index <- (.t - 1) * length(sample_sizes) + .n 
    
    s_est_mean[index] <- mean(foreach_result_matrix[1, ], na.rm = TRUE)
    s_est_sd[index]   <- sd(foreach_result_matrix[1, ], na.rm = TRUE)
    
    mise_vec          <- foreach_result_matrix[2, ]
    mise_mean[index]  <- mean(mise_vec[!is.infinite(mise_vec)], na.rm = TRUE)
    mise_sd[index]    <- sd(mise_vec[!is.infinite(mise_vec)],   na.rm = TRUE)
    
    s_0[index]        <- sum(is.na(foreach_result_matrix[1, ])) / n_sims
    
    taed_mean[index] <- mean(foreach_result_matrix[3, ], na.rm = TRUE)
    taed_sd[index] <- sd(foreach_result_matrix[3, ], na.rm = TRUE)
    
    cat(sprintf("%2.2f percent done\n", index / n_iter * 100))
    
  }
}
parallel::stopCluster(cl)

# save results to data frame then to csv file
params    <- list(N = sample_sizes, T = time_periods)
result_df <- expand.grid(params)[c(2, 1)]

result_df$s_est_mean <- s_est_mean
result_df$s_est_sd   <- s_est_sd
result_df$mise_mean  <- mise_mean
result_df$mise_sd    <- mise_sd
result_df$s_0        <- s_0
result_df$taed_mean <- taed_mean
result_df$taed_sd <- taed_sd

# write to file
date_time_str = substr(gsub(" ", "-", gsub(":", "-", as.character(Sys.time()))), 1, 16)
output_dir = file.path("bld", "R")
if (test_run) {
  file_name = file.path(output_dir, paste0("simulation_dgp6_", date_time_str, ".csv"))
} else{
  file_name = file.path(output_dir, "simulation_dgp6.csv")
}
write_csv(result_df, file_name)
