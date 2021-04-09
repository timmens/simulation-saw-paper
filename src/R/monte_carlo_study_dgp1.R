# dgp 1 (multiple regressors)

# devtools::install_github("https://github.com/timmens/sawr", force = TRUE)
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

s_est1_mean <- numeric(n_iter)
s_est1_sd   <- numeric(n_iter)
s_est2_mean <- numeric(n_iter)
s_est2_sd   <- numeric(n_iter)
mise1_mean  <- numeric(n_iter)
mise1_sd    <- numeric(n_iter)
mise2_mean  <- numeric(n_iter)
mise2_sd    <- numeric(n_iter)
hd1_mean    <- numeric(n_iter)
hd1_sd      <- numeric(n_iter)
hd2_mean    <- numeric(n_iter)
hd2_sd      <- numeric(n_iter)
hd_mean     <- numeric(n_iter)
hd_sd       <- numeric(n_iter)
s_01        <- numeric(n_iter)
s_02        <- numeric(n_iter)
taed_mean   <- numeric(n_iter)  # time average euclidian distance (taed)
taed_sd     <- numeric(n_iter)

seed <- 123
set.seed(seed)
starting_time <- Sys.time()

cl <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

for (t in time_periods) {
  
  for (n in sample_sizes) {
    
    cat(sprintf("dgp = %d; n = %d; t = %d\n", 1, n, t))
    
    foreach_result_matrix <- foreach::foreach(r = 1:n_sims, .combine = cbind) %dopar% {
      
      # create data
      data  <- dgp1(t, n)
      
      beta1 <- data$beta1 # true beta 1
      beta2 <- data$beta2 # true beta 2
      tau1  <- data$tau1  # true tau 1
      tau2  <- data$tau2  # true tau 2
      
      # apply method
      results  <- sawr::saw_fun(y=data$Y, X=data$X)
      tausList <- results$tausList
      
      ## evaluate metrics
      
      # counting predicted jumps
      s_est_tmp_1 <- sum(!is.na(tausList[[1]]))
      s_est_tmp_2 <- sum(!is.na(tausList[[2]]))
      
      # mean integrated squared error
      mise_tmp_1 <- mean((beta1 - results$betaMat[, 1])^2)
      mise_tmp_2 <- mean((beta2 - results$betaMat[, 2])^2)
      
      # (normalized) Hausdorff metrics
      hd_tmp1    <- dist_hausdorff(tau1, tausList[[1]]) 
      hd_tmp2    <- dist_hausdorff(tau2, tausList[[2]]) 
      hd_tmp     <- dist_hausdorff(c(tau1, tau2), unlist(tausList))
      
      # time-average of euclidian distance
      gamma_true <- beta_to_gamma(list(beta1, beta2))
      taed_tmp <- dist_euclidian_time_average(results$gamma, gamma_true)
      
      # return results from parallelized inner loop
      inner_loop_result <- c(
        s_est_tmp_1,
        s_est_tmp_2,
        mise_tmp_1,
        mise_tmp_2,
        ifelse(is.na(hd_tmp1), NA_real_, hd_tmp1),
        ifelse(is.na(hd_tmp2), NA_real_, hd_tmp2),
        ifelse(is.na(hd_tmp),  NA_real_, hd_tmp),
        taed_tmp
      )
      
      inner_loop_result  # return
    }
    
    .t <- which(t == time_periods)
    .n <- which(n == sample_sizes)
    index <- (.t - 1) * length(sample_sizes) + .n
    
    s_est1_mean[index] <- sum(foreach_result_matrix[1, ]) / n_sims
    s_est2_mean[index] <- sum(foreach_result_matrix[2, ]) / n_sims
    s_est1_sd[index]   <- sd(foreach_result_matrix[1, ])
    s_est2_sd[index]   <- sd(foreach_result_matrix[2, ])
    
    mise1_mean[index]  <- mean(foreach_result_matrix[3, ], na.rm = TRUE)
    mise1_sd[index]    <- sd(foreach_result_matrix[3, ], na.rm = TRUE)
    mise2_mean[index]  <- mean(foreach_result_matrix[4, ], na.rm = TRUE)
    mise2_sd[index]    <- sd(foreach_result_matrix[4, ], na.rm = TRUE)
    
    hd1_tmp            <- foreach_result_matrix[5, ]
    hd2_tmp            <- foreach_result_matrix[6, ]
    hd_tmp             <- foreach_result_matrix[7, ]
    
    hd1_mean[index]    <- mean(hd1_tmp[!is.infinite(hd1_tmp)], na.rm = TRUE)
    hd2_mean[index]    <- mean(hd2_tmp[!is.infinite(hd2_tmp)], na.rm = TRUE)
    hd1_sd[index]      <- sd(hd1_tmp[!is.infinite(hd1_tmp)], na.rm = TRUE)
    hd2_sd[index]      <- sd(hd2_tmp[!is.infinite(hd2_tmp)], na.rm = TRUE)
    hd_mean[index]     <- mean(hd_tmp[!is.infinite(hd_tmp)], na.rm = TRUE)
    hd_sd[index]       <- sd(hd_tmp[!is.infinite(hd_tmp)], na.rm = TRUE)
    
    s_01[index]        <- sum(is.na(foreach_result_matrix[5, ])) / n_sims
    s_02[index]        <- sum(is.na(foreach_result_matrix[6, ])) / n_sims
    
    taed_mean[index] <- mean(foreach_result_matrix[7, ], na.rm = TRUE)
    taed_sd[index] <- sd(foreach_result_matrix[7, ], na.rm = TRUE)
    
    cat(sprintf("%2.2f percent done\n", index / n_iter * 100))
    
  }
}
parallel::stopCluster(cl)

# save results to data frame then to csv file
params                <- list(N = sample_sizes, T = time_periods)
result_df             <- expand.grid(params)[c(2, 1)]

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
result_df$taed_mean <- taed_mean
result_df$taed_sd <- taed_sd

# write to file

date_time_str = substr(gsub(" ", "-", gsub(":", "-", as.character(Sys.time()))), 1, 16)
output_dir = file.path("bld", "R")
if (test_run) {
  file_name = file.path(output_dir, paste0("simulation_dgp1_", date_time_str, ".csv"))
} else{
  file_name = file.path(output_dir, "simulation_dgp1.csv")
}
write_csv(result_df, file_name)
