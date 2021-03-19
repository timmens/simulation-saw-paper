# dgp 2 (endogeneity)

# devtools::install_github("https://github.com/timmens/sawr", force = TRUE)
library("sawr")
library("readr")
library("foreach")
library("doParallel")

source("src/R/dgp.R")  # exports function dgp2

jumps <- c(1, 2, 3)  # S
sample_sizes <- c(30, 60, 120, 300)  # N
time_periods <- 2 ^ c(5, 6, 7) + 1  # T

n_sims <- 8  # 1000
if (n_sims != 1000) warning("Check n_sims.")

n_iter <- length(sample_sizes) * length(time_periods) * length(jumps)

s_est_mean  <- numeric(n_iter)
s_est_sd    <- numeric(n_iter)
mdcj_mean   <- numeric(n_iter)
mdcj_sd     <- numeric(n_iter)
mise_mean   <- numeric(n_iter)
mise_sd     <- numeric(n_iter)
hd_mean     <- numeric(n_iter)
hd_sd       <- numeric(n_iter)
s_0         <- numeric(n_iter)

seed <- 123
set.seed(seed)
starting_time <- Sys.time()

cl  <- parallel::makeCluster(4)
doParallel::registerDoParallel(cl)

for (t in time_periods) {
  
  for (s in jumps) { 
    
    true_beta_tau <- make_beta(t, s)
    
    for (n in sample_sizes) { 
      
      cat(sprintf("dgp = %d; n = %d; s = %d; t = %d\n", 2, n, s, t))
      
      tmp_result_matrix <- foreach::foreach(r = 1:n_sims, .combine = cbind) %dopar% {
        
        # create data
        data  <- dgp2(t, n, true_beta_tau$beta)
        Y <- data$Y
        Z <- data$Z
        
        # apply method
        stop("Instrument method not implemented in 'saw_fun' yet.")
        results <- sawr::saw_fun(Y ~ Z, dot=FALSE)
        
        # evaluate metrics
        estimated_taus <- results$tausList[[1]]
        s_est_mean_tmp <- sum(!is.na(estimated_taus))
        mdcj_tmp       <- MDCJ(true_beta_tau$tau, estimated_taus, s)
        mise_tmp       <- mean((true_beta_tau$beta - results$betaMat)^2)
        hd_mean_tmp    <- dist_hausdorff(true_beta_tau$tau, estimated_taus)
        
        inner_loop_results <- c(
          s_est_mean_tmp,
          ifelse(is.na(mdcj_tmp), NA_real_, mdcj_tmp),
          mise_tmp,
          ifelse(is.na(hd_mean_tmp), NA_real_, hd_mean_tmp)
        )
        inner_loop_results
      }
      
      .t <- which(t == time_periods)
      .n <- which(n == sample_sizes)
      .s <- which(s == jumps)
      index <- (.t - 1) * length(jumps) * length(sample_sizes) + 
        (.s - 1) * length(sample_sizes) + .n
      
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
      
      cat(sprintf("%2.2f percent done\n", index / n_iter * 100))
    }
  }
}
parallel::stopCluster(cl)

# save results to data frame then to csv file
params                    <- list(N = N, S = S, T = T)
result_df                 <- expand.grid(params)[c(3, 2, 1)]

result_df$s_est_mean      <- s_est_mean
result_df$s_est_sd        <- s_est_sd
result_df$mise_mean       <- mise_mean
result_df$mise_sd         <- mise_sd
result_df$mdcj_mean       <- mdcj_mean
result_df$mdcj_sd         <- mdcj_sd
result_df$hd_mean         <- hd_mean
result_df$hd_sd           <- hd_sd
result_df$s_0             <- s_0

# write to file
date_time_str = substr(gsub(" ", "-", gsub(":", "-", as.character(Sys.time()))), 1, 16)
output_dir = file.path("bld", "R")
file_name = file.path(output_dir, paste0("simulation_dgp2_", date_time_str, ".csv"))
write_csv(result_df, file_name)

# write additional info
additional_info <- c(
  paste0("nsim = ", n_sims), 
  paste0("seed = ", seed), 
  paste0("ellapsed time = ", Sys.time() - starting_time)
)
file_connection <- file(file.path(output_dir, "additional_info_dgp2.txt"))
writeLines(additional_info, file_connection)
close(file_connection)