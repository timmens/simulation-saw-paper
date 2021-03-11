#########################################################################################
###### FUNCTIONS FOR GENERATING AND EVALUTATING DATA IN THE MONTE CARLO SIMULATION ######
#########################################################################################

#########################################################################################
######################################## Metrics ########################################
#########################################################################################
MDCJ <- function(tau, tau_estimate, S) {
  if(length(tau_estimate)==0){return(Inf)}
  if(is.na(tau_estimate)) return(Inf)
  out <- 0L;
  for (j in 1:S) {
    out <- out + min(abs(tau_estimate - tau[j]))
  }
  out
}

dist_hausdorff <- function(set1, set2) {
  if (length(set1) == 0 && length(set2) == 0) return(0)
  if (is.na(set1) && is.na(set2)) return(0)
  if (length(set1) == 0 || length(set2) == 0) return(Inf)
  if (is.na(set1) || is.na(set2)) return(Inf)
  #if (any(is.na(set1), is.na(set2))) return(NA_real_)
  dist1 <- max(sapply(set1, function(y) min(abs(set2 - y))))
  dist2 <- max(sapply(set2, function(x) min(abs(set1 - x))))
  
  max(dist1, dist2)
}

#########################################################################################
################################ DGP Helper Functions ###################################
#########################################################################################
make_X <- function(T, N) {
  xi    <- rnorm(T * N)
  alpha <- rnorm(N)
  alpha <- rep(alpha, each = T)
  
  X     <- 0.5 * alpha + xi
  
  list(X = X, alpha = alpha)
}

make_Y <- function(X, beta, alpha, theta, e) {
  X * beta + alpha + sqrt(theta) * e
}

make_tau <- function(T, S) {
  tau <- numeric(S)
  for (j in 1:S) {
    tau[j] <- floor(j * (T - 1) / (S + 1))
  }
  tau
}

# make_tau_dgp1 <- function(T, plus=c(0,1)) {
#   tau    <- numeric(2)
#   values <- c(1, 3) + plus
#   
#   for (j in 1:2) {
#     tau[j] <- floor((values[j] * T) / 5)
#   }
#   tau
# }

# make_beta_dgp1 <- function(T, plus=c(0,1)) {
#   tau      <- make_tau_dgp1(T, plus)
#   rep_beta <- diff(c(0, tau, T))
#   betas    <- (3 / 2) * (-1) ^ (1 : (2 + 1))
#   
#   beta     <- rep(betas, times = rep_beta)
#   beta
#   
#   list(beta = beta, tau = tau)
# }

make_beta <- function(T, S) {
  if (S == 0) {
    beta <- rep(- 2 / 3, T)
    tau  <- list()
  }
  else {
    tau      <- make_tau(T, S)
    rep_beta <- diff(c(0, tau, T))
    betas    <- (3 / 2) * (-1) ^ (1 : (S + 1))
    
    beta     <- rep(betas, times = rep_beta)
    beta
  }
  list(beta = beta, tau = tau)
}

make_instruments <- function(X, k=3) {
  gamma <- eigen(var(X))$vectors
  if(k > ncol(gamma)) k <- ncol(gamma)
  gamma_k    <- gamma[, 1:k]
  
  means      <- colMeans(X)
  mean_mat   <- matrix(means, nrow=nrow(X), ncol=ncol(X), byrow=TRUE)
  X_centered <- X - mean_mat
  beta_mat   <- X_centered %*% gamma_k
  X_cent_hat <- beta_mat %*% t(gamma_k)
  residuals  <- X_centered - X_cent_hat
  
  residuals
}

#make_2SLS_regressors <- function(X, T=NULL, k=3) {
#  if(!is.null(T)) X <- matrix(X, ncol = T)
#  if(!is.matrix(X)) stop("<<X>> has to be a matrix or argument <<T>> has to be passed to function.") 
#  Z            <- make_instruments(X, k)
#  eig          <- eigen(var(Z))
#  gamma        <- eig$vectors
#  lambda_ind   <- eig$values > .Machine$double.eps
#  gamma        <- gamma[, lambda_ind]
#  X_projection <- lm(t(X) ~ gamma)$fitted.values
#  
#  #as.vector(X_projection)
#  X_projection
#}

#########################################################################################
######################################## DGPs ###########################################
#########################################################################################

dgp1 <- function(T, N, a11=0.5, a12=0.5, a21=0.5, a22=0.5, sderr_e=sqrt(.5)) {
  ##
  beta_tau1 <- make_beta(T, 2) # S_1==2
  beta_tau2 <- make_beta(T, 3) # S_2==3
  ##
  beta11 <- rep(beta_tau1$beta, N)
  beta22 <- rep(beta_tau2$beta, N)
  ##
  alpha <- rnorm(N)
  alpha <- rep(alpha, each=T)
  ##
  X1 <- a11 * rnorm(N * T) + a12 * alpha
  X2 <- a21 * rnorm(N * T) + a22 * alpha
  
  e  <- rnorm(T * N, 0, sderr_e)
  
  Y  <- alpha + X1 * beta11 + X2 * beta22 + e
  
  list(
       Y = matrix(Y, nrow=T), 
       X1 = matrix(X1, nrow=T), X2 = matrix(X2, nrow=T), 
       tau1 = beta_tau1$tau, tau2 = beta_tau2$tau, 
       beta1 = beta_tau1$beta, beta2 = beta_tau2$beta 
  )
}

# dgp1_autocorr <- function(T, N, a11=0.5, a12=0.5, a21=0.5, a22=0.5, sderr_e=sqrt(.75)) {
#   beta_tau1 <- make_beta_dgp1(T, 0)
#   beta_tau2 <- make_beta_dgp1(T, 1)
#   
#   beta11 <- rep(beta_tau1$beta, N)
#   beta22 <- rep(beta_tau2$beta, N)
#   
#   alpha <- rnorm(N)
#   alpha <- rep(alpha, each=T)
#   
#   X1 <- a11 * rnorm(N * T) + a12 * alpha
#   X2 <- a21 * rnorm(N * T) + a22 * alpha
#   
#   burn <- 50
#   zeta <- matrix(rnorm((burn + T) * N, 0, sqrt(sderr_e)), ncol = N)
#   rho  <- 0.5
#   e    <- matrix(NA, nrow = T + burn, ncol = N)
#   
#   e[1, ] <- zeta[1, ]
#   for (t in 2:(T + burn)) e[t, ] <- rho * e[t - 1, ] + zeta[t, ]
#   
#   e    <- e[(burn + 1):(burn + T), ]
#   e    <- as.vector(e)
#   
#   Y  <- alpha + X1 * beta11 + X2 * beta22 + e
#   
#   list(
#     Y = matrix(Y, nrow=T), 
#     X1 = matrix(X1, nrow=T), X2 = matrix(X2, nrow=T), 
#     tau1 = beta_tau1$tau, tau2 = beta_tau2$tau, 
#     beta1 = beta_tau1$beta, beta2 = beta_tau2$beta 
#   )    
# }

dgp2 <- function(T, N, beta) {
  theta <- rep(runif(N, 1, 2), each = T)
  e     <- rnorm(N * T, sd = sqrt(.75))
  beta  <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  Y <- make_Y(X, beta, alpha, theta, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}

dgp3 <- function(T, N, beta) {
  theta <- runif(N * T, 1, 2)
  e     <- rnorm(N * T, sd = sqrt(.75)) 
  beta  <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  Y <- make_Y(X, beta, alpha, theta, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}

dgp4 <- function(T, N, beta) {
  burn <- 50
  zeta <- matrix(rnorm((burn + T) * N, 0, sqrt(2)), ncol = N)
  rho  <- runif(N, 0, .5)
  e    <- matrix(NA, nrow = T + burn, ncol = N)
  
  e[1, ] <- zeta[1, ]
  for (t in 2:(T + burn)) e[t, ] <- rho * e[t - 1, ] + zeta[t, ]
  
  e    <- e[(burn + 1):(burn + T), ]
  e    <- as.vector(e)
  beta <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  Y     <- make_Y(X, beta, alpha, 1, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}

dgp5 <- function(T, N, beta) {
  lambda  <- rnorm(N, 0, sqrt(.5))
  f       <- rnorm(T, 0, sqrt(.5))
  epsilon <- rnorm(T * N, 0, sqrt(.5))
  alpha   <- rnorm(N)
  mu      <- rnorm(T * N)
  
  lambda  <- rep(lambda, each = T)
  f       <- rep(f, N)
  alpha   <- rep(alpha, each = T)
  
  e       <- lambda * f + epsilon
  beta    <- rep(beta, N)
  
  X <- .3 * alpha + .5 * lambda + .3 * f + mu
  Y <- make_Y(X, beta, alpha, 1, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}

dgp6 <- function(T, N, beta) {
  burn    <- 50
  zeta_mu <- matrix(rnorm((burn + T) * N, 0, sqrt(.5)), ncol = N)
  zeta_e  <- matrix(rnorm((burn + T) * N, 0, sqrt(.5)), ncol = N)
  rho_e   <- runif(N, 0, sqrt(.5))
  rho_mu  <- runif(N, 0, sqrt(.5))
  mu      <- matrix(NA, nrow = T + burn, ncol = N)
  epsilon <- matrix(NA, nrow = T + burn, ncol = N)
  
  mu[1, ]      <- zeta_mu[1, ]
  epsilon[1, ] <- zeta_e[1, ]
  for (t in 2:(T + burn)) {
    mu[t, ]      <- rho_mu * mu[t - 1, ] + zeta_mu[t, ]
    epsilon[t, ] <- rho_e * epsilon[t - 1, ] + zeta_e[t, ]
  }
  
  mu      <- mu[(burn + 1):(burn + T), ]
  mu      <- as.vector(mu)
  epsilon <- epsilon[(burn + 1):(burn + T), ]
  epsilon <- as.vector(epsilon)
  
  lambda  <- rnorm(N, 0, sqrt(.5))
  f       <- rnorm(T, 0, sqrt(.5))
  alpha   <- rnorm(N, 0, 1)
  
  lambda  <- rep(lambda, each = T)
  f       <- rep(f, N)
  alpha   <- rep(alpha, each = T)
  
  e       <- lambda * f + epsilon 
  beta    <- rep(beta, N)
  
  X <- .3 * alpha + .3 * lambda + .3 * f + mu
  Y <- make_Y(X, beta, alpha, 1, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}

dgp7 <- function(T, N) {
  beta <- rep(2, T)
  
  dgp6(T, N, beta)
}

DGP <- function(T, N, beta, index) {
  if (index <= 0) stop("index has to be a integer between 2 and 6")
  else if (index == 1) stop("dgp1 cannot be called through this wrapper function")
  else if (index == 2) return(dgp2(T, N, beta))
  else if (index == 3) return(dgp3(T, N, beta))
  else if (index == 4) return(dgp4(T, N, beta))
  else if (index == 5) return(dgp5(T, N, beta))
  else if (index == 6) return(dgp6(T, N, beta))
  else if (index == 7) return(dgp7(T, N)) 
  else DGP(0, 0, 0, 0)
}