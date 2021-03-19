################################################################################
## FUNCTIONS FOR GENERATING AND EVALUATING DATA IN THE MONTE CARLO SIMULATION ##
################################################################################


################################################################################
#################################### DGPs ######################################
################################################################################


DGP <- function(T, N, beta, index) {
  if (index <= 0)
    stop("index has to be a integer between 3 and 6")
  else if (index == 1)
    stop("dgp1 (multiple regrssors) cannot be called through this wrapper function")
  else if (index == 2)
    stop("dgp2 (endogenity) cannot be called through this wrapper function")
  else if (index == 3)
    return(dgp3(T, N, beta))
  else if (index == 4)
    return(dgp4(T, N, beta))
  else if (index == 5)
    return(dgp5(T, N, beta))
  else if (index == 6)
    return(dgp6(T, N, beta))
  else
    DGP(0, 0, 0, 0)
}


dgp1 <- function(T, N, a11 = 0.5, a12 = 0.5, a21 = 0.5, a22 = 0.5) {
    # dgp1 (multiple regressors)
    
    beta_tau1 <- make_beta(T, 2) # S_1==2
    beta_tau2 <- make_beta(T, 3) # S_2==3
    
    beta11 <- rep(beta_tau1$beta, N)
    beta22 <- rep(beta_tau2$beta, N)
    
    alpha <- rnorm(N)
    alpha <- rep(alpha, each = T)
    
    X1 <- a11 * rnorm(N * T) + a12 * alpha
    X2 <- a21 * rnorm(N * T) + a22 * alpha
    
    e  <- rnorm(T * N, 0, sd = sqrt(.5))
    Y  <- alpha + X1 * beta11 + X2 * beta22 + e
    
    list(
      Y = matrix(Y, nrow = T),
      X1 = matrix(X1, nrow = T),
      X2 = matrix(X2, nrow = T),
      tau1 = beta_tau1$tau,
      tau2 = beta_tau2$tau,
      beta1 = beta_tau1$beta,
      beta2 = beta_tau2$beta
    )
  }


dgp2 <- function(T, N, beta) {
  # dgp2 (endogeniety in the regressors)
  
  stop("This dgp needs to be checked by dominik.")

  tmp   <- make_X(T, N)
  alpha <- tmp$alpha

  e <- matrix(rnorm(2 * N * T), ncol=2)
  variance = matrix(c(1, 1/2, 1/2, 1), ncol=2)
  L <- chol(variance)
  e <- e %*% t(L)

  Z <- alpha / 2
  X <- Z + e[,1]
  
  beta  <- rep(beta, N)
  
  Y <- make_Y(X, beta, alpha, 1, e[,2])
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T), Z = matrix(Z, nrow = T))
}


dgp3 <- function(T, N, beta) {
  # dgp3 (heteroscedasticity in the time- and cross-section)
  
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
  # dgp4 (heteroscedasticity in the cross-section and serial correlation)
  
  burn <- 50
  zeta <- matrix(rnorm((burn + T) * N, 0, sqrt(2)), ncol = N)
  rho  <- runif(N, 0, .5)
  e    <- matrix(NA, nrow = T + burn, ncol = N)
  
  e[1,] <- zeta[1,]
  for (t in 2:(T + burn))
    e[t,] <- rho * e[t - 1,] + zeta[t,]
  
  e    <- e[(burn + 1):(burn + T),]
  e    <- as.vector(e)
  beta <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  Y     <- make_Y(X, beta, alpha, 1, e)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}


dgp5 <- function(T, N, beta) {
  # dgp5 (heterosc. in the time- and cross-section and time-fixed effect)
  
  gamma <- runif(N * T, 1, 2)
  e     <- rnorm(N * T)
  beta  <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  theta <- sin(seq(T))
  theta <- rep(theta, each = N)
  
  Y <- make_Y(X, beta, alpha, gamma, e) + theta  # please check if this is
                                                 # correctly added (theta)
  
  list(Y = matrix(Y, nrow = T), X = matrix(X, nrow = T))
}


dgp6 <- function(T, N, beta) {
  # dgp6 (no-jumps; equals dgp4 but without jumps)
  
  .beta <- rep(2, T)
  dgp4(T, N, .beta)
}


################################################################################
############################ DGP Helper Functions ##############################
################################################################################


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


make_beta <- function(T, S) {
  if (S == 0) {
    beta <- rep(-2 / 3, T)
    tau  <- list()
  }
  else {
    tau      <- make_tau(T, S)
    rep_beta <- diff(c(0, tau, T))
    betas    <- (3 / 2) * (-1) ^ (1:(S + 1))
    
    beta     <- rep(betas, times = rep_beta)
    beta
  }
  list(beta = beta, tau = tau)
}


################################################################################
#################################### Metrics ###################################
################################################################################


MDCJ <- function(tau, tau_estimate, S) {
  # Minimum distance CJ ?
  
  if (length(tau_estimate) == 0 || is.na(tau_estimate)) {
    dist = Inf
  } else {
    dist <- 0L
    
    for (j in 1:S) {
      dist <- dist + min(abs(tau_estimate - tau[j]))
    }
  }
  
  return(dist)
}


dist_hausdorff <- function(set1, set2) {
  # Hausdorff distance adjusted for our case.
  
  if (length(set1) == 0 && length(set2) == 0) {
    dist = 0
  } else if (is.na(set1) && is.na(set2)) {
    dist = 0
  } else if (length(set1) == 0 || length(set2) == 0) {
    dist = Inf
  } else if (is.na(set1) || is.na(set2)) {
    dist = Inf
  } else {
    dist = pracma::hausdorff_dist(set1, set2)
  }
  
  return(dist)
}
