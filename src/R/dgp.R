################################################################################
## FUNCTIONS FOR GENERATING AND EVALUATING DATA IN THE MONTE CARLO SIMULATION ##
################################################################################


################################################################################
#################################### DGPs ######################################
################################################################################


DGP <- function(T, N, beta, index) {
  if (index <= 0)
    stop("index has to be an integer between 3 and 6")
  else if (index == 1)
    stop("dgp1 (multiple regrssors) cannot be called through this wrapper function")
  else if (index == 2)
    return(dgp2(T, N, beta))
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


## Constants

ERROR_SD <- sqrt(0.5)


## Functions

dgp1 <- function(T, N) {
    # dgp1 (multiple regressors)
    
    beta_tau1 <- make_beta(T, 2, N, dyadic=FALSE) # S_1 = 2
    beta_tau2 <- make_beta(T, 3, N, dyadic=FALSE) # S_2 = 3
    
    beta11 <- rep(beta_tau1$beta, N)
    beta22 <- rep(beta_tau2$beta, N)
    
    alpha <- rnorm(N)
    alpha <- rep(alpha, each = T)
    
    X1 <- rnorm(N * T) + alpha / 2
    X2 <- rnorm(N * T) + alpha / 2
    
    ERROR_SD = sqrt(2)
    e  <- rnorm(T * N, ERROR_SD)
    
    Y  <- alpha + X1 * beta11 + X2 * beta22 + e
    
    list(
      Y = matrix(Y, nrow = T),
      X = list("x1"=matrix(X1, nrow = T), "x2"=matrix(X2, nrow = T)),
      tau1 = beta_tau1$tau,
      tau2 = beta_tau2$tau,
      beta1 = beta_tau1$beta,
      beta2 = beta_tau2$beta
    )
  }


dgp2 <- function(T, N, beta) {
  # dgp2 (endogeniety in the regressors)

  tmp   <- make_X(T, N)
  alpha <- tmp$alpha

  e <- rnorm(N * T, sd=ERROR_SD)

  Z <- alpha / 2 + rnorm(N * T)
  X <- 3 * Z + e
  
  beta <- rep(beta, N)
  
  Y <- make_Y(X, beta, alpha, 1, e)
  
  list(
    Y = matrix(Y, nrow = T),
    X = list(matrix(X, nrow = T)),
    Z = list(matrix(Z, nrow = T))
  )
}


dgp3 <- function(T, N, beta) {
  # dgp3 (heteroscedasticity in the time- and cross-section)
  
  sigma_sqrd <- runif(N * T, 1, 3)
  e     <- rnorm(N * T, ERROR_SD)
  beta  <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  Y <- make_Y(X, beta, alpha, sigma_sqrd, e)
  
  list(Y = matrix(Y, nrow = T), X = list(matrix(X, nrow = T)))
}


dgp4 <- function(T, N, beta, sd=sqrt(3)) {
  # dgp4 (heteroscedasticity in the cross-section and serial correlation)
  
  burn <- 100
  zeta <- matrix(rnorm((burn + T) * N, 0, sd), ncol = N)
  rho  <- runif(N, .25, .75)
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
  
  list(Y = matrix(Y, nrow = T), X = list(matrix(X, nrow = T)))
}


dgp5 <- function(T, N, beta) {
  # dgp5 (heterosc. in the time- and cross-section and time-fixed effect)
  
  sigma_sqrd <- runif(N * T, 1, 2)
  e     <- rnorm(N * T, ERROR_SD)
  beta  <- rep(beta, N)
  
  tmp   <- make_X(T, N)
  X     <- tmp$X
  alpha <- tmp$alpha
  
  theta <- make_time_effect(T)
  
  Y <- make_Y(X, beta, alpha, sigma_sqrd, e, theta)
  
  list(Y = matrix(Y, nrow = T), X = list(matrix(X, nrow = T)), time_effect = theta)
}


dgp6 <- function(T, N, beta) {
  # dgp6 (no-jumps; equals dgp4 but without jumps)
  
  .beta <- rep(1, T)
  dgp4(T, N, .beta, sd=2)
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


make_Y <- function(X, beta, alpha, sigma_sqrd, e, theta=NULL, mu=0) {
  #' Construct labels from parameters and features
  #' 
  #' @param X 1d feature vector of size (T x N)
  #' @param beta 1d slope paramater of size (T x N)
  #' @param alpha 1d individual fixed effects of size (T x N)
  #' @param sigma_sqrd 1d variance scaler of size (T x N)
  #' @param e 1d error terms of size (T x N)
  #' @param theta 1d time fixed effects of size (T) or scalar, defaults to NULL
  #' @param mu intercept, defaults to 0
  
  if (is.null(theta)) {
    theta <- 0
  } else {
    N <- as.integer(length(alpha) / length(theta))
    theta <- rep(theta, N)
  }
  
  y <- mu + X * beta + alpha + theta + sqrt(sigma_sqrd) * e
  return(y)
}


make_tau <- function(T, S, dyadic=FALSE) {
  tau <- numeric(S)
  for (j in 1:S) {
    if (dyadic) {
      tau[j] <- floor((T - 1) / 2 ** j) + 1
    } else {
      tau[j] <- floor(j * (T - 1) / (S + 1))   
    }
  }
  tau <- sort(tau)
  return(tau)
}


make_beta <- function(T, S, N, dyadic=FALSE) {
  if (N == 30) {
    magnitude <- 7
  } else if (N == 60) {
    magnitude <- 5
  } else if (N == 120) {
    magnitude <- 4
  } else if (N == 300) {
    magnitude <- 3
  } else {
    stop("N not in correct set.")
  }
  magnitude <- magnitude / 3

  if (S == 0) {
    beta <- rep(magnitude, T)
    tau  <- list()
  } else {
    tau <- make_tau(T, S, dyadic)
    rep_beta <- diff(c(0, tau, T))
    betas    <- magnitude * (-1) ^ (1:(S + 1))
    
    beta     <- rep(betas, times = rep_beta)
  }
  list(beta = beta, tau = tau)
}


make_time_effect <- function(T) {
  theta <- make_beta(T, floor(T / 10), 30)$beta
  return(theta)
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


beta_to_gamma <- function(beta, time_effect = NULL) {
  if (typeof(beta) == "list") {
    # multivariable case
    time_periods <- length(beta[[1]])
    gamma <- c()
    for (.beta in beta) {
      gamma <- cbind(gamma, .beta[-1], .beta[-time_periods])
    }
  } else {
    # case of one regressor
    time_periods <- length(beta)
    gamma <- cbind(beta[-1], beta[-time_periods])
  }
  
  if (is.null(time_effect)) {
    gamma <- cbind(gamma, 0)
  } else {
    delta_theta <- diff(time_effect)
    gamma <- cbind(gamma, delta_theta)
  }
  
  return(gamma)
}


dist_euclidean_time_average <- function(gamma_hat, gamma) {
  dist <- 0
  for (t in nrow(gamma)) {
    dist <- dist + sum((gamma_hat[t, ] - gamma[t, ])**2)
  }
  
  return(dist)
}
