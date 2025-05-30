
SIM_S_FMM(K,L,nit,noise_level,path,n_theta,n_phi,par_current,beta_hat){ 
  library(circular)
  library(ggplot2)
  library(imager)
  library(fields)
  
  # Load auxiliary functions
  source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM.R"), chdir = TRUE)
  
  # Classic MSE with angular wrap
  mse_wrap <- function(theta, theta_hat) {
    diffs <- wrap(theta - theta_hat)
    mean(diffs^2)
  }
  
  # Assume known alpha_P, omega_P, alpha_T, omega_T
  phi_mobius <- function(tau, alpha, omega) {
    2 * atan(omega * tan((tau - alpha) / 2))
  }
  
  # With blocks in both directions
  fit_fmm_linear5 <- function(theta, phi, f_data, 
                              alpha_P, omega_P, 
                              alpha_T, omega_T, 
                              M = 0) {
    K <- length(alpha_P)
    L <- length(alpha_T)
    
    n_theta <- length(theta)
    n_phi <- length(phi)
    
    # Compute phases
    Phi_P <- sapply(1:K, function(k) phi_mobius(theta, alpha_P[k], omega_P[k]))
    Phi_T <- sapply(1:L, function(l) phi_mobius(phi, alpha_T[l], omega_T[l]))
    
    # Compute total number of design columns
    ncols <- 4 * K * L + 2 * L + 2 * K
    
    X <- matrix(0, nrow = n_theta * n_phi, ncol = ncols)
    col_idx <- 1
    
    # Phantom block k = 0
    for (l in 1:L) {
      phi_l <- Phi_T[, l]
      base1 <- outer(rep(1, n_theta), cos(phi_l))
      base3 <- outer(rep(1, n_theta), sin(phi_l))
      
      X[, col_idx]     <- as.vector(base1)
      X[, col_idx + 1] <- as.vector(base3)
      col_idx <- col_idx + 2
    }
    
    # Block l = 0
    for (k in 1:K) {
      theta_k <- Phi_P[, k]
      base1 <- outer(cos(theta_k), rep(1, n_phi))
      base3 <- outer(sin(theta_k), rep(1, n_phi))
      
      X[, col_idx]     <- as.vector(base1)
      X[, col_idx + 1] <- as.vector(base3)
      col_idx <- col_idx + 2
    }
    
    # Blocks k = 1:K, l = 1:L
    for (k in 1:K) {
      theta_k <- Phi_P[, k]
      for (l in 1:L) {
        phi_l <- Phi_T[, l]
        base1 <- outer(cos(theta_k), cos(phi_l))
        base2 <- outer(sin(theta_k), cos(phi_l))
        base3 <- outer(cos(theta_k), sin(phi_l))
        base4 <- outer(sin(theta_k), sin(phi_l))
        
        X[, col_idx]     <- as.vector(base1)
        X[, col_idx + 1] <- as.vector(base2)
        X[, col_idx + 2] <- as.vector(base3)
        X[, col_idx + 3] <- as.vector(base4)
        col_idx <- col_idx + 4
      }
    }
    
    # Response vector
    y <- as.vector(f_data - M)
    X_int <- cbind(X, intercept = 1)
    
    # Fit model
    fit <- lm.fit(X_int, y)
    beta_hat <- fit$coefficients
    y_hat <- X_int %*% beta_hat + M
    
    # R²
    ss_res <- sum((y - (y_hat - M))^2)
    ss_tot <- sum((y - mean(y))^2)
    r2 <- 1 - ss_res / ss_tot
    
    # Output
    f_hat <- matrix(y_hat, nrow = n_theta, ncol = n_phi)
    
    list(fitted = f_hat, r2 = r2, coefficients = beta_hat)
  }
  
  # Function to refine parameters
  Ref <- function(alpha2, omega2, alpha1, omega1,
                  theta, phi, f_data,
                  tole = 0.05, it = 5,
                  fit_fmm_linear5, phi_mobius, M = 0) {
    
    # Rebuild initial parameter vector
    par_init <- c(alpha2, omega2, alpha1, omega1)
    par_current <- par_init
    r2_current <- -Inf
    K <- length(alpha2)
    L <- length(alpha1)
    
    # Generate pairs to optimize
    generate_pairs <- function(K, L) {
      pairs <- list()
      idx <- 1
      for (i in 1:K) {
        pairs[[idx]] <- c(i, i + K)
        idx <- idx + 1
      }
      for (j in 1:L) {
        start_index <- 1 + 2 * K
        pairs[[idx]] <- c(start_index + j - 1, start_index + j - 1 + L)
        idx <- idx + 1
      }
      return(pairs)
    }
    
    pairs_to_update <- generate_pairs(K, L)
    
    for (iter in 1:it) {
      for (pair in pairs_to_update) {
        mask <- rep(FALSE, length(par_init))
        mask[pair] <- TRUE
        
        make_obj_fmm <- function(theta, phi, f_data, par_init, mask, M = 0) {
          function(par_free) {
            full_par <- par_init
            full_par[mask] <- par_free
            alpha_P <- full_par[1:K]
            omega_P <- full_par[(K+1):(2*K)]
            alpha_T <- full_par[(2*K+1):(2*K+L)]
            omega_T <- full_par[(2*K+L+1):(2*K+2*L)]
            
            result <- tryCatch({
              res <- fit_fmm_linear5(theta, phi, f_data,
                                     alpha_P, omega_P,
                                     alpha_T, omega_T, M)
              -res$r2
            }, error = function(e) Inf)
            
            return(result)
          }
        }
        
        obj <- make_obj_fmm(theta, phi, f_data, par_current, mask)
        par_free_init <- par_current[mask]
        lower <- par_current * (1 - tole)
        upper <- par_current * (1 + tole)
        
        idx_alpha_P <- 1:K
        idx_alpha_T <- (2 * K + 1):(2 * K + L)
        idx_omega_P <- (K + 1):(2 * K)
        idx_omega_T <- (2 * K + L + 1):(2 * K + 2 * L)
        upper[idx_omega_P] <- pmin(1, upper[idx_omega_P])
        upper[idx_omega_T] <- pmin(1, upper[idx_omega_T])
        
        optim_res <- optim(par = par_free_init,
                           fn = obj,
                           method = "L-BFGS-B",
                           lower = lower[mask],
                           upper = upper[mask])
        
        par_new <- par_current
        par_new[mask] <- optim_res$par
        par_new[1:K] <- par_new[1:K] %% (2*pi)
        par_new[(2*K+1):(2*K+L)] <- par_new[(2*K+1):(2*K+L)] %% (2*pi)
        
        res_new <- fit_fmm_linear5(theta, phi, f_data,
                                   par_new[1:K], par_new[(K+1):(2*K)],
                                   par_new[(2*K+1):(2*K+L)],
                                   par_new[(2*K+L+1):(2*K+2*L)], M)
        r2_new <- res_new$r2
        
        if (r2_new > r2_current) {
          par_current <- par_new
          r2_current <- r2_new
        }
      }
    }
    
    # Build final design and fit model
    n_theta <- length(theta)
    n_phi <- length(phi)
    
    alpha_P <- par_current[1:K]
    omega_P <- par_current[(K+1):(2*K)]
    alpha_T <- par_current[(2*K+1):(2*K+L)]
    omega_T <- par_current[(2*K+L+1):(2*K+2*L)]
    
    Phi_P <- sapply(1:K, function(k) phi_mobius(theta, alpha_P[k], omega_P[k]))
    Phi_T <- sapply(1:L, function(l) phi_mobius(phi, alpha_T[l], omega_T[l]))
    
    X <- matrix(0, nrow = n_theta * n_phi, ncol = 4*K*L + 2*L + 2*K)
    col_idx <- 1
    
    for (l in 1:L) {
      phi_l <- Phi_T[, l]
      X[, col_idx]     <- as.vector(outer(rep(1, n_theta), cos(phi_l)))
      X[, col_idx + 1] <- as.vector(outer(rep(1, n_theta), sin(phi_l)))
      col_idx <- col_idx + 2
    }
    
    for (k in 1:K) {
      theta_k <- Phi_P[, k]
      X[, col_idx]     <- as.vector(outer(cos(theta_k), rep(1, n_phi)))
      X[, col_idx + 1] <- as.vector(outer(sin(theta_k), rep(1, n_phi)))
      col_idx <- col_idx + 2
    }
    
    for (k in 1:K) {
      for (l in 1:L) {
        phi_k <- Phi_P[, k]
        phi_l <- Phi_T[, l]
        X[, col_idx]     <- as.vector(outer(cos(phi_k), cos(phi_l)))
        X[, col_idx + 1] <- as.vector(outer(sin(phi_k), cos(phi_l)))
        X[, col_idx + 2] <- as.vector(outer(cos(phi_k), sin(phi_l)))
        X[, col_idx + 3] <- as.vector(outer(sin(phi_k), sin(phi_l)))
        col_idx <- col_idx + 4
      }
    }
    
    X <- cbind(X, intercept = 1)
    y <- as.vector(f_data - M)
    fit <- lm.fit(X, y)
    beta_hat <- fit$coefficients
    y_hat <- X %*% beta_hat + M
    
    return(list(
      par_current = par_current,
      X = X,
      r2 = r2_current,
      beta_hat = beta_hat,
      y_hat = matrix(y_hat, nrow = n_theta, ncol = n_phi)
    ))
  }
}
