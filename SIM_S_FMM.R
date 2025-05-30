
SIM_S_FMM(K,L,nit,noise_level,path,n_theta,n_phi,par_current,beta_hat){
  library(circular)
library(ggplot2)
library(imager)
library(fields)
  
  ################### FUNCTIONS
  # Load auxiliary functions
  source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM.R"), chdir = TRUE)
  
  # Classic MSE with angular wrap
  mse_wrap <- function(theta, theta_hat) {
    diffs <- wrap(theta - theta_hat)
    mean(diffs^2)
  }
  
  ### assume known alpha_P, omega_P, alpha_T, omega_T
  phi_mobius <- function(tau, alpha, omega) {
    2 * atan(omega * tan((tau - alpha) / 2))
  }

    #linear fit function
  fit_fmm_linear5 <- function(theta, phi, f_data, 
                              alpha_P, omega_P, 
                              alpha_T, omega_T, 
                              M = 0) {
    
    K <- length(alpha_P)
    L <- length(alpha_T)
    
    n_theta <- length(theta)
    n_phi <- length(phi)
    
    # Calculamos fases
    Phi_P <- sapply(1:K, function(k) phi_mobius(theta, alpha_P[k], omega_P[k]))  # n_theta x K
    Phi_T <- sapply(1:L, function(l) phi_mobius(phi, alpha_T[l], omega_T[l]))    # n_phi x L
    
    # Calcular número total de columnas de diseño
    ncols <- 4 * K * L+2*L+2*K
    
    
    X <- matrix(0, nrow = n_theta * n_phi, ncol = ncols)
    col_idx <- 1
    
    ## ----- BLOQUE k = 0 -----
    
    for (l in 1:L) {
      phi_l <- Phi_T[, l]
      base1 <- outer(rep(1, n_theta), cos(phi_l))
      base3 <- outer(rep(1, n_theta), sin(phi_l))
      
      X[, col_idx]     <- as.vector(base1)
      X[, col_idx + 1] <- as.vector(base3)
      
      col_idx <- col_idx + 2
    }
    
    
    ## ----- BLOQUE l = 0 -----
    
    for (k in 1:K) {
      theta_k <- Phi_P[, k]
      base1 <- outer(cos(theta_k), rep(1, n_phi))
      base3 <- outer(sin(theta_k), rep(1, n_phi))
      
      X[, col_idx]     <- as.vector(base1)
      X[, col_idx + 1] <- as.vector(base3)
      
      col_idx <- col_idx + 2
    }
    
    
    ## ----- BLOQUES k = 1:K, l = 1:L -----
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
    
    # Vector de respuesta
    y <- as.vector(f_data - M)
    #y <- as.vector(aperm(f_data, c(2,1)) - M)
    X_int <- cbind(X, intercept = 1)
    
    # Ajuste
    fit <- lm.fit(X_int, y)
    beta_hat <- fit$coefficients
    y_hat <- X_int %*% beta_hat + M
    
    # R²
    ss_res <- sum((y - (y_hat - M))^2)
    ss_tot <- sum((y - mean(y))^2)
    r2 <- 1 - ss_res / ss_tot
    
    # Salida
    f_hat <- matrix(y_hat, nrow = n_theta, ncol = n_phi)
    
    list(fitted = f_hat, r2 = r2, coefficients = beta_hat)
  }
  
  ###############refitting the parameters
  Ref<-  function(
    alpha2, omega2, alpha1, omega1,
    theta, phi, f_data,
    tole = 0.05, it = 5,
    fit_fmm_linear5, phi_mobius, M = 0
  ) {
    # Reconstruir vector inicial
    par_init <- c(alpha2, omega2, alpha1, omega1)
    par_current <- par_init
    r2_current <- -Inf
    K <- length(alpha2)
    L <- length(alpha1)
    
    # Generar pares a optimizar
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
            
            resultado <- tryCatch({
              res <- fit_fmm_linear5(theta, phi, f_data,
                                     alpha_P, omega_P,
                                     alpha_T, omega_T, M)
              -res$r2
            }, error = function(e) Inf)
            
            return(resultado)
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
    
    # Construcción del diseño y ajuste final
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
      X=X,
      r2 = r2_current,
      beta_hat = beta_hat,
      y_hat = matrix(y_hat, nrow = n_theta, ncol = n_phi)
    ))
  }
  
  
  sort_alpha <- function(alpha, omega, phi) {
    stopifnot(length(alpha) == length(omega), length(phi) == 1)
    
    ord <- order(alpha)
    alpha_sorted <- alpha[ord]
    omega_sorted <- omega[ord]
    
    delta <- (alpha_sorted - phi + 2 * pi) %% (2 * pi)
    delta[delta == 0] <- 2 * pi  # evitar que phi cuente como "siguiente"
    
    idx_start <- which.min(delta)
    if (length(idx_start) == 0) idx_start <- 1  # por seguridad
    
    # nueva forma segura de circular shift
    circ_order <- if (idx_start == 1) {
      1:length(alpha_sorted)
    } else {
      c(idx_start:length(alpha_sorted), 1:(idx_start - 1))
    }
    
    list(
      alpha = alpha_sorted[circ_order],
      omega = omega_sorted[circ_order]
    )
  }
  
  ###########################################################################
  alpha_P <- par_current[1:K]
  omega_P <- par_current[(K+1):(2*K)]
  alpha_T <- par_current[(2*K+1):(2*K+L)]
  omega_T <- par_current[(2*K+L+1):(2*K+2*L)]
theta <- seq(0, 2 * pi, length.out = n_theta)
phi   <- seq(0, 2 * pi, length.out = n_phi)
dim1<-n_theta
dim2<-n_phi
M<-0
# All the phases
Phi_P <- sapply(1:K, function(k) phi_mobius(theta, alpha_P[k], omega_P[k]))  # n_theta x K
Phi_T <- sapply(1:L, function(l) phi_mobius(phi, alpha_T[l], omega_T[l]))    # n_phi x L
#design-matrix
X <- matrix(0, nrow = n_theta * n_phi, ncol = 4 * K * L + 2 * L+2*K)
col_idx <- 1

## ----- BLOck k = 0 -----
for (l in 1:L) {
  phi_l <- Phi_T[, l]
  base1 <- outer(rep(1, n_theta), cos(phi_l))  # cos(Phi_0) * cos(Phi_l)
  base3 <- outer(rep(1, n_theta), sin(phi_l))  # cos(Phi_0) * sin(Phi_l)
  
  X[, col_idx]     <- as.vector(base1)
  X[, col_idx + 1] <- as.vector(base3)
  
  col_idx <- col_idx + 2
}


## ----- BLOck l = 0 -----

for (k in 1:K) {
  theta_k <- Phi_P[, k]
  base1 <- outer(cos(theta_k), rep(1, n_phi))
  base3 <- outer(sin(theta_k), rep(1, n_phi))
  
  X[, col_idx]     <- as.vector(base1)
  X[, col_idx + 1] <- as.vector(base3)
  
  col_idx <- col_idx + 2
}



for (k in 1:K) {
  for (l in 1:L) {
    phi_k <- Phi_P[, k]  # vector de n_theta
    phi_l <- Phi_T[, l]  # vector de n_phi
    
    # Creamos los 4 funciones base
    base1 <- outer(cos(phi_k), cos(phi_l))
    base2 <- outer(sin(phi_k), cos(phi_l))
    base3 <- outer(cos(phi_k), sin(phi_l))
    base4 <- outer(sin(phi_k), sin(phi_l))
    
    # Aplanamos y metemos en X
    X[, col_idx]     <- as.vector(base1)
    X[, col_idx + 1] <- as.vector(base2)
    X[, col_idx + 2] <- as.vector(base3)
    X[, col_idx + 3] <- as.vector(base4)
    
    col_idx <- col_idx + 4
  }
}


# Response vector

y_hat <- X %*% beta_hat 
signal <- matrix(y_hat, nrow = n_theta, ncol = n_phi)
y_hat_matrix <- matrix(y_hat, nrow = n_theta, ncol = n_phi)

image(matrix(signal, nrow = n_theta, ncol = n_phi), main = "Datos ajustados")


AL1<-array(0,c(nit,L))
OM1<-array(0,c(nit,L))
AL2<-array(0,c(nit,K))
OM2<-array(0,c(nit,K))
R2<-array(0,c(nit))

for (ii in 1:nit){
noise <- matrix(rnorm(n_theta * n_phi, 0, noise_level), nrow = n_theta, ncol = n_phi)
noise_img <- as.cimg(noise)
noise_smooth <- as.matrix(isoblur(noise_img, sigma = 1))  # sigma can be changed
# Rescaling
noise_smooth <- noise_smooth * (noise_level / sd(noise_smooth))
simul <- y_hat_matrix + noise_smooth
#simulated data
f_data<-matrix(simul, nrow = n_theta, ncol = n_phi)
#S-FMM fit
res<-S_FMM(K,L,f_data,path)

R2[ii]<-res$r2
AL2[,ii] <- res$par_current[1:K]
OM2[,ii] <- res$par_current[(K+1):(2*K)]
AL1[,ii] <- res$par_current[(2*K+1):(2*K+L)]
OM1[,ii] <- res$par_current[(2*K+L+1):(2*K+2*L)]

}
###############################################################################################
###########################################################################################


alpha_P_true <- par_current[1:K]
omega_P_true <- par_current[(K+1):(2*K)]
alpha_T_true <- par_current[(2*K+1):(2*K+L)]
omega_T_true <- par_current[(2*K+L+1):(2*K+2*L)]
# Inicialization of results vectors
mse_alpha_P <- numeric(K)
mse_omega_P <- numeric(K)
mse_alpha_T <- numeric(L)
mse_omega_T <- numeric(L)

phi <- pi
###parameter identifiability using order
res <- sort_alpha(alpha_P_true, omega_P_true, phi)
alpha_P_sorted <- res$alpha
omega_P_sorted <- res$omega


res_list <- lapply(1:nrow(AL1), function(i) {
  sort_alpha(AL1[i, ], OM1[i, ], phi)
})


AL1_sorted <- AL1 
OM1_sorted <- OM1

AL2_sorted <- AL2  
OM2_sorted <- OM2

for (i in 1:nrow(AL1)) {
  res <- sort_alpha(AL1[i,], OM1[i,], phi)
  AL1_sorted[i, ] <- res$alpha
  OM1_sorted[i, ] <- res$omega
}

for (i in 1:nrow(AL2)) {
  res <- sort_alpha(AL2[i,], OM2[i,], phi)
  AL2_sorted[i,] <- res$alpha
  OM2_sorted[i,] <- res$omega
}

#  MSE for alpha and omega
for (k in 1:K) {
  mse_alpha_P[k] <- mse_wrap(AL2_sorted[1:nit,k], alpha_P_sorted[k])
  mse_omega_P[k] <- mean((OM2_sorted[1:nit,k] - omega_P_sorted[k])^2)
  
}
for (k in 1:L) {
  mse_alpha_T[k] <- mse_wrap(AL1_sorted[1:nit,k], alpha_T_sorted[k])
  mse_omega_T[k] <- mean((OM1_sorted[1:nit,k] - omega_T_sorted[k])^2)
}


# mean R^2 
mean_r2 <- mean(R2[1:nit])
# Results
cat("MSE alpha_P:", mse_alpha_P, "\n")
cat("MSE omega_P:", mse_omega_P, "\n")
cat("MSE alpha_T:", mse_alpha_T, "\n")
cat("MSE omega_T:", mse_omega_T, "\n")
cat("Mean R2:", mean_r2, "\n")


################Confidende intervals
ci_alpha_P <- matrix(NA, K, 2)
ci_alpha_T <- matrix(NA, L, 2)
ci_omega_P <- matrix(NA, K, 2)
ci_omega_T <- matrix(NA, L, 2)

# Auxiliary function
circular_ci_around_median <- function(samples) {
  circ <- circular(samples, type = "angles", units = "radians", modulo = "2pi")
  med <- median.circular(circ)
  centered <- (circ - med) %% (2 * pi)
  centered <- circular(centered, type = "angles", units = "radians", modulo = "2pi")
  ic_centered <- quantile(centered, probs = c(0.025, 0.975))
  ic <- (ic_centered + med) %% (2 * pi)
  return(ic)
}
# ICs
for (k in 1:K) {
  # α_P (circular)
  ci_alpha_P[k, ]  <- circular_ci_around_median(AL2_sorted[1:nit, k])
  # ω_P (linear)
  ci_omega_P[k, ]  <- quantile(OM2_sorted[1:nit, k], probs = c(0.025, 0.975))
 }

for (k in 1:L) {
  # α_T (circular)
  ci_alpha_T[k, ]  <- circular_ci_around_median(AL1_sorted[1:nit, k])
  # ω_T (linear)
  ci_omega_T[k, ]  <- quantile(OM1_sorted[1:nit, k], probs = c(0.025, 0.975))
}



#Results Table
build_ci_table <- function(true_vals, ci_matrix, param_names, label) {
  df <- data.frame(
    Parameter = param_names,
    Noise_Level = label,
    True_Value = round(true_vals, 4),
    CI_Lower = round(ci_matrix[, 1], 4),
    CI_Upper = round(ci_matrix[, 2], 4)
  )
  return(df)
}


param_alpha_P <- paste0("alpha_P[", 1:K, "]")
param_omega_P <- paste0("omega_P[", 1:K, "]")
param_alpha_T <- paste0("alpha_T[", 1:L, "]")
param_omega_T <- paste0("omega_T[", 1:L, "]")

tab_alpha_P <- build_ci_table(alpha_P_sorted, ci_alpha_P, param_alpha_P)
tab_omega_P <- build_ci_table(omega_P_sorted, ci_omega_P, param_omega_P)
tab_alpha_T <- build_ci_table(alpha_T_sorted, ci_alpha_T, param_alpha_T)
tab_omega_T <- build_ci_table(omega_T_sorted, ci_omega_T, param_omega_T)

T1<- build_ci_table(alpha_P_sorted, ci_alpha_P, param_alpha_P)
T2<- build_ci_table(omega_P_sorted, ci_omega_P, param_alpha_P)
T3<-build_ci_table(alpha_T_sorted, ci_alpha_T, param_alpha_T)
T4<- build_ci_table(omega_T_sorted, ci_omega_T, param_omega_T)

rbind(T1,T2,T3,T4)
 
}