S_FMM <- function(K, L, vDataMatrix, path) {
  # Load required libraries
  library(ggplot2)
  library(FMM)
  library(tidyr)
  library(dplyr)
  library(readxl)
  library(akima)
  library(RcppCNPy)
  library(reticulate)
  library(fields)
  library(gplots)
  library(MASS)
  
  # Load auxiliary functions
  source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM.R"), chdir = TRUE)
  
  # Set up data
  ncom <- K
  ncom1 <- L
  Zmat <- vDataMatrix
  dim1 <- nrow(Zmat)
  dim2 <- ncol(Zmat)
  n_theta <- dim1
  n_phi <- dim2
  
  # SVD decomposition
  svd_result <- svd(Zmat)
  U <- svd_result$u
  V <- svd_result$v
  
  # Fit FMM to first singular vectors
  S1 <- fitFMM(vData = U[,1], nback = 3, omegaMin = 0.01, maxiter = 100)
  S2 <- fitFMM(vData = V[,1], nback = 3, omegaMin = 0.01, maxiter = 100)
  
  # Estimation
  ncan <- min(dim1, 40)
  M_time2 <- t(Zmat)
  ind <- seq(1, dim1, length.out = ncan)
  M_time2 <- M_time2[, round(ind)]
  
  fitml <- fitMultiFMM(vDataMatrix = M_time2[,1:ncan], nBack = ncom, maxIter = iter,
                       plotToFile = FALSE, filename = "", omegaMin = 0.01,
                       showPredeterminedPlot = FALSE)
  
  alpha1 <- fitml[[1]]$Alpha
  omega1 <- fitml[[1]]$Omega
  
  # Second stage
  ncoe <- ncom * 2 + 1
  coe <- matrix(0, nrow = dim1, ncol = ncoe)
  theta <- seq(0, 2 * pi, length.out = dim2 + 1)[-(dim2+1)]
  
  Xv <- array(0, c(ncom, 2, dim2))
  for (j in 1:ncom) {
    Xv[j,1,] <- cos(2 * atan2(omega1[j] * sin((theta - alpha1[j]) / 2), cos((theta - alpha1[j]) / 2)))
    Xv[j,2,] <- sin(2 * atan2(omega1[j] * sin((theta - alpha1[j]) / 2), cos((theta - alpha1[j]) / 2)))
  }
  
  Xmat <- matrix(0, nrow = dim2, ncol = ncom * 2)
  for (j in 1:ncom) {
    Xmat[, (j - 1) * 2 + 1] <- Xv[j,1,]
    Xmat[, (j - 1) * 2 + 2] <- Xv[j,2,]
  }
  
  rl2 <- numeric(dim1)
  for (i in 1:dim1) {
    y <- Zmat[i, ]
    ajus <- lm(y ~ Xmat)
    coe[i, ] <- coef(ajus)
    pred <- ajus$fitted.values
    ssres <- sum((Zmat[i,] - pred)^2)
    sstot <- sum((Zmat[i,] - mean(Zmat[i,]))^2)
    rl2[i] <- 1 - ssres / sstot
  }
  
  fitml1 <- fitMultiFMM(vDataMatrix = coe[,1:ncoe], nBack = ncom1, maxIter = iter,
                        plotToFile = FALSE, filename = "", showPredeterminedPlot = FALSE,
                        omegaMin = 0.01)
  
  alpha2 <- fitml1[[1]]$Alpha
  omega2 <- fitml1[[1]]$Omega
  
  # Final fit
  theta <- seq(0, 2 * pi, length.out = dim1)
  phi <- seq(0, 2 * pi, length.out = dim2)
  f_data <- Zmat
  
  final_fit <- fit_fmm_linear5(theta, phi, f_data, alpha2, omega2, alpha1, omega1, M = 0)
  
  # Refinement
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.9, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  par_current <- res$par_current
  K <- length(alpha2)
  L <- length(alpha1)
  alpha2 <- par_current[1:K]
  omega2 <- par_current[(K+1):(2*K)]
  alpha1 <- par_current[(2*K+1):(2*K+L)]
  omega1 <- par_current[(2*K+L+1):(2*K+2*L)]
  
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.05, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  par_current <- res$par_current
  alpha1 <- par_current[1:K]
  omega1 <- par_current[(K+1):(2*K)]
  alpha2 <- par_current[(2*K+1):(2*K+L)]
  omega2 <- par_current[(2*K+L+1):(2*K+2*L)]
  
  res <- Ref(alpha2, omega2, alpha1, omega1, theta, phi, f_data, tole = 0.9, it = 2,
             fit_fmm_linear5 = fit_fmm_linear5, phi_mobius = phi_mobius)
  
  len <- ncol(res$X) - 1
  r2 <- res$r2
  beta_hat <- res$beta_hat
  y_hat <- res$X %*% beta_hat + res$M
  PredFMM <- matrix(y_hat, nrow = n_theta, ncol = n_phi)
  
  return(list(
    PredFMM = PredFMM,
    r2 = r2,
    alpha1 = alpha1,
    omega1 = omega1,
    alpha2 = alpha2,
    omega2 = omega2,
    beta_hat = beta_hat
  ))
}

