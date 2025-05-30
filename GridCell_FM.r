GridCell_FM(archivo,path){ 
  library(FMM)
  library(readxl)
  library(reticulate)
  library(fields)
  # Load auxiliary functions
  source(paste0(path, "/FMM_ECG3D_Codes/auxMultiFMM.R"), chdir = TRUE)
  
  np <- import("numpy")
  datos <- np$load(archivo, allow_pickle = TRUE)
  mod2 <- datos$f[["spikes_mod2"]]
  ###Occupancy: same across neurons###
  heatmapsT <- array(0,c(168,200,200))  # here we store the data
  Z_mat <- array(0,c(168,200,200))     # here the predicted values
  
  ###analysis per neuron###
  for (exp in 1:168) {
    
    t <- py_to_r(datos$f[["t"]])
    x <- py_to_r(datos$f[["x"]])
    y <- py_to_r(datos$f[["y"]])
    # Convert to R object (if from Python)
    t_spike <- py_to_r(mod2[[1]][[exp]])
    # Filter for valid time range
    t_spike <- t_spike[t_spike >= min(t) & t_spike <= max(t)]
    # Interpolation to obtain (x, y) positions at those times
    x_spike <- approx(t, x, xout = t_spike)$y
    y_spike <- approx(t, y, xout = t_spike)$y
    
    # Verify sufficient data to estimate density
    if (length(x_spike) > 1 && length(y_spike) > 1) {
      heatmap <- kde2d(x_spike, y_spike, n = 200)
      heatmapsT[exp,,] <- heatmap$z
    } else {
      heatmaps[exp,,] <- NA  # or store NA, as preferred
    }
  }
  
  dim1 <- 200
  dim2 <- 200
  kde_occupancy <- kde2d(x, y, n = 200)
  Z_mat[exp,,] <- (heatmapsT[exp,,]/(kde_occupancy$z + 5))
}
