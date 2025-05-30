# S-FMM
Codes to run Spatial-FMM
## How to use

We describe below three functions:  
- `S_FMM`  
- `SIM_S_FMM`  
- `GridCell_FM`  

The function `S_FMM` executes the S-FMM estimation pipeline, including SVD preprocessing, two-layer FMM estimation, linear fit, and refinement. It returns the fitted surface, parameter estimates, and fit quality.

---

### `S_FMM` function description

`S_FMM` performs a full end-to-end estimation of S-FMM from a 2D data matrix, using two-stage FMM modeling and spatial/temporal decomposition.

**Arguments**  
- `K`: integer. Number of FMM components for the first stage (temporal or per-channel).  
- `L`: integer. Number of FMM components for the second stage (spatial or across channels).  
- `vDataMatrix`: numeric matrix. A 2D data matrix where rows are spatial samples (e.g., angle θ) and columns are temporal samples (e.g., angle φ).  
- `path`: character. Path to the folder containing the script `auxMultiFMM.R`, if needed.

**Return values**  
A named R list containing:  
- `PredFMM`: matrix. The reconstructed 2D surface after final FMM modeling and refinement.  
- `r2`: numeric. The global R² value measuring model fit quality.  
- `alpha1`, `omega1`: vectors. Estimated FMM phase and frequency parameters for the first stage.  
- `alpha2`, `omega2`: vectors. Estimated FMM phase and frequency parameters for the second stage.  
- `beta_hat`: numeric vector. Final estimated regression coefficients from the 3DFMM linear model.

---

### Summary of processing steps

1. SVD decomposition of the input data matrix to extract dominant spatial/temporal modes.  
2. Stage 1 FMM fit: The leading singular vectors are fitted using the `fitMultiFMM` function.  
3. Stage 2 FMM fit: Outputs from Stage 1 are used to construct predictors for a second layer of FMM fits via `fitMultiFMM`.  
4. Stage 3: Global model fitting using the custom function `fit_fmm_linear5`, based on both `alpha` and `omega` parameters from Stage 1 and 2.  
5. Stage 4: Parameter refinement using the `Ref` function, which optimizes phase and frequency estimates for improved model accuracy.  
6. Output reconstruction: The final surface is reconstructed and the global R² is computed.

---

### Auxiliary functions

The script includes three internal auxiliary functions that support the core modeling process:  
- `phi_mobius`: Computes Möbius transformations, which are essential for phase adjustment and normalization on the unit circle.  
- `fit_fmm_linear5`: Performs global linear modeling with known `alpha` and `omega` parameters from dual-layer FMM decomposition.  
- `Ref`: Refines the parameter estimates (`alpha` and `omega`) through an iterative optimization process, enhancing model fidelity.  

These functions are embedded in the script to ensure self-containment and do not require external dependencies.
