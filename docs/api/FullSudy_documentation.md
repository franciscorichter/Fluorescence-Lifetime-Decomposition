# Main Analysis Script Documentation (FullSudy.R)

## Overview
This is the primary analysis script for comprehensive fluorescence lifetime decomposition across multiple model complexities (k=1 to k=5). It provides a complete workflow from data loading through model comparison to visualization.

## Purpose
- Load and prepare real fluorescence lifetime data
- Run EM algorithm for multiple component numbers (k=1-5)
- Compare models using BIC and log-likelihood criteria
- Generate frame-by-frame reconstructions
- Create diagnostic and comparison plots

## Workflow Overview

### 1. Data Loading and Preparation
```r
# Load preprocessed lifetime data
load("df_lifetimes.RData")  # Contains df_lifetimes with columns: tau, i, j

# Prepare data structure
df_events <- df_lifetimes
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))
```

### 2. EM Algorithm Setup
- Source C++ implementation: `Rcpp::sourceCpp("em_algorithm.cpp")`
- Define wrapper function for EM estimation
- Set up parameter initialization and iteration control

### 3. Multi-Model Analysis
```r
# Run EM for k=1 to k=5
results_list <- list()
BIC_values <- numeric()
for(k in 1:5){
  res <- em_estimation(df_events, k, n_rows, n_cols, max_iter = 10)
  results_list[[k]] <- res

  # Calculate BIC
  p_params <- length(unique_pixels) * (k - 1) + k
  BIC <- -2 * res$loglikelihoods[length(res$loglikelihoods)] + p_params * log(nrow(df_events))
  BIC_values[k] <- BIC
}
```

### 4. Model Comparison and Selection
- **BIC Analysis**: Lower BIC indicates better model fit
- **Log-likelihood comparison**: Higher values indicate better fit
- **Parameter count penalty**: More complex models are penalized

### 5. Reconstruction and Visualization
- Frame-by-frame intensity prediction for each model
- Component-wise decomposition (individual dye contributions)
- Spatial mapping of estimated parameters

## Key Functions

### EM Estimation Wrapper
```r
em_estimation <- function(df_events, k, n_rows, n_cols, max_iter = 30) {
  # Data preparation
  # Parameter initialization
  # C++ EM execution
  # Return results
}
```

**Purpose**: Interface between R data structures and C++ EM implementation

**Key Operations**:
- Data structure conversion for C++ compatibility
- Parameter initialization (uniform π, random λ)
- EM execution via C++ function
- Result packaging for R analysis

### BIC Calculation
```r
# Bayesian Information Criterion
p_params <- length(unique_pixels) * (k - 1) + k  # Free parameters
BIC <- -2 * loglik + p_params * log(n_observations)
```

**Purpose**: Model selection accounting for both fit quality and complexity

**Interpretation**:
- Lower BIC = Better model
- ΔBIC > 10 = Strong evidence against higher-k model
- Accounts for overfitting in complex models

### Reconstruction Engine
```r
# For each time frame t:
for(t in 1:n_frames){
  for(each pixel){
    # Predict intensity using EM parameters
    intensity = n_ij * Σ_d π_d * exp(-λ_d * t)

    # Store composite and component-specific predictions
    composite_img[i,j] <- intensity
    component_imgs[d,i,j] <- n_ij * π_d * exp(-λ_d * t)
  }
}
```

**Purpose**: Generate model predictions for visualization and validation

## Output Components

### 1. Model Comparison Results
```r
# Example output:
# Final estimated decay rates (k = 2):
# [1] 0.52 1.48
# BIC for k = 2 model: 15432.1

# BIC vs. Number of Components plot
# Log-Likelihood vs. Number of Components plot
```

### 2. Reconstruction Visualizations
- **Composite images**: Model-predicted total intensity per frame
- **Component images**: Individual dye contributions per frame
- **Spatial parameter maps**: Estimated decay rates and mixing proportions

### 3. Diagnostic Information
- **Convergence tracking**: Log-likelihood evolution
- **Parameter history**: Decay rate estimates over iterations
- **Pixel statistics**: Coverage and data quality metrics

## Usage Instructions

### Basic Execution
```r
# Ensure all files are in working directory:
# - FullSudy.R
# - em_algorithm.cpp
# - df_lifetimes.RData

# Source and run
source("FullSudy.R")

# Results automatically saved to:
# - em_results_multi_k.RData
# - Generated plots displayed interactively
```

### Parameter Modification
```r
# Adjust at top of script:
max_iter <- 30      # Increase for better convergence
n_frames <- 10      # Adjust based on data
min_photons <- 5    # Data quality threshold
```

### Customization Options
```r
# Model range
k_values <- 1:5     # Can change to 2:4 for faster execution

# Initialization
lambda_range <- c(0.1, 5.0)  # Fluorescence lifetime range
pi_initialization <- "uniform"  # or "random"

# Output control
save_plots <- TRUE
verbose <- TRUE
```

## Data Requirements

### Input Data Format
- **df_lifetimes**: Data frame with columns:
  - `tau`: Lifetime values (discrete time units)
  - `i`: Row coordinates (1-based)
  - `j`: Column coordinates (1-based)

### Data Quality Expectations
- **Photon coverage**: >100 photons per pixel for reliable estimates
- **Lifetime range**: Typically 1-20 time units
- **Spatial coverage**: Complete grid or well-defined regions

## Performance Characteristics

### Computational Requirements
- **Time**: 5-15 minutes for full k=1-5 analysis (50×50 grid)
- **Memory**: ~100-200MB for typical datasets
- **Storage**: ~50MB for saved results

### Scalability
- **Grid size**: Linear scaling with number of pixels
- **Photon count**: Linear scaling with data size
- **Model complexity**: Linear scaling with k

## Output Files

### Saved Results (`em_results_multi_k.RData`)
```r
# Contains:
results_list       # EM results for each k
BIC_values         # BIC scores for model comparison
recon_results      # Reconstruction data for visualization
# Usage: load("em_results_multi_k.RData")
```

### Generated Plots
1. **BIC vs. k**: Model selection plot
2. **Log-likelihood vs. k**: Fit quality comparison
3. **Reconstruction frames**: Visual validation
4. **Parameter evolution**: Convergence diagnostics

## Interpretation Guide

### Model Selection
```r
# Example BIC comparison:
BIC_values <- c(18000, 15432, 15678, 15901, 16123)
# k=2 has lowest BIC → Best model
best_k <- which.min(BIC_values)
```

### Parameter Validation
```r
# Check parameter reasonableness:
lambda_ranges <- c(0.1, 5.0)  # Typical fluorescence lifetimes
pi_ranges <- c(0, 1)          # Valid mixing proportions
convergence_check <- loglik_diff < 0.01  # Convergence criterion
```

### Spatial Validation
- **Mixing proportion maps**: Should show coherent spatial patterns
- **Reconstruction quality**: Predicted vs. observed decay curves
- **Edge effects**: Check boundary pixel behavior

## Troubleshooting

### Common Issues

1. **Convergence Problems**
   - Increase `max_iter`
   - Check initialization parameters
   - Verify data quality (photon counts)

2. **Memory Issues**
   - Reduce grid size for large images
   - Use data subsampling for initial analysis
   - Clear intermediate variables

3. **Poor Model Fit**
   - Check data preprocessing
   - Verify lifetime calculation
   - Consider different k range

### Diagnostic Checks
```r
# Convergence check
final_ll <- results_list[[k]]$loglikelihoods[max_iter]
previous_ll <- results_list[[k]]$loglikelihoods[max_iter-1]
convergence <- abs(final_ll - previous_ll) / abs(final_ll) < 0.01

# Parameter stability
lambda_stable <- sd(results_list[[k]]$lambda_hist[max_iter-5:max_iter,]) < 0.1
```

## Integration with Other Scripts

### Complementary Scripts
- **`simulation2k.R`**: Use for validation with known ground truth
- **`twoDyeDataAnalysis.R`**: Alternative pipeline with different preprocessing

### Workflow Integration
```r
# Typical analysis sequence:
# 1. Run simulation2k.R for method validation
# 2. Run twoDyeDataAnalysis.R for real data preprocessing
# 3. Run FullSudy.R for comprehensive multi-model analysis
```

## Extensions and Customization

### Adding New Model Selection Criteria
```r
# Example: AIC calculation
aic_values <- numeric(5)
for(k in 1:5){
  p_params <- length(unique_pixels) * (k - 1) + k
  aic_values[k] <- -2 * final_ll + 2 * p_params
}
```

### Custom Visualization
```r
# Add spatial smoothing to results
smooth_lambda_map <- gaussian_smooth(estimated_lambda_map, sigma = 1.0)
```

### Parameter Constraints
```r
# Add bounds to decay rates
lambda_bounds <- c(0.1, 10.0)  # Physically reasonable range
lambda_hat <- pmax(pmin(lambda_hat, lambda_bounds[2]), lambda_bounds[1])
```

## References and Further Reading

- **Model Selection**: Burnham, K. P. & Anderson, D. R. "Model Selection and Multimodel Inference"
- **EM Algorithm**: Bilmes, J. A. "A Gentle Tutorial of the EM Algorithm"
- **Fluorescence Lifetime Imaging**: Becker, W. "Fluorescence Lifetime Imaging - Techniques and Applications"
- **Bayesian Information Criterion**: Schwarz, G. "Estimating the Dimension of a Model"
