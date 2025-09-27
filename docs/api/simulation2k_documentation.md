# Simulation Script Documentation (simulation2k.R)

## Overview
This script demonstrates the complete workflow for fluorescence lifetime decomposition using simulated data with two molecular species (dyes). It includes data generation, EM estimation, model comparison, and comprehensive visualization.

## Purpose
- Generate synthetic fluorescence lifetime data with known ground truth
- Apply EM algorithm to recover decay parameters
- Compare k=2 vs k=3 model performance
- Visualize results including animations and diagnostic plots

## Key Components

### 1. Simulation Setup
```r
# Simulation parameters
H <- 50       # Grid dimensions (50x50 pixels)
W <- 50
m <- 2        # Number of dyes
lambda_true <- c(0.5, 1.5)  # Known decay rates
```

### 2. Spatial Mask Definition
- **Dye 1**: Circular region (radius 15) centered at (25,25)
- **Dye 2**: Rectangular region covering rows 20-40, columns 10-30

### 3. Data Generation Process
```r
# For each dye in each pixel:
# 1. Sample photon count (5-30 photons)
# 2. Generate exponential lifetimes: τ ~ Exponential(λ_d)
# 3. Combine all photons into single dataset
```

### 4. Photon Count Animation
- Visualizes temporal decay as discrete time frames
- Shows photon survival over time (τ > t)
- Creates animated GIF showing decay dynamics

### 5. EM Algorithm Implementation
```r
# Key features:
# - Pixel-specific mixing proportions (π_ij)
# - Global decay rates (λ_d)
# - 3x3 spatial neighborhood regularization
# - Iterative parameter updates with convergence monitoring
```

### 6. Model Comparison
- Fits both k=2 and k=3 models to same data
- Compares BIC values for model selection
- Evaluates reconstruction accuracy

### 7. Visualization Components

#### A. Composite λ Maps
- **True map**: Ground truth weighted by dye presence
- **Estimated map**: EM recovered parameters

#### B. Mixing Proportion Decomposition
- **True π**: Calculated from actual photon counts
- **Estimated π**: EM recovered mixing proportions
- **Error analysis**: Mean absolute error per dye

#### C. Convergence Diagnostics
- **Log-likelihood evolution**: Monitors EM convergence
- **Parameter evolution**: Decay rate estimates over iterations

#### D. Reconstruction Animation
- **Three-panel view**: Composite + individual dye contributions
- **Time evolution**: Frame-by-frame decay visualization
- **Fixed color scale**: Consistent visualization across frames

## Usage Instructions

### Basic Execution
```r
# Set working directory to script location
setwd("/path/to/Fluorescence-Lifetime-Decomposition")

# Source the simulation script
source("simulation2k.R")

# Script runs automatically and generates:
# - Animated plots (if interactive session)
# - Static comparison plots
# - Convergence diagnostics
# - Model comparison statistics
```

### Parameter Modification
```r
# Adjust simulation parameters at top of script:
H <- 100          # Larger grid (100x100)
W <- 100
lambda_true <- c(0.3, 2.0)  # Different decay rates
T_frames <- 15    # More time frames
```

### Output Interpretation

#### Expected Results (Typical)
```
Final Comparison of Decay Rates:
  Dye True_Rate Estimated_Rate Absolute_Error
1 Dye1       0.5          0.52          0.02
2 Dye2       1.5          1.48          0.02

Mixing Proportions Decomposition Error:
  Dye  MAE
1 Dye1 0.05
2 Dye2 0.04

BIC for k=2 model: 12456.7
BIC for k=3 model: 12892.1  # Higher BIC indicates worse fit
```

#### Visualization Guidelines
1. **λ Maps**: Look for spatial correspondence between true/estimated patterns
2. **Convergence Plots**: Check for log-likelihood plateau (good convergence)
3. **Parameter Evolution**: Monitor stabilization of decay rate estimates
4. **Reconstruction Animation**: Verify temporal decay matches expectations

## Technical Details

### Algorithm Parameters
- **max_iter**: 30 iterations (typically sufficient)
- **Neighborhood size**: 3×3 pixels for spatial smoothing
- **Initialization**: Uniform mixing proportions, random decay rates in [1,3]

### Computational Complexity
- **Time**: O(H×W×m×iterations×photons_per_pixel)
- **Memory**: O(H×W×m) for parameter storage
- **Typical runtime**: 2-5 minutes for 50×50 grid

### Convergence Monitoring
```r
# Check convergence by examining:
loglikelihoods    # Should plateau
lambda_hist       # Should stabilize
final_vs_initial  # Parameter improvement
```

## Troubleshooting

### Common Issues

1. **Poor Convergence**
   - Increase `max_iter`
   - Check initialization range for `lambda_hat`
   - Verify data preprocessing

2. **Spatial Artifacts**
   - Adjust neighborhood size
   - Check for edge effects in boundary pixels

3. **Memory Issues**
   - Reduce grid size (H, W)
   - Decrease number of photons per pixel

### Parameter Sensitivity
- **Initialization range**: [0.1, 5.0] for fluorescence lifetimes
- **Neighborhood size**: 1×1 (no smoothing) to 5×5 (heavy smoothing)
- **Convergence threshold**: Monitor log-likelihood change < 0.01

## Extensions and Modifications

### Adding New Features
```r
# Example: Add regularization parameter
pi_update <- function(pi_hat, resp, neighbors, alpha = 0.1) {
  # alpha controls smoothing strength
  smoothed_pi <- (1-alpha) * pi_hat + alpha * neighborhood_average
  return(smoothed_pi)
}
```

### Custom Dye Patterns
```r
# Define new spatial masks
mask3 <- matrix(0, H, W)
mask3[10:20, 30:40] <- 1  # Third dye region
masks <- list(mask1, mask2, mask3)
```

## Integration with Other Scripts

This simulation script demonstrates the same algorithms used in:
- `FullSudy.R`: Multi-k analysis with real data
- `twoDyeDataAnalysis.R`: Complete pipeline with preprocessing

The core EM functions are identical across scripts, ensuring consistency.

## References

- **EM Algorithm**: Dempster, A. P., et al. "Maximum likelihood from incomplete data via the EM algorithm." Journal of the Royal Statistical Society. 1977.
- **Fluorescence Lifetime Theory**: Lakowicz, J. R. "Principles of Fluorescence Spectroscopy." Springer, 2006.
- **Spatial Regularization**: Besag, J. "Spatial interaction and the statistical analysis of lattice systems." Journal of the Royal Statistical Society. 1974.
