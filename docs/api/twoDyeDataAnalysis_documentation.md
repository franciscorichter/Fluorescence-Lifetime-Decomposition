# Data Analysis Script Documentation (twoDyeDataAnalysis.R)

## Overview
This script provides a complete pipeline for analyzing real fluorescence lifetime imaging data, from raw time series through preprocessing, lifetime conversion, EM estimation, and visualization.

## Purpose
- Load and preprocess raw fluorescence time series data
- Convert time series to discrete lifetime events
- Apply EM algorithm for multi-exponential decomposition
- Compare k=2 and k=3 models
- Generate comprehensive visualizations and diagnostics

## Workflow Architecture

### 1. Data Loading
```r
# Load NumPy data file
np <- import("numpy")
data <- np$load("best_fit.npy")
data1 <- data[,,1,,1,1]  # Extract relevant dimensions
```

**Input Format**: Multi-dimensional NumPy array (typically: rows × cols × time × other_dims)

### 2. Time Series Preprocessing
```r
preprocess_pixel_ts <- function(ts, min_nonzero = 3) {
  # 1. Find first frame with sufficient photons
  # 2. Extract decay phase from peak
  # 3. Enforce non-increasing behavior
  # 4. Truncate at first zero
}
```

**Purpose**: Clean and prepare raw time series for lifetime analysis

**Key Operations**:
- **Threshold filtering**: Remove low-signal regions
- **Peak detection**: Identify decay starting point
- **Monotonicity enforcement**: Ensure physically valid decay
- **Zero truncation**: Remove noise floor

### 3. Lifetime Event Conversion
```r
calculate_lifetimes_from_frames <- function(frames) {
  # Convert time series to discrete lifetime events
  # Compare consecutive frames to detect photon decay
  # Return data frame with (tau, i, j) for each decay event
}
```

**Algorithm**:
- For each time step t: compare frame[t] vs frame[t+1]
- Decay events occur where intensity decreases
- Lifetime τ = t for each decayed photon
- Accumulate events across all pixels and times

### 4. EM Algorithm Implementation

#### k=2 Model
```r
# Initialize parameters
pi_hat <- matrix(1/2, nrow = N_pix, ncol = 2)
lambda_hat <- runif(2, min = 1, max = 3)

# EM iterations with spatial regularization
for(iter in 1:max_iter){
  # E-step: compute responsibilities
  # M-step: update lambda and pi with 3x3 neighborhoods
}
```

#### k=3 Model
- Similar structure with additional component
- Automatic model comparison via BIC

### 5. Results and Visualization

#### Parameter Estimates
```r
# Example output:
# Final Parameter Estimates (k=2):
#   Dye Estimated_Rate
# 1 Dye1          0.52
# 2 Dye2          1.48

# BIC Comparison:
# BIC for k=2 model: 15432.1
# BIC for k=3 model: 15678.3
```

#### Spatial Maps
- **Composite λ maps**: Spatially resolved effective decay rates
- **Mixing proportion maps**: Relative abundance of each dye

#### Frame Reconstructions
- **Time evolution**: Frame-by-frame decay visualization
- **Component separation**: Individual dye contributions
- **Composite prediction**: Total predicted intensity

## Key Functions

### Time Series Preprocessing
```r
preprocess_pixel_ts(ts, min_nonzero = 3)
```

**Parameters**:
- `ts`: Time series vector for single pixel
- `min_nonzero`: Minimum photon threshold for valid data

**Returns**: Cleaned time series or empty vector if insufficient data

**Algorithm Details**:
1. Find first frame meeting photon threshold
2. Extract from start frame to end
3. Find intensity peak (maximum)
4. Take decay phase from peak onward
5. Enforce non-increasing constraint
6. Truncate at first zero

### Lifetime Conversion
```r
calculate_lifetimes_from_frames(frames)
```

**Input**: List of intensity matrices (one per time frame)

**Output**: Data frame with columns:
- `tau`: Lifetime (time of decay)
- `i`: Row coordinate
- `j`: Column coordinate

**Algorithm**:
- Compare consecutive frames: `diff = frame[t] - frame[t+1]`
- Decay events where `diff > 0`
- Multiple photons per pixel-frame possible
- Accumulate across all pixels and times

### EM Algorithm
- **E-step**: Soft assignment using current parameters
- **M-step**: Update global λ and local π with spatial regularization
- **Neighborhood**: 3×3 pixels for mixing proportion smoothing

## Usage Instructions

### Basic Execution
```r
# Ensure data file is available:
# best_fit.npy in working directory

# Source the analysis script:
source("twoDyeDataAnalysis.R")

# Results include:
# - Parameter estimates for k=2 and k=3
# - BIC comparison
# - Spatial parameter maps
# - Frame-by-frame visualizations
```

### Parameter Tuning
```r
# Adjust preprocessing parameters:
min_photons <- 5      # Higher = stricter filtering
max_iter <- 30        # More iterations for convergence
lambda_range <- c(0.5, 3.0)  # Adjust for expected lifetimes
```

### Data Requirements
- **Format**: NumPy .npy file with time series data
- **Dimensions**: Typically (rows, cols, time, ...)
- **Data type**: Integer photon counts
- **Quality**: Signal-to-noise ratio > 5 recommended

## Output Interpretation

### Model Selection
```r
# BIC comparison example:
BIC_2 <- 15432.1  # Lower is better
BIC_3 <- 15678.3
delta_BIC <- BIC_3 - BIC_2  # >10 indicates strong preference for k=2
```

### Parameter Validation
```r
# Reasonable parameter ranges:
lambda_range <- c(0.1, 5.0)    # Typical fluorescence lifetimes
pi_range <- c(0, 1)           # Valid mixing proportions
convergence_threshold <- 0.01  # Log-likelihood change
```

### Spatial Validation
- **Coherent patterns**: Parameter maps should show spatial structure
- **Boundary effects**: Check edge pixel behavior
- **Reconstruction quality**: Predicted vs. observed decay

## Performance and Scalability

### Computational Complexity
- **Preprocessing**: O(rows × cols × time)
- **Lifetime conversion**: O(rows × cols × time)
- **EM algorithm**: O(photons × components × iterations)
- **Visualization**: O(rows × cols × frames)

### Typical Performance
- **Preprocessing**: <1 minute
- **EM fitting**: 2-5 minutes per model
- **Visualization**: 1-3 minutes
- **Total runtime**: 10-20 minutes

### Memory Usage
- **Peak usage**: ~500MB for 100×100×50 dataset
- **Bottlenecks**: EM responsibility matrices, reconstruction arrays

## Troubleshooting

### Common Issues

1. **No Lifetime Events Detected**
   ```r
   # Solutions:
   # - Lower min_photons threshold
   # - Check data file format and dimensions
   # - Verify data preprocessing steps
   ```

2. **Poor EM Convergence**
   ```r
   # Solutions:
   # - Increase max_iter
   # - Adjust initialization parameters
   # - Check for data quality issues
   ```

3. **Memory Issues**
   ```r
   # Solutions:
   # - Process data in chunks
   # - Reduce spatial resolution
   # - Use data subsampling
   ```

### Diagnostic Checks
```r
# Data quality checks:
n_pixels_total <- n_rows * n_cols
n_pixels_with_data <- length(valid_pixel_ts)
coverage_ratio <- n_pixels_with_data / n_pixels_total

# Convergence check:
loglik_diff <- abs(loglikelihoods[max_iter] - loglikelihoods[max_iter-1])
is_converged <- loglik_diff / abs(loglikelihoods[max_iter]) < 0.01
```

## Integration and Workflow

### Complete Analysis Pipeline
```r
# Typical workflow:
# 1. Run twoDyeDataAnalysis.R for initial analysis
# 2. Run FullSudy.R for comprehensive model comparison
# 3. Use simulation2k.R for validation with known parameters
```

### Output Compatibility
- **Results saving**: Compatible with RData format
- **Plot formats**: Base R graphics, easily exported
- **Parameter sharing**: Results can be used across scripts

## Extensions

### Adding New Preprocessing Steps
```r
# Example: Add background subtraction
preprocess_with_background <- function(ts, bg_percentile = 10) {
  bg_level <- quantile(ts, bg_percentile/100)
  ts_corrected <- ts - bg_level
  ts_corrected <- pmax(ts_corrected, 0)  # Non-negative
  return(ts_corrected)
}
```

### Custom Lifetime Models
```r
# Example: Non-exponential decay model
# (Would require modification to EM algorithm)
custom_likelihood <- function(tau, params) {
  # Implement custom decay function
}
```

### Batch Processing
```r
# Process multiple files:
file_list <- c("data1.npy", "data2.npy", "data3.npy")
results_list <- list()
for(file in file_list){
  data <- np$load(file)
  results <- analyze_flim_data(data)
  results_list[[file]] <- results
}
```

## References

- **FLIM Analysis**: SwissFluorescence. "Fluorescence Lifetime Imaging Microscopy"
- **Time Series Processing**: Cowpertwait, P. S. P. "Introductory Time Series with R"
- **EM Algorithm**: McLachlan, G. J. "Finite Mixture Models"
- **Model Selection**: Claeskens, G. "Model Selection and Model Averaging"
