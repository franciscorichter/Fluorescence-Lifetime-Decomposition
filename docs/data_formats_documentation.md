# Data Formats and Preprocessing Documentation

## Overview
This document describes the data formats used in the fluorescence lifetime decomposition toolkit, including input requirements, preprocessing steps, and output formats.

## Input Data Formats

### 1. Lifetime Events Format (df_lifetimes)

**Primary input format for EM analysis**

**Structure**: R data frame with columns:
- `tau`: Numeric, lifetime values (discrete time units)
- `i`: Integer, row coordinates (1-based indexing)
- `j`: Integer, column coordinates (1-based indexing)

**Example**:
```r
# Load example data
load("df_lifetimes.RData")

# Structure
head(df_lifetimes)
#   tau i j
# 1   1 1 1
# 2   2 1 1
# 3   1 1 2
# 4   3 1 2
# ...

# Summary statistics
summary(df_lifetimes)
#      tau           i              j
# Min. : 1   Min. :  1.0   Min. :  1.0
# Max. :20   Max. :100.0   Max. :100.0
# Mean : 8.5 Mean : 50.5   Mean : 50.5
```

**Generation Methods**:
- From raw time series (see `twoDyeDataAnalysis.R`)
- From simulation (see `simulation2k.R`)
- From experimental FLIM data preprocessing

### 2. Raw Time Series Format

**Input format for preprocessing pipeline**

**Structure**: NumPy .npy file with dimensions typically:
- `data[rows, cols, time]` - Basic format
- `data[rows, cols, time, channels, z_slices, ...]` - Extended format

**Example**:
```r
# Load in R (requires reticulate)
library(reticulate)
np <- import("numpy")
data <- np$load("best_fit.npy")

# Extract relevant dimensions
data1 <- data[,,1,,1,1]  # rows × cols × time
dim(data1)  # e.g., 100 100 50

# Convert to R array
time_series <- array(data1, dim = dim(data1))
```

**Data Characteristics**:
- **Data type**: Integer photon counts
- **Range**: 0 to thousands of photons per pixel-frame
- **Temporal structure**: Decay pattern over time frames

### 3. Simulation Data Format

**Generated internally for validation**

**Structure**: List with components:
- `masks`: Binary matrices defining dye regions
- `n_photons`: Photon count arrays [dye, row, col]
- `df_photons`: Data frame with photon lifetimes and positions

## Data Preprocessing Pipeline

### 1. Raw Time Series Preprocessing

#### Purpose
Convert raw photon count time series to analyzable format by:
- Removing noise and background
- Identifying valid decay phases
- Ensuring physical constraints

#### Algorithm Steps
```r
preprocess_pixel_ts <- function(ts, min_nonzero = 3) {
  # 1. Find first valid starting point
  start_idx <- which(ts >= min_nonzero)[1]

  # 2. Extract valid portion
  ts_valid <- ts[start_idx:length(ts)]

  # 3. Find intensity peak
  peak_idx <- which.max(ts_valid)

  # 4. Take decay phase
  decay_ts <- ts_valid[peak_idx:length(ts_valid)]

  # 5. Enforce monotonicity
  for(t in 2:length(decay_ts)) {
    decay_ts[t] <- min(decay_ts[t-1], decay_ts[t])
  }

  # 6. Truncate at zero
  zero_idx <- which(decay_ts == 0)[1]
  if(!is.na(zero_idx) && zero_idx > 1) {
    decay_ts <- decay_ts[1:(zero_idx-1)]
  }

  return(decay_ts)
}
```

#### Parameters
- `min_nonzero`: Minimum photon threshold (default: 3)
- `ts`: Input time series vector

#### Output
- Cleaned time series or empty vector if insufficient data

### 2. Lifetime Event Conversion

#### Purpose
Convert preprocessed time series to discrete lifetime events for EM analysis

#### Algorithm
```r
calculate_lifetimes_from_frames <- function(frames) {
  # Append zero frame
  frames_ext <- c(frames, list(matrix(0, nrow(frames[[1]]), ncol(frames[[1]]))))

  # Detect decay events
  lifetimes_df <- data.frame(tau = integer(0), i = integer(0), j = integer(0))
  for(t in 1:length(frames)) {
    diff_mat <- frames_ext[[t]] - frames_ext[[t+1]]
    idx <- which(diff_mat > 0, arr.ind = TRUE)
    for(k in 1:nrow(idx)) {
      count <- diff_mat[idx[k,"row"], idx[k,"col"]]
      lifetimes_df <- rbind(lifetimes_df,
        data.frame(tau = rep(t, count),
                  i = rep(idx[k,"row"], count),
                  j = rep(idx[k,"col"], count)))
    }
  }
  return(lifetimes_df)
}
```

#### Input Format
- `frames`: List of intensity matrices (one per time frame)
- Each matrix: Numeric, dimensions (rows × cols)

#### Output Format
- Data frame with decay events
- Each row represents one photon decay event

## Data Quality Requirements

### Minimum Requirements
- **Photon coverage**: >5 photons per pixel for reliable estimation
- **Lifetime range**: Typically 1-50 time units
- **Spatial coverage**: >50% of pixels should have valid data

### Quality Metrics
```r
# Calculate quality metrics
total_photons <- nrow(df_lifetimes)
pixels_with_data <- length(unique(paste(df_lifetimes$i, df_lifetimes$j)))
grid_size <- n_rows * n_cols
coverage_ratio <- pixels_with_data / grid_size

# Lifetime distribution
tau_range <- range(df_lifetimes$tau)
tau_mean <- mean(df_lifetimes$tau)
tau_cv <- sd(df_lifetimes$tau) / tau_mean  # Coefficient of variation
```

### Data Validation Checks
```r
# Essential checks
check_data_quality <- function(df_lifetimes, n_rows, n_cols) {
  # 1. Check dimensions
  if(max(df_lifetimes$i) > n_rows || max(df_lifetimes$j) > n_cols) {
    warning("Lifetime coordinates exceed grid dimensions")
  }

  # 2. Check photon coverage
  pixel_counts <- table(paste(df_lifetimes$i, df_lifetimes$j))
  if(mean(pixel_counts) < 10) {
    warning("Low average photon count per pixel")
  }

  # 3. Check lifetime range
  tau_range <- range(df_lifetimes$tau)
  if(tau_range[2] - tau_range[1] < 5) {
    warning("Limited lifetime range may affect parameter estimation")
  }

  return(list(coverage = length(pixel_counts)/(n_rows*n_cols),
              mean_photons = mean(pixel_counts),
              tau_range = tau_range))
}
```

## Output Data Formats

### 1. EM Results Format

**Structure**: R list with components:
```r
em_results <- list(
  pi_hat = matrix(),           # Mixing proportions (N_pix × k)
  lambda_hat = numeric(),      # Decay rates (length k)
  loglikelihoods = numeric(),  # Log-likelihood history
  lambda_hist = matrix()       # Parameter evolution
)
```

**Interpretation**:
- `pi_hat[i,j]`: Proportion of dye j at pixel i
- `lambda_hat[j]`: Decay rate of dye j
- `loglikelihoods[t]`: Model fit quality at iteration t

### 2. Reconstruction Results

**Structure**: List with frame-by-frame predictions:
```r
recon_results <- list(
  frames = list(),      # List of prediction matrices per frame
  global_max = numeric() # Maximum intensity for scaling
)

# Each frame contains:
frame_data <- list(
  composite = matrix(),      # Total predicted intensity
  components = matrix()      # Individual dye contributions
)
```

### 3. Model Comparison Results

**Structure**: Data frames for analysis:
```r
# BIC comparison
bic_df <- data.frame(k = 1:5, BIC = bic_values)

# Parameter estimates
params_df <- data.frame(
  Model = paste("k=", 1:5),
  LogLik = final_logliks,
  BIC = bic_values,
  Parameters = param_counts
)
```

## Data Storage and Management

### Saving Results
```r
# Save EM results
save(results_list, BIC_values, recon_results,
     file = "em_results_multi_k.RData")

# Save processed data
save(df_lifetimes, file = "processed_lifetimes.RData")
```

### Loading Results
```r
# Load for further analysis
load("em_results_multi_k.RData")
load("df_lifetimes.RData")
```

## Data Format Conversion

### Time Series to Lifetime Events
```r
# Convert preprocessed time series to lifetime format
frames_list <- lapply(1:n_frames, function(t) data1_preproc[,,t])
df_lifetimes <- calculate_lifetimes_from_frames(frames_list)
```

### Lifetime Events to Analysis Format
```r
# Prepare for EM analysis
df_events <- df_lifetimes
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
unique_pixels <- unique(df_events$pixel)
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))
```

## Example Datasets

### 1. Synthetic Two-Dye Dataset
```r
# Generated by simulation2k.R
# Characteristics:
# - 50×50 pixel grid
# - Two dye regions with known λ = [0.5, 1.5]
# - ~25,000 photon events
# - Ground truth for validation
```

### 2. Experimental FLIM Dataset
```r
# Typical experimental data
# Characteristics:
# - 100×100 pixel grid
# - 50 time frames
# - Multiple fluorophores
# - Real photon statistics
```

## Data Visualization Formats

### Spatial Maps
- **Matrix format**: (rows × cols) for 2D visualization
- **Color mapping**: Viridis scale for intensity representation
- **Scaling**: Consistent color limits across frames

### Temporal Data
- **Animation frames**: List of matrices for gganimate
- **Time evolution**: Frame index as animation variable

### Statistical Plots
- **Data frames**: Compatible with ggplot2
- **Long format**: For faceted visualizations

## Performance Considerations

### Memory Usage by Data Type
- **Raw time series**: High (rows × cols × time × bytes_per_count)
- **Lifetime events**: Medium (n_events × 3 × bytes_per_value)
- **EM parameters**: Low (n_pixels × k × bytes_per_parameter)

### Processing Time by Stage
- **Preprocessing**: O(rows × cols × time)
- **Lifetime conversion**: O(rows × cols × time)
- **EM fitting**: O(n_events × k × iterations)
- **Visualization**: O(rows × cols × frames × k)

## Troubleshooting Data Issues

### Common Data Problems

1. **Insufficient Photon Counts**
   - Symptom: Poor EM convergence, high BIC values
   - Solution: Lower `min_nonzero` threshold or increase data collection time

2. **Coordinate Out of Bounds**
   - Symptom: Index errors in processing
   - Solution: Verify data dimensions and coordinate ranges

3. **Non-monotonic Time Series**
   - Symptom: Negative lifetimes or invalid decay patterns
   - Solution: Check preprocessing monotonicity enforcement

4. **Memory Issues**
   - Symptom: Out of memory errors
   - Solution: Process in chunks or subsample data

### Data Validation Checklist
- [ ] Coordinate ranges match grid dimensions
- [ ] Lifetime values are positive integers
- [ ] Sufficient photon coverage per pixel
- [ ] Reasonable lifetime distribution
- [ ] No obvious data corruption or artifacts

## References

- **FLIM Data Formats**: FLI standard formats and conventions
- **Data Preprocessing**: Image processing techniques for time series
- **Quality Metrics**: Statistical measures for photon counting data
