# Usage Examples and Tutorials

This document provides practical examples and tutorials for using the Fluorescence Lifetime Decomposition toolkit, from basic analysis to advanced applications.

## Quick Start Guide

### 1. Basic Model Fitting

```r
# Load required libraries
library(Rcpp)
library(ggplot2)
library(dplyr)
library(viridis)

# Load your lifetime data
load("df_lifetimes.RData")

# Prepare data for analysis
df_events <- df_lifetimes
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
unique_pixels <- unique(df_events$pixel)
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

n_rows <- max(df_events$i)
n_cols <- max(df_events$j)

# Source the C++ EM implementation
Rcpp::sourceCpp("em_algorithm.cpp")

# Define EM wrapper function
em_estimation <- function(df_events, k, n_rows, n_cols, max_iter = 30) {
  df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
  unique_pixels <- unique(df_events$pixel)
  df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

  N_pix <- length(unique_pixels)
  coords <- do.call(rbind, strsplit(unique_pixels, "_"))
  coords <- matrix(as.integer(coords), ncol = 2)
  colnames(coords) <- c("i", "j")

  pi_init <- matrix(1/k, nrow = N_pix, ncol = k)
  lambda_init <- runif(k, min = 1, max = 3)

  res <- EMAlgorithmCpp(df_events$pixel_idx, df_events$tau, pi_init, lambda_init,
                       coords, n_rows, n_cols, max_iter)
  return(res)
}

# Run EM for k=2
k <- 2
result <- em_estimation(df_events, k, n_rows, n_cols)

# Display results
print("Estimated decay rates:")
print(result$lambda_hat)

print("Final log-likelihood:")
print(result$loglikelihoods[length(result$loglikelihoods)])
```

### 2. Model Comparison

```r
# Compare multiple models (k=1 to 5)
results_list <- list()
BIC_values <- numeric(5)

for(k in 1:5) {
  cat("\nFitting model with k =", k, "\n")
  res <- em_estimation(df_events, k, n_rows, n_cols, max_iter = 20)
  results_list[[k]] <- res

  # Calculate BIC
  p_params <- length(unique_pixels) * (k - 1) + k
  BIC <- -2 * res$loglikelihoods[length(res$loglikelihoods)] + p_params * log(nrow(df_events))
  BIC_values[k] <- BIC

  cat("BIC for k =", k, ":", BIC, "\n")
}

# Plot BIC comparison
library(ggplot2)
df_bic <- data.frame(k = 1:5, BIC = BIC_values)
ggplot(df_bic, aes(x = k, y = BIC)) +
  geom_line(color = "red", size = 1.5) +
  geom_point(color = "blue", size = 3) +
  labs(title = "BIC vs. Number of Components",
       x = "Number of Components (k)", y = "BIC") +
  theme_minimal()

# Find best model
best_k <- which.min(BIC_values)
cat("Best model: k =", best_k, "with BIC =", BIC_values[best_k])
```

### 3. Visualization and Reconstruction

```r
# Create reconstruction for best model
best_result <- results_list[[best_k]]
final_pi <- best_result$pi_hat
final_lambda <- best_result$lambda_hat

# Compute photon counts per pixel
pixel_counts_df <- aggregate(list(count = df_events$pixel),
                            by = list(pixel = df_events$pixel), FUN = length)
counts_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
for(i in 1:nrow(pixel_counts_df)){
  parts <- as.integer(unlist(strsplit(pixel_counts_df$pixel[i], "_")))
  counts_matrix[parts[1], parts[2]] <- pixel_counts_df$count[i]
}

# Create mapping
pixel_index_map <- setNames(seq_along(unique_pixels), unique_pixels)

# Reconstruct first few frames
n_frames <- max(df_events$tau)
for(t in 1:min(5, n_frames)) {
  composite_img <- matrix(0, nrow = n_rows, ncol = n_cols)

  for(pid in unique_pixels){
    idx <- pixel_index_map[pid]
    parts <- as.integer(unlist(strsplit(pid, "_")))
    i_idx <- parts[1]
    j_idx <- parts[2]
    n_ij <- counts_matrix[i_idx, j_idx]

    # Predict intensity
    intensity <- n_ij * sum(final_pi[idx, ] * exp(-final_lambda * t))
    composite_img[i_idx, j_idx] <- intensity
  }

  # Plot
  image(t(composite_img[n_rows:1, ]), col = viridis(100),
        main = paste("Frame", t, "Reconstruction"),
        axes = FALSE)
}
```

## Tutorial: Complete Analysis Workflow

### Step 1: Data Loading and Validation

```r
# Load your data
load("df_lifetimes.RData")

# Basic data validation
cat("Dataset summary:")
cat("Total photon events:", nrow(df_lifetimes), "\n")
cat("Grid dimensions:", max(df_lifetimes$i), "x", max(df_lifetimes$j), "\n")
cat("Lifetime range:", min(df_lifetimes$tau), "-", max(df_lifetimes$tau), "\n")

# Check photon distribution
library(ggplot2)
ggplot(df_lifetimes, aes(x = tau)) +
  geom_histogram(binwidth = 1, fill = "blue", alpha = 0.7) +
  labs(title = "Lifetime Distribution", x = "Lifetime", y = "Count") +
  theme_minimal()
```

### Step 2: Model Selection

```r
# Run analysis for multiple k values
source("FullSudy.R")  # This runs the complete analysis

# Or run individual components:
k_values <- 1:5
bic_results <- data.frame(k = k_values, BIC = NA, LogLik = NA)

for(k in k_values) {
  res <- em_estimation(df_events, k, n_rows, n_cols)
  bic_results$BIC[k] <- -2 * res$loglikelihoods[length(res$loglikelihoods)] +
    (length(unique_pixels) * (k-1) + k) * log(nrow(df_events))
  bic_results$LogLik[k] <- res$loglikelihoods[length(res$loglikelihoods)]
}

# Plot comparison
library(tidyr)
bic_long <- gather(bic_results, key = "Metric", value = "Value", -k)
ggplot(bic_long, aes(x = k, y = Value, color = Metric)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  labs(title = "Model Comparison", x = "Number of Components (k)") +
  theme_minimal() +
  facet_wrap(~Metric, scales = "free_y")
```

### Step 3: Parameter Interpretation

```r
# Extract best model results
best_k <- which.min(bic_results$BIC)
best_params <- em_estimation(df_events, best_k, n_rows, n_cols)

# Create spatial maps
pi_hat <- best_params$pi_hat
lambda_hat <- best_params$lambda_hat

# Mixing proportion maps for each component
for(comp in 1:best_k) {
  pi_map <- matrix(0, nrow = n_rows, ncol = n_cols)
  for(i in 1:length(unique_pixels)) {
    coords <- as.integer(unlist(strsplit(unique_pixels[i], "_")))
    pi_map[coords[1], coords[2]] <- pi_hat[i, comp]
  }

  image(t(pi_map[n_rows:1, ]), col = viridis(100),
        main = paste("Mixing Proportion - Component", comp),
        axes = FALSE)
}

# Composite decay rate map
lambda_map <- matrix(0, nrow = n_rows, ncol = n_cols)
for(i in 1:length(unique_pixels)) {
  coords <- as.integer(unlist(strsplit(unique_pixels[i], "_")))
  lambda_map[coords[1], coords[2]] <- sum(pi_hat[i, ] * lambda_hat)
}

image(t(lambda_map[n_rows:1, ]), col = viridis(100),
      main = "Composite Decay Rate Map",
      axes = FALSE)
```

### Step 4: Validation and Diagnostics

```r
# Convergence analysis
loglik_history <- best_params$loglikelihoods
plot(loglik_history, type = "l", col = "blue",
     xlab = "Iteration", ylab = "Log-Likelihood",
     main = "Convergence Plot")
abline(h = tail(loglik_history, 1), lty = 2, col = "red")

# Parameter evolution
lambda_history <- best_params$lambda_hist
matplot(t(lambda_history), type = "l", lty = 1,
        xlab = "Iteration", ylab = "Decay Rate",
        main = "Parameter Evolution")
legend("topright", legend = paste("Component", 1:best_k),
       col = 1:best_k, lty = 1)

# Model fit quality
final_loglik <- tail(loglik_history, 1)
aic <- -2 * final_loglik + 2 * (length(unique_pixels) * (best_k - 1) + best_k)
cat("Model fit metrics:\n")
cat("Final log-likelihood:", final_loglik, "\n")
cat("AIC:", aic, "\n")
cat("BIC:", bic_results$BIC[best_k], "\n")
```

## Advanced Examples

### 1. Custom Initialization

```r
# Use informed initialization based on data
custom_initialization <- function(df_events, k) {
  # Simple initialization based on lifetime quantiles
  tau_sorted <- sort(df_events$tau)
  quantiles <- quantile(tau_sorted, probs = seq(0, 1, length.out = k+1))

  lambda_init <- 1 / diff(quantiles)  # Convert quantiles to rates
  pi_init <- matrix(1/k, nrow = length(unique_pixels), ncol = k)

  return(list(lambda = lambda_init, pi = pi_init))
}

# Use custom initialization
init_params <- custom_initialization(df_events, k)
result_custom <- EMAlgorithmCpp(df_events$pixel_idx, df_events$tau,
                               init_params$pi, init_params$lambda,
                               coords, n_rows, n_cols, max_iter)
```

### 2. Spatial Smoothing Control

```r
# Modify neighborhood size in C++ code (requires recompilation)
# In em_algorithm.cpp, change:
# int i_min = std::max(1, i0 - 1);  // 3x3 neighborhood
# int i_max = std::min(n_rows, i0 + 1);
#
# To 5x5:
# int i_min = std::max(1, i0 - 2);
# int i_max = std::min(n_rows, i0 + 2);
```

### 3. Batch Processing Multiple Files

```r
# Process multiple data files
file_list <- c("data1.RData", "data2.RData", "data3.RData")
results_summary <- data.frame()

for(file in file_list) {
  load(file)

  # Run analysis
  result <- em_estimation(df_events, 2, n_rows, n_cols)

  # Store summary
  summary_row <- data.frame(
    file = file,
    lambda1 = result$lambda_hat[1],
    lambda2 = result$lambda_hat[2],
    final_loglik = result$loglikelihoods[length(result$loglikelihoods)]
  )
  results_summary <- rbind(results_summary, summary_row)
}

print(results_summary)
```

### 4. Bootstrap Uncertainty Estimation

```r
# Bootstrap analysis for parameter uncertainty
n_bootstrap <- 50
bootstrap_results <- matrix(NA, nrow = n_bootstrap, ncol = 2)

for(b in 1:n_bootstrap) {
  # Resample photons with replacement
  sample_idx <- sample(1:nrow(df_events), replace = TRUE)
  df_bootstrap <- df_events[sample_idx, ]

  # Fit model
  res_bootstrap <- em_estimation(df_bootstrap, 2, n_rows, n_cols, max_iter = 15)

  # Store results
  bootstrap_results[b, ] <- res_bootstrap$lambda_hat
}

# Calculate confidence intervals
ci_lambda1 <- quantile(bootstrap_results[,1], c(0.025, 0.975))
ci_lambda2 <- quantile(bootstrap_results[,2], c(0.025, 0.975))

cat("95% CI for lambda1:", ci_lambda1, "\n")
cat("95% CI for lambda2:", ci_lambda2, "\n")
```

### 5. Integration with Experimental Design

```r
# Analyze multiple experimental conditions
conditions <- c("control", "treatment1", "treatment2")
condition_results <- list()

for(condition in conditions) {
  # Load condition-specific data
  load(paste0(condition, "_lifetimes.RData"))

  # Analyze
  result <- em_estimation(df_events, 2, n_rows, n_cols)
  condition_results[[condition]] <- result
}

# Compare conditions
comparison_df <- data.frame()
for(condition in conditions) {
  res <- condition_results[[condition]]
  row <- data.frame(
    condition = condition,
    lambda1 = res$lambda_hat[1],
    lambda2 = res$lambda_hat[2]
  )
  comparison_df <- rbind(comparison_df, row)
}

print(comparison_df)
```

## Troubleshooting Common Issues

### Issue 1: Poor Convergence

```r
# Diagnostic function
diagnose_convergence <- function(result, threshold = 0.01) {
  loglik <- result$loglikelihoods
  final_ll <- tail(loglik, 1)
  prev_ll <- tail(loglik, 2)[1]

  rel_change <- abs(final_ll - prev_ll) / abs(final_ll)

  cat("Convergence check:\n")
  cat("Final log-likelihood:", final_ll, "\n")
  cat("Relative change:", rel_change, "\n")
  cat("Converged:", rel_change < threshold, "\n")

  return(rel_change < threshold)
}

# Usage
converged <- diagnose_convergence(result)
if(!converged) {
  # Try with more iterations
  result <- em_estimation(df_events, k, n_rows, n_cols, max_iter = 50)
}
```

### Issue 2: Unrealistic Parameter Estimates

```r
# Parameter validation
validate_parameters <- function(lambda_hat, pi_hat) {
  issues <- c()

  # Check decay rates
  if(any(lambda_hat < 0.1) || any(lambda_hat > 10)) {
    issues <- c(issues, "Decay rates outside typical range [0.1, 10]")
  }

  # Check mixing proportions
  if(any(pi_hat < 0) || any(pi_hat > 1)) {
    issues <- c(issues, "Mixing proportions outside [0,1] range")
  }

  # Check normalization
  row_sums <- rowSums(pi_hat)
  if(any(abs(row_sums - 1) > 0.01)) {
    issues <- c(issues, "Mixing proportions not properly normalized")
  }

  return(issues)
}

# Usage
problems <- validate_parameters(result$lambda_hat, result$pi_hat)
if(length(problems) > 0) {
  cat("Parameter issues found:\n")
  print(problems)
}
```

### Issue 3: Memory Issues with Large Datasets

```r
# Memory-efficient processing
process_large_dataset <- function(df_events, chunk_size = 10000) {
  # Split data into chunks
  n_chunks <- ceiling(nrow(df_events) / chunk_size)
  chunk_results <- list()

  for(chunk in 1:n_chunks) {
    start_idx <- (chunk - 1) * chunk_size + 1
    end_idx <- min(chunk * chunk_size, nrow(df_events))
    df_chunk <- df_events[start_idx:end_idx, ]

    # Process chunk
    res_chunk <- em_estimation(df_chunk, 2, n_rows, n_cols, max_iter = 10)
    chunk_results[[chunk]] <- res_chunk
  }

  # Combine results (simplified - would need proper aggregation)
  return(chunk_results)
}

# Usage for large datasets
if(nrow(df_events) > 50000) {
  results <- process_large_dataset(df_events)
}
```

## Performance Optimization Tips

### 1. Compilation and Caching

```r
# Pre-compile C++ code for faster subsequent runs
Rcpp::sourceCpp("em_algorithm.cpp")

# Cache compiled function
em_func <- EMAlgorithmCpp
```

### 2. Parallel Processing (where applicable)

```r
# For multiple k values (requires parallel package)
library(parallel)
library(doParallel)

# Setup parallel backend
cl <- makeCluster(detectCores() - 1)
registerDoParallel(cl)

# Parallel model fitting
k_values <- 1:5
results_parallel <- foreach(k = k_values) %dopar% {
  em_estimation(df_events, k, n_rows, n_cols)
}

# Stop cluster
stopCluster(cl)
```

### 3. Memory Management

```r
# Clear large intermediate objects
gc()  # Garbage collection

# Remove temporary variables
rm(large_data_frame)
gc()

# For very large datasets, consider data.table
# library(data.table)
# df_events <- as.data.table(df_events)
```

## Output and Reporting

### 1. Generate Analysis Report

```r
# Create comprehensive report
generate_report <- function(results_list, bic_values, best_k) {
  cat("=== Fluorescence Lifetime Analysis Report ===\n\n")

  cat("1. Model Comparison:\n")
  for(k in 1:length(bic_values)) {
    cat("k =", k, ": BIC =", round(bic_values[k], 2), "\n")
  }
  cat("Selected model: k =", best_k, "\n\n")

  cat("2. Parameter Estimates:\n")
  best_result <- results_list[[best_k]]
  for(i in 1:best_k) {
    cat("Component", i, ": λ =", round(best_result$lambda_hat[i], 3), "\n")
  }

  cat("\n3. Convergence Quality:\n")
  final_ll <- best_result$loglikelihoods[length(best_result$loglikelihoods)]
  cat("Final log-likelihood:", round(final_ll, 2), "\n")
}

# Usage
generate_report(results_list, BIC_values, best_k)
```

### 2. Export Results

```r
# Save complete results
save(results_list, BIC_values, df_events,
     file = paste0("analysis_results_", format(Sys.Date(), "%Y%m%d"), ".RData"))

# Export parameter summary
params_summary <- data.frame(
  Component = 1:best_k,
  Decay_Rate = best_result$lambda_hat,
  Final_LogLik = final_ll
)
write.csv(params_summary, "parameter_estimates.csv", row.names = FALSE)
```

## Getting Help

### Common Questions

1. **"My model isn't converging"**
   - Check data quality and photon counts
   - Try different initialization
   - Increase max_iter

2. **"Results don't make physical sense"**
   - Validate input data format
   - Check parameter ranges
   - Compare with simulation results

3. **"The analysis is too slow"**
   - Use C++ implementation
   - Reduce max_iter for initial analysis
   - Consider data subsampling

### Example Gallery

The `docs/` folder contains example outputs showing:
- Model comparison plots
- Spatial parameter maps
- Convergence diagnostics
- Reconstruction animations

### Further Reading

- **EM Algorithm Theory**: Bishop, C. M. "Pattern Recognition and Machine Learning"
- **Fluorescence Lifetime Imaging**: Becker, W. "Fluorescence Lifetime Imaging Techniques and Applications"
- **Rcpp Performance**: Eddelbuettel, D. "Seamless R and C++ Integration with Rcpp"
