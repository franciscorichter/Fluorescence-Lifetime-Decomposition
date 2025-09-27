##############################################
# real_photon_decay_analysis.R
# Processing real photon decay data using EM and per‐pixel local estimation.
# Loads data from a .npy file, pre‐processes the time‐series,
# converts to lifetime events, runs global EM, and reconstructs decomposed images.
##############################################

rm(list = ls())
library(reticulate)

# -------------------------------
# 1) LOAD REAL DATA
# -------------------------------
np <- import("numpy")
# Adjust the file path if necessary:
data <- np$load("best_fit.npy")
# Assume that data dimensions are: Rows x Cols x Time (among other dimensions)
data1 <- data[,,1,,1,1]
cat("Data dimensions (Rows x Cols x Time):", paste(dim(data1), collapse = " x "), "\n")

# -------------------------------
# 2) PREPROCESS THE TIME SERIES
# -------------------------------
# We use a minimum threshold (min_photons) so that only pixels that at some point have 
# at least that many photons are considered. (If you see very few pixels, consider lowering this.)
min_photons <- 3  

# Function: Preprocess a single pixel’s time series:
#   – If the maximum is below the threshold, return an empty vector.
#   – Otherwise, truncate the series at the first occurrence of zero.
preprocess_pixel_ts <- function(ts, min_nonzero = min_photons) {
  # Find the first time where the count is at least min_nonzero.
  start_idx <- which(ts >= min_nonzero)[1]
  if (is.na(start_idx)) return(numeric(0))
  
  # Take the portion of the series starting at start_idx.
  ts_valid <- ts[start_idx:length(ts)]
  
  # Find the peak index (maximum value) in ts_valid.
  peak_idx <- which.max(ts_valid)
  
  # Use the decay phase starting from the peak.
  decay_ts <- ts_valid[peak_idx:length(ts_valid)]
  
  # Enforce non-increasing behavior:
  if (length(decay_ts) > 1) {
    for (t in 2:length(decay_ts)) {
      decay_ts[t] <- min(decay_ts[t - 1], decay_ts[t])
    }
  }
  
  # Truncate at the first occurrence of 0 (if any)
  zero_idx <- which(decay_ts == 0)
  if (length(zero_idx) > 0 && zero_idx[1] > 1) {
    decay_ts <- decay_ts[1:(zero_idx[1] - 1)]
  } else if (length(zero_idx) > 0 && zero_idx[1] == 1) {
    decay_ts <- numeric(0)
  }
  
  return(decay_ts)
}


# Process each pixel
n_rows <- dim(data1)[1]
n_cols <- dim(data1)[2]
n_frames <- dim(data1)[3]

valid_pixel_ts <- list()  # Will store series for pixels with sufficient data
for (i in 1:n_rows) {
  for (j in 1:n_cols) {
    ts <- preprocess_pixel_ts(data1[i, j, ])
    if (length(ts) > 0) {
      valid_pixel_ts[[paste(i, j, sep = "_")]] <- ts
    }
  }
}
cat("Number of pixels scanned:", n_rows * n_cols, "\n")
cat("Number of pixels with data:", length(valid_pixel_ts), "\n")

# Build a 3D array (n_rows x n_cols x n_frames) for the preprocessed data.
# For pixels with no valid data, the entries remain NA.
data1_preproc <- array(NA, dim = c(n_rows, n_cols, n_frames))
for (i in 1:n_rows) {
  for (j in 1:n_cols) {
    pid <- paste(i, j, sep = "_")
    if (!is.null(valid_pixel_ts[[pid]])) {
      ts <- valid_pixel_ts[[pid]]
      data1_preproc[i, j, 1:length(ts)] <- ts
    }
  }
}

# -------------------------------
# 3) CONVERT TIME SERIES TO DISCRETE LIFETIMES
# -------------------------------
# This function works by comparing each frame to the next (after appending a final zero frame)
# and recording at which time a drop (i.e. decay event) occurs.
calculate_lifetimes_from_frames <- function(frames) {
  n_frames <- length(frames)
  n_rows <- nrow(frames[[1]])
  n_cols <- ncol(frames[[1]])
  
  # Append a final zero frame (using zeros; NAs are replaced with 0)
  frames_ext <- frames
  frames_ext[[n_frames + 1]] <- matrix(0, n_rows, n_cols)
  
  lifetimes <- integer(0)
  for (t in seq_len(n_frames)) {
    diff_mat <- frames_ext[[t]] - frames_ext[[t + 1]]
    diff_mat[is.na(diff_mat)] <- 0
    died_counts <- diff_mat[diff_mat > 0]
    if (length(died_counts) > 0) {
      lifetimes <- c(lifetimes, rep(t, sum(died_counts)))
    }
  }
  return(lifetimes)
}

# Build the list of frames (one matrix per time step)
time_series <- lapply(1:n_frames, function(t) data1_preproc[,,t])
global_lifetimes <- calculate_lifetimes_from_frames(time_series)
cat("Total photon decay events (from valid pixels):", length(global_lifetimes), "\n")
if(length(global_lifetimes) == 0){
  stop("No lifetime events detected. Consider adjusting the min_photons threshold or inspecting the data.")
}



##############################################
# real_photon_decay_analysis.R (continued)
##############################################

# -------------------------------
# 3) CONVERT TIME SERIES TO DISCRETE LIFETIMES
# -------------------------------
# Modified function: returns a data frame with columns tau, i, j
calculate_lifetimes_from_frames <- function(frames) {
  n_frames <- length(frames)
  n_rows <- nrow(frames[[1]])
  n_cols <- ncol(frames[[1]])
  
  # Append a final zero frame (replace NAs with 0)
  frames_ext <- frames
  frames_ext[[n_frames + 1]] <- matrix(0, n_rows, n_cols)
  
  lifetimes_df <- data.frame(tau = integer(0), i = integer(0), j = integer(0))
  for (t in seq_len(n_frames)) {
    diff_mat <- frames_ext[[t]] - frames_ext[[t + 1]]
    diff_mat[is.na(diff_mat)] <- 0
    idx <- which(diff_mat > 0, arr.ind = TRUE)
    if (nrow(idx) > 0) {
      for (k in 1:nrow(idx)) {
        count <- diff_mat[idx[k, "row"], idx[k, "col"]]
        if (count > 0) {
          lifetimes_df <- rbind(lifetimes_df,
                                data.frame(tau = rep(t, count),
                                           i = rep(idx[k, "row"], count),
                                           j = rep(idx[k, "col"], count)))
        }
      }
    }
  }
  return(lifetimes_df)
}

# Build list of frames from preprocessed data
frames_list <- lapply(1:n_frames, function(t) data1_preproc[,,t])
df_lifetimes <- calculate_lifetimes_from_frames(frames_list)
cat("Total photon decay events (from valid pixels):", nrow(df_lifetimes), "\n")
if(nrow(df_lifetimes) == 0){
  stop("No lifetime events detected. Consider adjusting the min_photons threshold or inspecting the data.")
}

# -------------------------------
# Prepare Data for EM Analysis
# -------------------------------
df_events <- df_lifetimes  # rename for clarity
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
unique_pixels <- unique(df_events$pixel)
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

# -------------------------------
# EM ALGORITHM FOR k=2 MODEL
# -------------------------------
m <- 2  # k=2 model
N <- nrow(df_events)
N_pix <- length(unique_pixels)
pi_hat <- matrix(1/m, nrow = N_pix, ncol = m)   # initialize mixing proportions
lambda_hat <- runif(m, min = 1, max = 3)          # initialize global decay rates

max_iter <- 30
tau_vals <- df_events$tau
pixel_indices <- df_events$pixel_idx

loglikelihoods <- numeric(max_iter)
lambda_hist <- matrix(0, nrow = max_iter, ncol = m)

for (iter in 1:max_iter) {
  numerators <- matrix(0, nrow = N, ncol = m)
  for (d in 1:m) {
    numerators[, d] <- pi_hat[pixel_indices, d] * lambda_hat[d] * exp(-lambda_hat[d] * tau_vals)
  }
  resp <- numerators / rowSums(numerators)
  
  loglikelihood <- sum(log(rowSums(numerators)))
  loglikelihoods[iter] <- loglikelihood
  
  for (d in 1:m) {
    lambda_hat[d] <- sum(resp[, d]) / sum(resp[, d] * tau_vals)
  }
  
  for (pix in 1:N_pix) {
    coords <- as.numeric(unlist(strsplit(unique_pixels[pix], "_")))
    i0 <- coords[1]
    j0 <- coords[2]
    i_range <- max(1, i0 - 1):min(n_rows, i0 + 1)
    j_range <- max(1, j0 - 1):min(n_cols, j0 + 1)
    neighbors <- expand.grid(i = i_range, j = j_range)
    neighbor_ids <- paste(neighbors$i, neighbors$j, sep = "_")
    idx <- which(df_events$pixel %in% neighbor_ids)
    if (length(idx) > 0) {
      pi_hat[pix, ] <- colMeans(resp[idx, , drop = FALSE])
    }
  }
  pi_hat <- pmax(pi_hat, 1e-6)
  pi_hat <- pi_hat / rowSums(pi_hat)
  
  lambda_hist[iter, ] <- lambda_hat
  cat("k=2, Iteration", iter, ": Log-Likelihood =", round(loglikelihood,2),
      "| Estimated Rates =", paste(round(lambda_hat,3), collapse=", "), "\n")
}

# Composite λ Map for k=2 model
estimated_lambda_map <- matrix(0, nrow = n_rows, ncol = n_cols)
for (k in 1:N_pix) {
  coords <- as.numeric(unlist(strsplit(unique_pixels[k], "_")))
  i_coord <- coords[1]
  j_coord <- coords[2]
  estimated_lambda_map[i_coord, j_coord] <- sum(pi_hat[k, ] * lambda_hat)
}

# BIC for k=2 model
n_obs <- N  # number of photon events
p2 <- N_pix * (m - 1) + m
BIC_2 <- -2 * loglikelihoods[max_iter] + p2 * log(n_obs)
cat("BIC for k=2 model:", BIC_2, "\n")

# -------------------------------
# EM ALGORITHM FOR k=3 MODEL
# -------------------------------
m3 <- 3  # k=3 model
pi_hat_3 <- matrix(1/m3, nrow = N_pix, ncol = m3)
lambda_hat_3 <- runif(m3, min = 1, max = 3)
max_iter <- 30
loglikelihoods_3 <- numeric(max_iter)
lambda_hist_3 <- matrix(0, nrow = max_iter, ncol = m3)

for (iter in 1:max_iter) {
  numerators <- matrix(0, nrow = N, ncol = m3)
  for (d in 1:m3) {
    numerators[, d] <- pi_hat_3[pixel_indices, d] * lambda_hat_3[d] * exp(-lambda_hat_3[d] * tau_vals)
  }
  resp3 <- numerators / rowSums(numerators)
  
  loglikelihood_3 <- sum(log(rowSums(numerators)))
  loglikelihoods_3[iter] <- loglikelihood_3
  
  for (d in 1:m3) {
    lambda_hat_3[d] <- sum(resp3[, d]) / sum(resp3[, d] * tau_vals)
  }
  
  for (pix in 1:N_pix) {
    coords <- as.numeric(unlist(strsplit(unique_pixels[pix], "_")))
    i0 <- coords[1]
    j0 <- coords[2]
    i_range <- max(1, i0 - 1):min(n_rows, i0 + 1)
    j_range <- max(1, j0 - 1):min(n_cols, j0 + 1)
    neighbors <- expand.grid(i = i_range, j = j_range)
    neighbor_ids <- paste(neighbors$i, neighbors$j, sep = "_")
    idx <- which(df_events$pixel %in% neighbor_ids)
    if (length(idx) > 0) {
      pi_hat_3[pix, ] <- colMeans(resp3[idx, , drop = FALSE])
    }
  }
  pi_hat_3 <- pmax(pi_hat_3, 1e-6)
  pi_hat_3 <- pi_hat_3 / rowSums(pi_hat_3)
  
  lambda_hist_3[iter, ] <- lambda_hat_3
  cat("k=3, Iteration", iter, ": Log-Likelihood =", round(loglikelihood_3,2),
      "| Estimated Rates =", paste(round(lambda_hat_3,3), collapse=", "), "\n")
}

# Composite λ Map for k=3 model
estimated_lambda_map_3 <- matrix(0, nrow = n_rows, ncol = n_cols)
for (k in 1:N_pix) {
  coords <- as.numeric(unlist(strsplit(unique_pixels[k], "_")))
  i_coord <- coords[1]
  j_coord <- coords[2]
  estimated_lambda_map_3[i_coord, j_coord] <- sum(pi_hat_3[k, ] * lambda_hat_3)
}

# BIC for k=3 model
p3 <- N_pix * (m3 - 1) + m3
BIC_3 <- -2 * loglikelihoods_3[max_iter] + p3 * log(n_obs) / 9
cat("BIC for k=3 model:", BIC_3, "\n")

# -------------------------------
# RECONSTRUCTION VISUALIZATION
# -------------------------------
# Display the composite lambda maps from the two models side by side.
library(viridis)
par(mfrow = c(1,2))
image(t(estimated_lambda_map[n_rows:1, ]), col = viridis(100),
      main = "Estimated Composite Decay Rates", axes = FALSE)
image(t(estimated_lambda_map_3[n_rows:1, ]), col = viridis(100),
      main = "Estimated Composite $\\lambda$ Map (k=3)", axes = FALSE)
par(mfrow = c(1,1))

# -------------------------------
# REPORT FINAL RESULTS
# -------------------------------
cat("\nFinal Parameter Estimates (k=2):\n")
print(data.frame(Dye = paste("Dye", 1:m),
                 Estimated_Rate = lambda_hat))
cat("\nFinal Parameter Estimates (k=3):\n")
print(data.frame(Dye = paste("Dye", 1:m3),
                 Estimated_Rate = lambda_hat_3))
cat("\nBIC Comparison:\n")
cat("BIC for k=2 model:", BIC_2, "\n")
cat("BIC for k=3 model:", BIC_3, "\n")
















# -------------------------------
# FRAME-BY-FRAME RECONSTRUCTION (k=2)
# -------------------------------
library(viridis)

# (If not already done, compute counts_matrix from df_events)
pixel_counts_df <- aggregate(list(count = df_events$pixel), 
                             by = list(pixel = df_events$pixel), FUN = length)
counts_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
for(i in 1:nrow(pixel_counts_df)){
  parts <- as.integer(unlist(strsplit(pixel_counts_df$pixel[i], "_")))
  counts_matrix[parts[1], parts[2]] <- pixel_counts_df$count[i]
}

# Create a mapping from pixel id to row index in the final estimates matrix.
pixel_index_map <- setNames(seq_along(unique_pixels), unique_pixels)

# Assume final EM estimates (for k=2) are stored in:
final_pi <- pi_hat         # Matrix with rows corresponding to unique_pixels.
final_lambda <- lambda_hat # Vector of length 2.

# Determine global maximum intensity (across all frames and panels) for consistent color scaling.
global_max <- 0
for(t in 1:n_frames){
  for(pid in unique_pixels){
    idx <- pixel_index_map[pid]
    parts <- as.integer(unlist(strsplit(pid, "_")))
    i_idx <- parts[1]
    j_idx <- parts[2]
    n_ij <- counts_matrix[i_idx, j_idx]
    intensity <- n_ij * ( final_pi[idx, 1] * exp(-final_lambda[1] * t) +
                            final_pi[idx, 2] * exp(-final_lambda[2] * t) )
    if(intensity > global_max) global_max <- intensity
  }
}

# For each frame, compute the composite and decomposed reconstructions.
for(t in 1:n_frames){
  composite_mat <- matrix(0, nrow = n_rows, ncol = n_cols)
  dye1_mat <- matrix(0, nrow = n_rows, ncol = n_cols)
  dye2_mat <- matrix(0, nrow = n_rows, ncol = n_cols)
  
  # For each pixel with valid data (in unique_pixels), compute predictions.
  for(pid in unique_pixels){
    idx <- pixel_index_map[pid]
    parts <- as.integer(unlist(strsplit(pid, "_")))
    i_idx <- parts[1]
    j_idx <- parts[2]
    n_ij <- counts_matrix[i_idx, j_idx]
    comp_val <- n_ij * (final_pi[idx,1] * exp(-final_lambda[1] * t) +
                          final_pi[idx,2] * exp(-final_lambda[2] * t))
    d1_val <- n_ij * (final_pi[idx,1] * exp(-final_lambda[1] * t))
    d2_val <- n_ij * (final_pi[idx,2] * exp(-final_lambda[2] * t))
    composite_mat[i_idx, j_idx] <- comp_val
    dye1_mat[i_idx, j_idx] <- d1_val
    dye2_mat[i_idx, j_idx] <- d2_val
  }
  
  # Plot the three images side by side using the same color scale.
  par(mfrow = c(1,3), mar = c(2,2,2,2))
  image(t(composite_mat[n_rows:1, ]), col = viridis(100), zlim = c(0, global_max),
        main = paste("Composite Frame", t), axes = FALSE)
  image(t(dye1_mat[n_rows:1, ]), col = viridis(100), zlim = c(0, global_max),
        main = paste("Dye 1 Frame", t), axes = FALSE)
  image(t(dye2_mat[n_rows:1, ]), col = viridis(100), zlim = c(0, global_max),
        main = paste("Dye 2 Frame", t), axes = FALSE)
  
  # Pause briefly so the plots can be visually inspected.
  Sys.sleep(0.5)
}

