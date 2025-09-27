##############################################
# real_photon_decay_analysis_multi_k.R
# Analysis using C++ EM Estimation for k = 1,2,3,4,5.
# Loads saved lifetime events, runs EM via a C++ function for each k,
# reconstructs frame-by-frame decompositions, prints diagnostics and BIC values,
# and produces summary plots.
##############################################

rm(list = ls())
library(reticulate)
library(Rcpp)
library(viridis)
library(ggplot2)
library(dplyr)
library(tidyr)

# -------------------------------
# 0) Load preprocessed lifetime data (df_lifetimes)
# -------------------------------
load("df_lifetimes.RData")  # loads df_lifetimes (columns: tau, i, j)
cat("Loaded df_lifetimes with", nrow(df_lifetimes), "events.\n")

# -------------------------------
# 1) Prepare Data for EM Analysis
# -------------------------------
df_events <- df_lifetimes
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
unique_pixels <- unique(df_events$pixel)
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))
n_rows <- max(df_events$i)
n_cols <- max(df_events$j)
tau_vals <- df_events$tau
pixel_indices <- df_events$pixel_idx

# -------------------------------
# 2) Source the C++ EM Function
# -------------------------------
Rcpp::sourceCpp("em_algorithm.cpp")
# The C++ function EMAlgorithmCpp is now available.

# -------------------------------
# 3) Wrapper Function for EM Estimation
# -------------------------------
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

  res <- EMAlgorithmCpp(df_events$pixel_idx, df_events$tau, pi_init, lambda_init, coords, n_rows, n_cols, max_iter)
  return(res)
}

# -------------------------------
# 4) Run EM Estimation for k = 1,2,3,4,5 and store results
# -------------------------------
results_list <- list()
BIC_values <- numeric()
for(k in 1:5){
  cat("\n******** Running EM for k =", k, "********\n")
  res <- em_estimation(df_events, k = k, n_rows, n_cols, max_iter = 10)
  results_list[[as.character(k)]] <- res
  # Number of free parameters: each pixel has (k-1) free mixing params, plus k global rates
  p_params <- length(unique_pixels) * (k - 1) + k
  BIC <- -2 * res$loglikelihoods[length(res$loglikelihoods)] + p_params * log(nrow(df_events)) /9
  BIC_values[k] <- BIC
  cat("Final estimated decay rates (k =", k, "):\n")
  print(res$lambda_hat)
  cat("BIC for k =", k, "model:", BIC, "\n")
}

# Plot BIC vs. k
df_BIC <- data.frame(k = 1:5, BIC = BIC_values)
ggplot(df_BIC, aes(x = k, y = BIC)) +
  geom_line(col = "red", size = 1.5) +
  geom_point(col = "blue", size = 3) +
  labs(title = "BIC vs. Number of Components", x = "Number of Components (k)", y = "BIC") +
  theme_minimal()

# -------------------------------
# 5) Reconstruction: Frame-by-Frame Decomposition for each k
# -------------------------------
# We'll do reconstruction for each k and save the plots.
# The reconstruction formula (for each pixel at frame t) is:
# I_pred(t; i,j) = n_{ij} * sum_{d=1}^k [pi_d(i,j) * exp(-lambda_d * t)]
# Additionally, we plot the individual contributions for each component.
# Compute photon counts per pixel from df_events:
pixel_counts_df <- aggregate(list(count = df_events$pixel),
                             by = list(pixel = df_events$pixel),
                             FUN = length)
counts_matrix <- matrix(0, nrow = n_rows, ncol = n_cols)
for(i in 1:nrow(pixel_counts_df)){
  parts <- as.integer(unlist(strsplit(pixel_counts_df$pixel[i], "_")))
  counts_matrix[parts[1], parts[2]] <- pixel_counts_df$count[i]
}

# Create mapping from unique_pixels to indices
pixel_index_map <- setNames(seq_along(unique_pixels), unique_pixels)

# Determine global maximum intensity for fixed color scale across frames
# (We do this for each k separately.)
recon_results <- list()
for(k in 1:5){
  cat("\n--- Reconstruction for k =", k, "---\n")
  res <- results_list[[as.character(k)]]
  final_pi <- res$pi_hat
  final_lambda <- res$lambda_hat
  # Compute global max intensity across all frames
  global_max <- 0
  n_frames <- max(df_events$tau)
  for(t in 1:n_frames){
    for(pid in unique_pixels){
      idx <- pixel_index_map[pid]
      parts <- as.integer(unlist(strsplit(pid, "_")))
      i_idx <- parts[1]
      j_idx <- parts[2]
      n_ij <- counts_matrix[i_idx, j_idx]
      intensity <- n_ij * sum(final_pi[idx, ] * exp(-final_lambda * t))
      global_max <- max(global_max, intensity)
    }
  }
  cat("Global maximum intensity (k =", k, "):", global_max, "\n")

  # For each frame, compute the composite and individual component images.
  frames_list_recon <- list()
  for(t in 1:n_frames){
    composite_img <- matrix(0, nrow = n_rows, ncol = n_cols)
    comp_components <- matrix(0, nrow = n_rows, ncol = n_cols * k)
    # We'll arrange the k component images side-by-side after the composite.
    for(pid in unique_pixels){
      idx <- pixel_index_map[pid]
      parts <- as.integer(unlist(strsplit(pid, "_")))
      i_idx <- parts[1]
      j_idx <- parts[2]
      n_ij <- counts_matrix[i_idx, j_idx]
      comp_val <- n_ij * sum(final_pi[idx, ] * exp(-final_lambda * t))
      composite_img[i_idx, j_idx] <- comp_val
      for(d in 1:k){
        comp_components[i_idx, (d-1)*n_cols + j_idx] <- n_ij * final_pi[idx, d] * exp(-final_lambda[d] * t)
      }
    }
    # Combine composite and component images side by side
    full_img <- cbind(composite_img, comp_components)
    frames_list_recon[[t]] <- full_img
  }
  recon_results[[as.character(k)]] <- list(frames = frames_list_recon, global_max = global_max)

  # Plot each frame for this k
  for(t in 1:3){
    par(mfrow = c(1,1), mar = c(2,2,2,2))
    image(t(frames_list_recon[[t]][n_rows:1, ]), col = viridis(100),
          zlim = c(0, recon_results[[as.character(k)]]$global_max),
          main = paste("Reconstruction for k =", k, "Frame", t), axes = FALSE)
    Sys.sleep(0.5)
  }
}


# Optionally, save the BIC summary and reconstruction summaries for later use:
save(results_list, BIC_values, recon_results, file = "em_results_multi_k.RData")



# Plot Final Log-Likelihood vs. Number of Components
df_lik <- data.frame(k = 1:5, FinalLogLik = sapply(results_list, function(x) x$loglikelihoods[x$loglikelihoods]))
ggplot(df_lik, aes(x = k, y = FinalLogLik)) +
  geom_line(color = "darkgreen", size = 1.2) +
  geom_point(color = "orange", size = 3) +
  labs(title = "Final Log-Likelihood vs. Number of Components", 
       x = "Number of Components (k)", 
       y = "Final Log-Likelihood") +
  theme_minimal()

