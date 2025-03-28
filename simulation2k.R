# ----------------------------
# Full R Script: Data Simulation, EM Estimation with Neighborhood Update,
# Animation of Photon Counts, and Plots of Likelihood & Rates (k = 2)
# ----------------------------

# Load required libraries
library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(gganimate)

# Set seed for reproducibility
set.seed(42)

# ----------------------------
# Simulation Parameters
# ----------------------------
H <- 50       # number of rows
W <- 50       # number of columns
m <- 2        # number of dyes (k = 2)
lambda_true <- c(0.5, 1.5)  # true decay rates for the two dyes

# ----------------------------
# Define Spatial Masks
# ----------------------------
grid <- expand.grid(i = 1:H, j = 1:W)

# Dye 1: Circle centered at (25,25) with radius 15
mask1 <- with(grid, ifelse(sqrt((i - 25)^2 + (j - 25)^2) <= 15, 1, 0))
mask1 <- matrix(mask1, nrow = H, ncol = W, byrow = FALSE)

# Dye 2: Rectangle covering rows 20 to 40 and columns 10 to 30
mask2 <- matrix(0, nrow = H, ncol = W)
mask2[20:40, 10:30] <- 1

masks <- list(mask1, mask2)

# ----------------------------
# Simulate Photon Counts & Lifetimes
# ----------------------------
photons_list <- list()
n_photons <- array(0, dim = c(m, H, W))

for (d in 1:m) {
  for (i in 1:H) {
    for (j in 1:W) {
      if (masks[[d]][i, j] == 1) {
        n_ij <- sample(5:30, 1)  # photon count uniformly between 5 and 30
        n_photons[d, i, j] <- n_ij
        lifetimes <- rexp(n_ij, rate = lambda_true[d])
        df_temp <- data.frame(tau = lifetimes, i = i, j = j, true_d = d)
        photons_list[[length(photons_list) + 1]] <- df_temp
      }
    }
  }
}
df_photons <- do.call(rbind, photons_list)
df_photons <- df_photons[sample(nrow(df_photons)), ]  # shuffle rows

# Create a pixel identifier (as "i_j")
df_photons <- df_photons %>% mutate(pixel = paste(i, j, sep = "_"))
unique_pixels <- unique(df_photons$pixel)
df_photons$pixel_idx <- as.integer(factor(df_photons$pixel, levels = unique_pixels))

# ----------------------------
# Animate the Photon Counts Over Time
# ----------------------------
# We define "frames" as discrete time points.
# At each time t, a photon is counted if its lifetime (tau) is greater than t.
T_frames <- 10  # number of frames
grid_pixels <- expand.grid(i = 1:H, j = 1:W)
intensity_list <- list()
for(t in 0:(T_frames - 1)){
  df_temp <- df_photons %>%
    filter(tau > t) %>%
    group_by(i, j) %>%
    summarise(intensity = n(), .groups = "drop")
  # Join with the full grid to include pixels with zero counts
  df_temp <- full_join(grid_pixels, df_temp, by = c("i", "j"))
  df_temp$intensity[is.na(df_temp$intensity)] <- 0
  df_temp$frame <- t
  intensity_list[[length(intensity_list) + 1]] <- df_temp
}
df_intensity <- bind_rows(intensity_list)
# Set fixed limits based on time 0 (maximum intensity)
max_intensity <- max(df_intensity$intensity[df_intensity$frame == 0])

# Create the animated plot using gganimate
p_anim <- ggplot(df_intensity, aes(x = j, y = i, fill = intensity)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, max_intensity)) +
  labs(title = 'Photon Counts per Pixel at time: {frame_time}', x = "Column", y = "Row") +
  theme_minimal() +
  scale_y_reverse() +
  transition_manual(frame) +
  ease_aes('linear')

# To view the animation (if running in an interactive session), use:
animate(p_anim, nframes = T_frames, fps = 2, width = 400, height = 400)

# ----------------------------
# EM Algorithm Initialization
# ----------------------------
N <- nrow(df_photons)
N_pix <- length(unique_pixels)
pi_hat <- matrix(1/m, nrow = N_pix, ncol = m)  # pixel-specific mixing proportions
lambda_hat <- runif(m, min = 1, max = 3)         # global decay rates

max_iter <- 30
tau_vals <- df_photons$tau
pixel_indices <- df_photons$pixel_idx

# Prepare storage for log-likelihood and decay rate history
loglikelihoods <- numeric(max_iter)
lambda_hist <- matrix(0, nrow = max_iter, ncol = m)

# ----------------------------
# EM Algorithm Iterations with Neighborhood Update and Printing
# ----------------------------
for (iter in 1:max_iter) {
  # E-step: Compute responsibilities for each photon
  numerators <- matrix(0, nrow = N, ncol = m)
  for (d in 1:m) {
    numerators[, d] <- pi_hat[pixel_indices, d] * lambda_hat[d] * exp(-lambda_hat[d] * tau_vals)
  }
  resp <- numerators / rowSums(numerators)

  # Compute log-likelihood for this iteration
  loglikelihood <- sum(log(rowSums(numerators)))
  loglikelihoods[iter] <- loglikelihood

  # M-step: Update global decay rates
  for (d in 1:m) {
    lambda_hat[d] <- sum(resp[, d]) / sum(resp[, d] * tau_vals)
  }

  # M-step: Update pixel-specific mixing proportions using 3x3 neighborhood
  for (pix in 1:N_pix) {
    coords <- as.numeric(unlist(strsplit(unique_pixels[pix], "_")))
    i0 <- coords[1]
    j0 <- coords[2]

    # Determine neighbors (3x3 window, taking grid boundaries into account)
    i_range <- max(1, i0 - 1):min(H, i0 + 1)
    j_range <- max(1, j0 - 1):min(W, j0 + 1)
    neighbors <- expand.grid(i = i_range, j = j_range)
    neighbor_ids <- paste(neighbors$i, neighbors$j, sep = "_")

    idx <- which(df_photons$pixel %in% neighbor_ids)
    if (length(idx) > 0) {
      pi_hat[pix, ] <- colMeans(resp[idx, , drop = FALSE])
    }
  }
  # Normalize for numerical stability
  pi_hat <- pmax(pi_hat, 1e-6)
  pi_hat <- pi_hat / rowSums(pi_hat)

  # Store current lambda estimates
  lambda_hist[iter, ] <- lambda_hat

  # Print progress for this iteration
  cat("Iteration", iter,
      ": Log-Likelihood =", round(loglikelihood, 2),
      "| Estimated Rates =", paste(round(lambda_hat, 3), collapse = ", "), "\n")
}

# ----------------------------
# Final Comparison of Decay Rates
# ----------------------------
rates_comparison <- data.frame(Dye = paste("Dye", 1:m),
                               True_Rate = lambda_true,
                               Estimated_Rate = lambda_hat,
                               Absolute_Error = abs(lambda_true - lambda_hat))
cat("\nFinal Comparison of Decay Rates:\n")
print(rates_comparison)

# ----------------------------
# Composite λ Maps (Weighted by Mixing Proportions)
# ----------------------------
# Estimated composite λ map
estimated_lambda_map <- matrix(0, nrow = H, ncol = W)
for (k in 1:N_pix) {
  coords <- as.numeric(unlist(strsplit(unique_pixels[k], "_")))
  i_coord <- coords[1]
  j_coord <- coords[2]
  estimated_lambda_map[i_coord, j_coord] <- sum(pi_hat[k, ] * lambda_hat)
}

# True composite λ map (weighted by mask activation)
true_lambda_map <- matrix(0, nrow = H, ncol = W)
for (d in 1:m) {
  true_lambda_map <- true_lambda_map + masks[[d]] * lambda_true[d]
}

# Plot composite λ maps using base R image
par(mfrow = c(1, 2))
image(t(true_lambda_map[nrow(true_lambda_map):1, ]), col = viridis(100),
      main = "True Composite λ Map", axes = FALSE)
image(t(estimated_lambda_map[nrow(estimated_lambda_map):1, ]), col = viridis(100),
      main = "Estimated Composite λ Map", axes = FALSE)
par(mfrow = c(1, 1))

# ----------------------------
# Decomposition: Compare Mixing Proportions (π)
# ----------------------------
# Compute true mixing proportions from photon counts
true_pi_grid <- array(0, dim = c(m, H, W))
total_counts <- apply(n_photons, c(2, 3), sum)
for (d in 1:m) {
  true_pi_grid[d, , ] <- ifelse(total_counts > 0, n_photons[d, , ] / total_counts, 0)
}

# Map estimated mixing proportions back to grid
est_pi_grid <- array(0, dim = c(m, H, W))
for (k in 1:N_pix) {
  coords <- as.numeric(unlist(strsplit(unique_pixels[k], "_")))
  i_coord <- coords[1]
  j_coord <- coords[2]
  est_pi_grid[, i_coord, j_coord] <- pi_hat[k, ]
}

# Compute Mean Absolute Error (MAE) per dye (only for pixels with photons)
mae_per_dye <- numeric(m)
for (d in 1:m) {
  mask <- total_counts > 0
  err <- abs(true_pi_grid[d, , ][mask] - est_pi_grid[d, , ][mask])
  mae_per_dye[d] <- mean(err)
}
decomp_results <- data.frame(Dye = paste("Dye", 1:m), MAE = mae_per_dye)
cat("\nMixing Proportions Decomposition Error:\n")
print(decomp_results)

# ----------------------------
# Plot: Log-Likelihood over Iterations
# ----------------------------
df_llik <- data.frame(Iteration = 1:max_iter, LogLikelihood = loglikelihoods)
p_llik <- ggplot(df_llik, aes(x = Iteration, y = LogLikelihood)) +
  geom_line(color = "blue", size = 1.2) +
  geom_point(color = "red", size = 2) +
  labs(title = "Log-Likelihood over EM Iterations",
       x = "Iteration", y = "Log-Likelihood") +
  theme_minimal()
print(p_llik)

# ----------------------------
# Plot: Evolution of Decay Rates over Iterations
# ----------------------------
df_lambda <- as.data.frame(lambda_hist)
df_lambda$Iteration <- 1:max_iter
df_lambda_long <- melt(df_lambda, id.vars = "Iteration", variable.name = "Dye", value.name = "Rate")

p_rates <- ggplot(df_lambda_long, aes(x = Iteration, y = Rate, color = Dye)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(title = "Evolution of Estimated Decay Rates over Iterations",
       x = "Iteration", y = "Decay Rate") +
  theme_minimal() +
  scale_color_viridis_d()
print(p_rates)






library(ggplot2)
library(dplyr)
library(tidyr)
library(reshape2)
library(viridis)
library(gganimate)

# ----------------------------
# Reconstruction Animation with Zero Background
# ----------------------------

# Define number of frames for the animation (each frame corresponds to a time t)
T_frames <- 10

# Create a full grid of pixel coordinates and a pixel identifier
grid_pixels <- expand.grid(i = 1:H, j = 1:W) %>%
  mutate(pixel = paste(i, j, sep = "_"))

# For each time frame t, compute the composite intensity (number of photons with tau > t)
intensity_list <- list()
for(t in 0:(T_frames - 1)){
  df_temp <- df_photons %>%
    filter(tau > t) %>%
    group_by(i, j) %>%
    summarise(intensity = n(), .groups = "drop")
  # Join with full grid so that pixels with zero counts are included
  df_temp <- full_join(grid_pixels, df_temp, by = c("i", "j"))
  df_temp$intensity[is.na(df_temp$intensity)] <- 0
  df_temp$frame <- t
  df_temp$pixel <- paste(df_temp$i, df_temp$j, sep = "_")
  intensity_list[[length(intensity_list) + 1]] <- df_temp
}
df_intensity <- bind_rows(intensity_list)

# Create a data frame from the estimated mixing proportions (pi_hat)
# (Assumes pi_hat has two columns for Dye 1 and Dye 2 and rows corresponding to unique_pixels)
df_pi <- data.frame(pixel = unique_pixels,
                    pi1 = pi_hat[,1],
                    pi2 = pi_hat[,2])
# For background pixels (not present in unique_pixels), assign zeros.
# Join and replace NAs with 0.
df_intensity <- left_join(df_intensity, df_pi, by = "pixel") %>%
  mutate(pi1 = ifelse(is.na(pi1), 0, pi1),
         pi2 = ifelse(is.na(pi2), 0, pi2))

# Compute the reconstructed intensities for Dye 1 and Dye 2:
#   Estimated Dye d intensity = Composite intensity * estimated π_d at that pixel.
df_intensity <- df_intensity %>% mutate(
  dye1 = intensity * pi1,
  dye2 = intensity * pi2
)

# Pivot the data to long format for faceting the animation into three panels.
df_long <- df_intensity %>% pivot_longer(
  cols = c("intensity", "dye1", "dye2"),
  names_to = "Component",
  values_to = "Value"
)

# Rename components for nicer facet labels
df_long$Component <- recode(df_long$Component,
                            intensity = "Composite",
                            dye1 = "Dye 1",
                            dye2 = "Dye 2")

# Determine an overall maximum value across all components (for fixed color scale)
max_value <- max(df_long$Value, na.rm = TRUE)

# Create the animated plot (3 panels) using gganimate
p_recon <- ggplot(df_long, aes(x = j, y = i, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, max_value)) +
  labs(title = 'Frame: {current_frame}', x = "Column", y = "Row") +
  theme_minimal() +
  scale_y_reverse() +
  facet_wrap(~ Component, ncol = 3) +
  transition_manual(frame) +
  ease_aes('linear')

# To view the animation interactively, use:
animate(p_recon, nframes = T_frames, fps = 2, width = 800, height = 400)

# ----------------------------
# Static Plot for the First Frame (Frame 0)
# ----------------------------

df_first <- df_long %>% filter(frame == 0)

p_first <- ggplot(df_first, aes(x = j, y = i, fill = Value)) +
  geom_tile() +
  scale_fill_viridis_c(limits = c(0, max_value)) +
  labs(title = 'Reconstruction at Frame 0', x = "Column", y = "Row") +
  theme_minimal() +
  scale_y_reverse() +
  facet_wrap(~ Component, ncol = 3)

# Print the static plot for the first frame
print(p_first)


# ----------------------------
# Fit k=3 Model to the Same Data and Compare with k=2
# ----------------------------

# (We assume that the following objects are already defined from the previous simulation/EM code:
#  df_photons, unique_pixels, pixel_indices, tau_vals, N (number of photons), H, W,
#  true_lambda_map, loglikelihoods, lambda_hat, pi_hat, and estimated_lambda_map from the k=2 model.)

# --- Fit k=3 Model ---
m3 <- 3
N_pix <- length(unique_pixels)
pi_hat_3 <- matrix(1/m3, nrow = N_pix, ncol = m3)   # initialize mixing proportions uniformly
lambda_hat_3 <- runif(m3, min = 1, max = 3)          # initialize global decay rates

max_iter <- 30
loglikelihoods_3 <- numeric(max_iter)
lambda_hist_3 <- matrix(0, nrow = max_iter, ncol = m3)

for (iter in 1:max_iter) {
  numerators <- matrix(0, nrow = N, ncol = m3)
  for (d in 1:m3) {
    numerators[, d] <- pi_hat_3[pixel_indices, d] * lambda_hat_3[d] * exp(-lambda_hat_3[d] * tau_vals)
  }
  resp3 <- numerators / rowSums(numerators)

  # Compute log-likelihood for the iteration
  loglikelihood_3 <- sum(log(rowSums(numerators)))
  loglikelihoods_3[iter] <- loglikelihood_3

  # Update global decay rates
  for (d in 1:m3) {
    lambda_hat_3[d] <- sum(resp3[, d]) / sum(resp3[, d] * tau_vals)
  }

  # Update pixel-specific mixing proportions using 3x3 neighborhood
  for (pix in 1:N_pix) {
    coords <- as.numeric(unlist(strsplit(unique_pixels[pix], "_")))
    i0 <- coords[1]
    j0 <- coords[2]

    # Determine neighbors (3x3 window, respecting grid boundaries)
    i_range <- max(1, i0 - 1):min(H, i0 + 1)
    j_range <- max(1, j0 - 1):min(W, j0 + 1)
    neighbors <- expand.grid(i = i_range, j = j_range)
    neighbor_ids <- paste(neighbors$i, neighbors$j, sep = "_")

    idx <- which(df_photons$pixel %in% neighbor_ids)
    if (length(idx) > 0) {
      pi_hat_3[pix, ] <- colMeans(resp3[idx, , drop = FALSE])
    }
  }
  # Normalize for numerical stability
  pi_hat_3 <- pmax(pi_hat_3, 1e-6)
  pi_hat_3 <- pi_hat_3 / rowSums(pi_hat_3)

  lambda_hist_3[iter, ] <- lambda_hat_3

  cat("k=3, Iteration", iter,
      ": Log-Likelihood =", round(loglikelihood_3, 2),
      "| Estimated Rates =", paste(round(lambda_hat_3, 3), collapse = ", "), "\n")
}

# --- Compute BIC for k=3 Model ---
# BIC = -2*log-likelihood + p*log(n)
# Here, n = number of photons, and p = free parameters = N_pix*(m3-1) + m3.
n_obs <- N  # number of photons
p3 <- N_pix * (m3 - 1) + m3
BIC_3 <- -2 * loglikelihood_3 + p3 * log(n_obs)
cat("\nBIC for k=3 model:", BIC_3, "\n")

# --- For k=2 Model, Retrieve BIC (from earlier run) ---
# (Assume k=2 model was fitted and final loglikelihood is stored in loglikelihoods[max_iter],
#  with m2 = 2, and free parameters p2 = N_pix*(2-1) + 2.)
m2 <- 2
p2 <- N_pix * (m2 - 1) + m2
loglikelihood_k2 <- loglikelihoods[max_iter]
BIC_2 <- -2 * loglikelihood_k2 + p2 * log(n_obs)
cat("BIC for k=2 model:", BIC_2, "\n")

# --- Compute Reconstruction (Composite λ Map) for k=3 ---
estimated_lambda_map_3 <- matrix(0, nrow = H, ncol = W)
for (k in 1:N_pix) {
  coords <- as.numeric(unlist(strsplit(unique_pixels[k], "_")))
  i_coord <- coords[1]
  j_coord <- coords[2]
  # Composite estimate: weighted average of the k=3 rates at that pixel.
  estimated_lambda_map_3[i_coord, j_coord] <- sum(pi_hat_3[k, ] * lambda_hat_3)
}

# --- Compute Reconstruction Errors ---
# Here, we compare against the true composite λ map (from simulation, based on k=2 true parameters).
reconstruction_error_3 <- mean(abs(true_lambda_map - estimated_lambda_map_3))
cat("Reconstruction MAE for k=3 model:", reconstruction_error_3, "\n")

# For k=2 model (if not already computed), assume estimated_lambda_map is from k=2 run.
reconstruction_error_2 <- mean(abs(true_lambda_map - estimated_lambda_map))
cat("Reconstruction MAE for k=2 model:", reconstruction_error_2, "\n")

# --- Plot Composite λ Maps for Visual Comparison ---
par(mfrow = c(1,2))
image(t(estimated_lambda_map[nrow(estimated_lambda_map):1, ]), col = viridis(100),
      main = "Estimated Composite λ Map (k=2)", axes = FALSE)
image(t(estimated_lambda_map_3[nrow(estimated_lambda_map_3):1, ]), col = viridis(100),
      main = "Estimated Composite λ Map (k=3)", axes = FALSE)
par(mfrow = c(1,1))

