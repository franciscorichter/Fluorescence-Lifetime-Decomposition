# Fluorescence Lifetime Decomposition

A comprehensive toolkit for analyzing fluorescence lifetime imaging microscopy (FLIM) data using Expectation-Maximization (EM) algorithms to decompose multi-exponential decay curves from multiple fluorophores or dyes.

## Overview

This project implements advanced statistical methods for decomposing fluorescence lifetime data into individual components corresponding to different molecular species or dyes. The core algorithm uses an EM approach with spatial regularization to estimate:

- **Global decay rates** (λ) for each molecular species
- **Pixel-specific mixing proportions** (π) indicating the contribution of each species at each spatial location
- **Model comparison** using Bayesian Information Criterion (BIC)

## Key Features

- **Multi-exponential decomposition**: Supports 1-5 component models
- **Spatial regularization**: Uses 3x3 neighborhood averaging for mixing proportions
- **C++ optimization**: Fast EM implementation via Rcpp
- **Real-time visualization**: Frame-by-frame reconstruction animations
- **Comprehensive diagnostics**: Log-likelihood tracking, BIC analysis, and parameter evolution plots

## Mathematical Foundation

### Fluorescence Decay Model

The fluorescence intensity at pixel (i,j) and time t follows:

```
I(t; i,j) = n_ij * Σ_d π_d(i,j) * λ_d * exp(-λ_d * t)
```

Where:
- `n_ij`: Total photon count at pixel (i,j)
- `π_d(i,j)`: Mixing proportion of dye d at pixel (i,j)
- `λ_d`: Global decay rate for dye d
- `t`: Time frame

### EM Algorithm

The algorithm iteratively estimates parameters by:

1. **E-step**: Compute responsibilities for each photon belonging to each dye
2. **M-step**: Update decay rates and mixing proportions using neighborhood regularization

## Installation

### Prerequisites

- R (version 4.0+)
- R packages: `Rcpp`, `ggplot2`, `dplyr`, `tidyr`, `viridis`, `gganimate`, `reticulate`
- C++ compiler (for Rcpp)

### Installation Steps

```bash
# Install required R packages
install.packages(c("Rcpp", "ggplot2", "dplyr", "tidyr", "viridis", "gganimate", "reticulate"))

# Clone or download the repository
git clone https://github.com/franciscorichter/Fluorescence-Lifetime-Decomposition.git
cd Fluorescence-Lifetime-Decomposition

# Compile C++ code (automatically done when sourcing in R)
```

### Data Requirements

The package expects fluorescence lifetime data in one of these formats:

1. **Preprocessed lifetime events**: RData file with `df_lifetimes` containing columns `tau`, `i`, `j`
2. **Raw time series**: NumPy .npy files with dimensions (rows × cols × time)
3. **Simulation data**: Generated internally for testing and validation

## Project Structure

This project follows a clean, organized structure for easy navigation and maintenance:

```
Fluorescence-Lifetime-Decomposition/
├── src/                 # Source code
│   ├── FullSudy.R              # Multi-k model analysis
│   ├── simulation2k.R          # Simulation study
│   ├── twoDyeDataAnalysis.R    # Data preprocessing pipeline
│   ├── em_algorithm.cpp        # C++ EM implementation
│   └── README.md              # Source code documentation
├── data/               # Data files
│   ├── df_lifetimes.RData     # Preprocessed lifetime events
│   └── README.md             # Data documentation
├── docs/               # Documentation
│   ├── README.md             # Main documentation
│   ├── api/                  # Technical documentation
│   │   ├── FullSudy_documentation.md
│   │   ├── simulation2k_documentation.md
│   │   ├── twoDyeDataAnalysis_documentation.md
│   │   ├── em_algorithm_documentation.md
│   │   └── data_formats_documentation.md
├── examples/           # Usage examples and tutorials
│   ├── usage_examples.md     # Comprehensive tutorials
│   └── README.md            # Examples documentation
├── tests/              # Test suite
│   ├── test_basic_functionality.R  # Basic functionality tests
│   └── README.md           # Test documentation
├── output/             # Generated results (auto-created)
├── scripts/            # Utility scripts (future use)
├── .gitignore         # Git ignore rules
├── Makefile           # Build and test automation
└── README.md          # This file
```

### Getting Started

1. **Quick Analysis**: Run `make test` to verify everything works
2. **Basic Usage**: See `examples/usage_examples.md` for tutorials
3. **Development**: Check `src/README.md` for source code details
4. **Documentation**: Browse `docs/api/` for technical details

### Quick Start Example

```r
# Load required libraries
library(Rcpp)
library(ggplot2)
library(dplyr)
library(viridis)

# Load preprocessed lifetime data
load("data/df_lifetimes.RData")

# Source the C++ EM implementation
Rcpp::sourceCpp("src/em_algorithm.cpp")

# Prepare data for analysis
df_events <- df_lifetimes
df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
unique_pixels <- unique(df_events$pixel)
df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

# Run EM for k=2 model
source("src/twoDyeDataAnalysis.R")  # Contains the main analysis functions

# Visualize results
# (Results will be displayed as plots and saved to workspace)
```

### Advanced Usage

For custom analysis with different parameters:

```r
# Custom EM wrapper function
em_estimation <- function(df_events, k, n_rows, n_cols, max_iter = 30) {
  # Prepare data
  df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
  unique_pixels <- unique(df_events$pixel)
  df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

  N_pix <- length(unique_pixels)
  coords <- do.call(rbind, strsplit(unique_pixels, "_"))
  coords <- matrix(as.integer(coords), ncol = 2)
  colnames(coords) <- c("i", "j")

  # Initialize parameters
  pi_init <- matrix(1/k, nrow = N_pix, ncol = k)
  lambda_init <- runif(k, min = 1, max = 3)

  # Run EM
  res <- EMAlgorithmCpp(df_events$pixel_idx, df_events$tau,
                       pi_init, lambda_init, coords, n_rows, n_cols, max_iter)
  return(res)
}

# Run for multiple k values
results_list <- list()
for(k in 1:5){
  res <- em_estimation(df_events, k, n_rows, n_cols)
  results_list[[k]] <- res
}

# Compare using BIC
bic_values <- numeric(5)
for(k in 1:5){
  res <- results_list[[k]]
  p_params <- length(unique_pixels) * (k - 1) + k
  bic_values[k] <- -2 * res$loglikelihoods[length(res$loglikelihoods)] + p_params * log(nrow(df_events))
}
```

## File Descriptions

### Core Scripts

- **`src/FullSudy.R`**: Main analysis script for multi-k model comparison (k=1-5)
- **`src/twoDyeDataAnalysis.R`**: Complete pipeline for real photon decay data analysis
- **`src/simulation2k.R`**: Simulation study with k=2 and k=3 model comparison
- **`src/em_algorithm.cpp`**: C++ implementation of the EM algorithm with spatial regularization

### Data Files

- **`data/df_lifetimes.RData`**: Preprocessed lifetime events (τ, i, j coordinates)
- **`best_fit.npy`**: Example raw fluorescence time series data

### Documentation

- **`docs/README.md`**: This main documentation file
- **`docs/api/`**: Detailed API and technical documentation
- **`examples/usage_examples.md`**: Comprehensive tutorials and examples
- **`tests/test_basic_functionality.R`**: Test suite for validation

## Output and Visualization

### Model Comparison

The toolkit provides several diagnostic plots:

- **BIC vs. Number of Components**: Model selection criterion
- **Log-Likelihood Evolution**: Convergence monitoring
- **Parameter Evolution**: Decay rate estimates over EM iterations
- **Composite λ Maps**: Spatially resolved decay rate estimates

### Reconstruction Visualizations

- **Frame-by-frame animations**: Time evolution of fluorescence decay
- **Component decomposition**: Individual dye contributions
- **Composite reconstructions**: Model-predicted intensity maps

### Example Output Interpretation

```
Final Comparison of Decay Rates:
  Dye True_Rate Estimated_Rate Absolute_Error
1 Dye1       0.5           0.52           0.02
2 Dye2       1.5           1.48           0.02

BIC for k=2 model: 15432.1
BIC for k=3 model: 15678.3  # Higher BIC indicates overparameterization
```

## Algorithm Details

### EM Implementation

The EM algorithm uses:

- **E-step**: Soft assignment of photons to dyes based on current parameters
- **M-step**: 
  - Update global decay rates using weighted averages
  - Update pixel-specific mixing proportions using 3×3 spatial neighborhoods
  - Normalize proportions for numerical stability

### Spatial Regularization

Mixing proportions are smoothed using local neighborhood averaging:

```cpp
// 3x3 neighborhood update in C++
for each pixel (i0,j0):
  for each neighbor (i,j) in 3x3 window:
    accumulate responsibilities from photons in neighborhood
  average to get new mixing proportions
```

### Convergence Criteria

- Maximum iterations: 30 (default)
- Log-likelihood change monitoring
- Parameter stabilization

## Performance Considerations

- **Computational complexity**: O(n_pixels × n_components × n_iterations)
- **Memory usage**: Scales with number of photons and pixels
- **C++ acceleration**: ~10-50x speedup over pure R implementation

## Troubleshooting

### Common Issues

1. **No lifetime events detected**: Adjust `min_photons` threshold in preprocessing
2. **EM convergence issues**: Increase `max_iter` or adjust initialization
3. **Memory errors**: Reduce spatial resolution or use data subsampling

### Parameter Tuning

- **Initialization**: Decay rates typically between 0.1-5.0 for fluorescence
- **Neighborhood size**: 3×3 provides good balance of smoothness vs. detail
- **Convergence**: Monitor log-likelihood plateau for convergence

## Citation

If you use this software in your research, please cite:

```
Richter, F. (2024). Fluorescence Lifetime Decomposition: EM-based Analysis of Multi-exponential Decay Data.
GitHub repository: https://github.com/franciscorichter/Fluorescence-Lifetime-Decomposition
```

## Contributing

Contributions are welcome! Areas for improvement:

- GPU acceleration
- Alternative regularization methods
- Additional model selection criteria
- Integration with common FLIM analysis packages

## License

This project is released under the MIT License. See LICENSE file for details.

## Support

For questions or issues, please:

1. Check the troubleshooting section
2. Review example scripts for usage patterns
3. Open an issue on the GitHub repository
4. Contact the maintainer: franciscorichter [at] github

---

**Note**: This toolkit is designed for research use. Always validate results with known standards and consider the physical constraints of your fluorescence system.
