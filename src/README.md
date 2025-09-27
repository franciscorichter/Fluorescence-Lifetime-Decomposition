# Fluorescence Lifetime Decomposition - Source Code
# This folder contains all the source code for the project

## Contents

### R Scripts
- **FullSudy.R**: Main analysis script for multi-k model comparison (k=1-5)
- **simulation2k.R**: Simulation study with k=2 and k=3 model comparison and visualization
- **twoDyeDataAnalysis.R**: Complete pipeline for real photon decay data analysis

### C++ Implementation
- **em_algorithm.cpp**: High-performance C++ implementation of EM algorithm with spatial regularization

## Usage

All scripts can be run from the project root directory:

```bash
# Run main analysis
Rscript src/FullSudy.R

# Run simulation study
Rscript src/simulation2k.R

# Run data analysis
Rscript src/twoDyeDataAnalysis.R
```

## Dependencies

- R with packages: Rcpp, ggplot2, dplyr, tidyr, viridis, gganimate, reticulate
- C++ compiler for Rcpp compilation
- Python (via reticulate) for NumPy data loading
