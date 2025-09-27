# Data Files
# This folder contains data files used by the fluorescence lifetime decomposition project

## Contents

### Sample Datasets
- **df_lifetimes.RData**: Preprocessed lifetime events data frame with columns:
  - `tau`: Lifetime values (discrete time units)
  - `i`: Row coordinates (1-based)
  - `j`: Column coordinates (1-based)

## Data Format

### Lifetime Events Format
The primary data format used for EM analysis:

```r
# Example structure
df_lifetimes <- data.frame(
  tau = c(1, 2, 1, 3, ...),    # Lifetime values
  i = c(1, 1, 1, 1, ...),       # Row coordinates
  j = c(1, 1, 2, 2, ...)        # Column coordinates
)
```

### Data Requirements
- **Photon coverage**: >100 photons per pixel for reliable estimation
- **Lifetime range**: Typically 1-20 time units
- **Spatial coverage**: Complete grid or well-defined regions

## Usage

Load data in R scripts:

```r
# Load lifetime data
load("data/df_lifetimes.RData")

# Use in analysis
source("src/twoDyeDataAnalysis.R")  # Processes the loaded data
```

## Adding New Data

1. **Preprocessing**: Use `src/twoDyeDataAnalysis.R` to preprocess raw time series
2. **Format**: Ensure data follows the lifetime events format
3. **Save**: Store as .RData files in this directory

## Data Quality

Before using new datasets:
1. Verify coordinate ranges match expected grid dimensions
2. Check photon count distribution
3. Validate lifetime value ranges
4. Ensure sufficient spatial coverage
