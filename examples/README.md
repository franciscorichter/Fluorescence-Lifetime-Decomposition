# Examples and Tutorials
# This folder contains practical examples and tutorials for using the toolkit

## Contents

### Usage Examples
- **usage_examples.md**: Comprehensive tutorials and examples covering:
  - Basic model fitting and comparison
  - Visualization and reconstruction
  - Advanced usage patterns
  - Troubleshooting common issues
  - Performance optimization tips

## Tutorial Topics

### Quick Start
- Basic EM fitting with k=2 model
- Model comparison using BIC
- Simple visualization examples

### Advanced Usage
- Custom parameter initialization
- Bootstrap uncertainty estimation
- Batch processing multiple files
- Integration with experimental workflows

### Troubleshooting
- Convergence issues and solutions
- Memory optimization for large datasets
- Parameter validation techniques
- Performance debugging

## Running Examples

### Basic Example
```r
# From project root directory
cd /path/to/Fluorescence-Lifetime-Decomposition

# Load data and run basic analysis
Rscript -e "
load('data/df_lifetimes.RData')
source('src/twoDyeDataAnalysis.R')
"
```

### Tutorial Scripts
The examples can be adapted into runnable R scripts:

```r
# Save as tutorial_1_basic_analysis.R
source('src/twoDyeDataAnalysis.R')  # Your basic analysis here
```

## Learning Path

1. **Start Here**: usage_examples.md - Overview and quick start
2. **Core Concepts**: Understand EM algorithm and fluorescence lifetime theory
3. **Basic Usage**: Run simple examples with provided data
4. **Advanced Topics**: Customize for your specific use case
5. **Troubleshooting**: Use diagnostic tools and optimization tips

## Contributing Examples

When adding new examples:
1. Include clear problem statement
2. Provide step-by-step instructions
3. Add expected outputs/results
4. Include troubleshooting notes
5. Update this README

## Documentation

For detailed documentation, see:
- Main project docs: `docs/README.md`
- API documentation: `docs/api/`
- Technical details: Individual script documentations
