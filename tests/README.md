# Test Suite
# This folder contains tests to validate the fluorescence lifetime decomposition project

## Contents

### Basic Functionality Test
- **test_basic_functionality.R**: Comprehensive test suite covering:
  - Project structure validation
  - C++ code compilation
  - Data loading and validation
  - R dependency checks
  - Basic functionality verification

## Running Tests

### Automated Testing
```bash
# Run all tests using the Makefile
make test

# Or run directly with R
Rscript tests/test_basic_functionality.R
```

### Individual Test Components

#### 1. Project Structure Test
- Verifies all required directories exist
- Checks for essential files
- Validates file organization

#### 2. C++ Compilation Test
- Tests if C++ code compiles successfully
- Verifies Rcpp integration
- Checks for compiler availability

#### 3. Data Loading Test
- Validates sample data can be loaded
- Checks data format and structure
- Verifies data quality metrics

#### 4. Dependency Test
- Confirms all required R packages are installed
- Checks package versions compatibility

#### 5. Functionality Test
- Tests basic data preparation pipeline
- Validates analysis setup
- Checks for common configuration issues

## Test Results

### Success Indicators
- ✅ All tests pass: Project is ready for use
- ⚠️ Some tests fail: Check specific error messages
- ❌ Critical failures: Project needs fixes before use

### Common Issues and Solutions

#### C++ Compilation Fails
```bash
# Install Xcode command line tools (macOS)
xcode-select --install

# Install build essentials (Ubuntu)
sudo apt-get install build-essential
```

#### Missing R Packages
```r
# Install required packages
install.packages(c("Rcpp", "ggplot2", "dplyr", "tidyr", "viridis"))
```

#### Data Loading Issues
- Check that `data/df_lifetimes.RData` exists
- Verify data file format and structure
- Ensure sufficient disk space

## Adding New Tests

### Test File Structure
```r
# Template for new test files
test_new_feature <- function() {
  cat("Test: New Feature\n")
  cat("================\n")

  tryCatch({
    # Test implementation
    cat("✅ New feature test passed\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ New feature test failed:", e$message, "\n")
    return(FALSE)
  })
}
```

### Test Best Practices
1. **Clear naming**: Use descriptive test function names
2. **Error handling**: Use tryCatch for robust error reporting
3. **Informative output**: Provide helpful error messages
4. **Independence**: Tests should not depend on each other
5. **Documentation**: Comment complex test logic

## Continuous Integration

### Automated Testing Setup
```bash
# Add to .github/workflows/ if using GitHub Actions
# Run tests on push/PR
# Notify on failures
```

### Pre-commit Hooks
```bash
# Run tests before commits
# Fail commit if tests don't pass
```

## Troubleshooting

### Test Environment Issues
- Ensure R and Rscript are in PATH
- Check R package installation permissions
- Verify C++ compiler installation

### Data-Specific Issues
- Check file permissions on data files
- Verify data file integrity
- Ensure data format compatibility

### Performance Issues
- Large datasets may slow tests
- Consider subsampling for quick validation
- Monitor memory usage for large tests

## Test Coverage

Current test coverage includes:
- ✅ Project structure and organization
- ✅ C++ compilation and linking
- ✅ Data loading and validation
- ✅ R dependency management
- ✅ Basic functionality setup

Future test additions:
- 🔄 Algorithm accuracy validation
- 🔄 Performance benchmarking
- 🔄 Integration tests with real data
- 🔄 Cross-platform compatibility
