#!/usr/bin/env Rscript
# Basic functionality tests for Fluorescence Lifetime Decomposition

# Load required libraries
suppressPackageStartupMessages({
  library(Rcpp)
  library(ggplot2)
  library(dplyr)
})

# Test 1: Check project structure
test_project_structure <- function() {
  cat("Test 1: Project Structure\n")
  cat("=========================\n")

  # Check required directories
  required_dirs <- c("src", "data", "docs", "examples")
  missing_dirs <- c()

  for(dir in required_dirs) {
    if(!dir.exists(dir)) {
      missing_dirs <- c(missing_dirs, dir)
    }
  }

  if(length(missing_dirs) > 0) {
    cat("❌ Missing directories:", paste(missing_dirs, collapse = ", "), "\n")
    return(FALSE)
  } else {
    cat("✅ All required directories exist\n")
  }

  # Check required files
  required_files <- c("src/em_algorithm.cpp", "data/df_lifetimes.RData")
  missing_files <- c()

  for(file in required_files) {
    if(!file.exists(file)) {
      missing_files <- c(missing_files, file)
    }
  }

  if(length(missing_files) > 0) {
    cat("❌ Missing files:", paste(missing_files, collapse = ", "), "\n")
    return(FALSE)
  } else {
    cat("✅ All required files exist\n")
  }

  cat("\n")
  return(TRUE)
}

# Test 2: Check C++ compilation
test_cpp_compilation <- function() {
  cat("Test 2: C++ Compilation\n")
  cat("=======================\n")

  tryCatch({
    # Try to compile C++ code
    Rcpp::sourceCpp("src/em_algorithm.cpp")
    cat("✅ C++ compilation successful\n")
    cat("\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ C++ compilation failed:", e$message, "\n")
    cat("💡 Make sure you have a C++ compiler installed\n")
    cat("   On macOS: install Xcode command line tools\n")
    cat("   On Ubuntu: install build-essential\n")
    cat("\n")
    return(FALSE)
  })
}

# Test 3: Check data loading
test_data_loading <- function() {
  cat("Test 3: Data Loading\n")
  cat("===================\n")

  tryCatch({
    # Try to load sample data
    load("data/df_lifetimes.RData")

    # Basic data validation
    if(!exists("df_lifetimes")) {
      cat("❌ Data loaded but df_lifetimes not found\n")
      return(FALSE)
    }

    cat("✅ Data loaded successfully\n")
    cat("   Dimensions:", nrow(df_lifetimes), "x", ncol(df_lifetimes), "\n")
    cat("   Columns:", paste(colnames(df_lifetimes), collapse = ", "), "\n")

    # Check data quality
    if(nrow(df_lifetimes) < 100) {
      cat("⚠️  Warning: Small dataset (", nrow(df_lifetimes), "events)\n")
    } else {
      cat("✅ Sufficient data for analysis\n")
    }

    cat("\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Data loading failed:", e$message, "\n")
    cat("💡 Check that data/df_lifetimes.RData exists\n")
    cat("\n")
    return(FALSE)
  })
}

# Test 4: Check R dependencies
test_dependencies <- function() {
  cat("Test 4: R Dependencies\n")
  cat("=====================\n")

  required_packages <- c("Rcpp", "ggplot2", "dplyr", "tidyr", "viridis")
  missing_packages <- c()

  for(pkg in required_packages) {
    if(!requireNamespace(pkg, quietly = TRUE)) {
      missing_packages <- c(missing_packages, pkg)
    }
  }

  if(length(missing_packages) > 0) {
    cat("❌ Missing packages:", paste(missing_packages, collapse = ", "), "\n")
    cat("💡 Install with: install.packages(c('", paste(missing_packages, collapse = "', '"), "'))\n")
    cat("\n")
    return(FALSE)
  } else {
    cat("✅ All required packages are installed\n")
    cat("\n")
    return(TRUE)
  }
}

# Test 5: Basic functionality
test_basic_functionality <- function() {
  cat("Test 5: Basic Functionality\n")
  cat("===========================\n")

  tryCatch({
    # Load data
    load("data/df_lifetimes.RData")

    # Prepare data for analysis
    df_events <- df_lifetimes
    df_events$pixel <- with(df_events, paste(i, j, sep = "_"))
    unique_pixels <- unique(df_events$pixel)
    df_events$pixel_idx <- as.integer(factor(df_events$pixel, levels = unique_pixels))

    n_rows <- max(df_events$i)
    n_cols <- max(df_events$j)

    # Test EM function exists (without running full analysis)
    if(!file.exists("src/em_algorithm.cpp")) {
      cat("❌ EM algorithm source not found\n")
      return(FALSE)
    }

    # Test data preparation
    if(length(unique_pixels) < 10) {
      cat("⚠️  Warning: Very few pixels with data (", length(unique_pixels), ")\n")
    } else {
      cat("✅ Data preparation looks good\n")
      cat("   Pixels with data:", length(unique_pixels), "\n")
      cat("   Grid size:", n_rows, "x", n_cols, "\n")
    }

    cat("\n")
    return(TRUE)
  }, error = function(e) {
    cat("❌ Basic functionality test failed:", e$message, "\n")
    cat("\n")
    return(FALSE)
  })
}

# Run all tests
run_all_tests <- function() {
  cat("Fluorescence Lifetime Decomposition - Test Suite\n")
  cat("===============================================\n\n")

  tests <- list(
    test_project_structure,
    test_cpp_compilation,
    test_data_loading,
    test_dependencies,
    test_basic_functionality
  )

  results <- logical(length(tests))

  for(i in seq_along(tests)) {
    results[i] <- tests[[i]]()
  }

  cat("Test Summary:\n")
  cat("=============\n")
  passed <- sum(results)
  total <- length(results)

  cat("Passed:", passed, "/", total, "\n")

  if(passed == total) {
    cat("🎉 All tests passed! Project is ready to use.\n")
    return(0)
  } else {
    cat("⚠️  Some tests failed. Check output above for details.\n")
    return(1)
  }
}

# Main execution
if(!interactive()) {
  exit_code <- run_all_tests()
  quit(save = "no", status = exit_code)
} else {
  run_all_tests()
}
