# Makefile for Fluorescence Lifetime Decomposition Project

.PHONY: help clean build test docs examples data

# Default target
help:
	@echo "Fluorescence Lifetime Decomposition Project"
	@echo "========================================"
	@echo ""
	@echo "Available targets:"
	@echo "  help       - Show this help message"
	@echo "  build      - Compile C++ code"
	@echo "  test       - Run basic functionality tests"
	@echo "  examples   - Run example analyses"
	@echo "  clean      - Remove generated files"
	@echo "  install    - Install R dependencies"
	@echo "  docs       - Show documentation structure"
	@echo "  structure  - Show project structure"
	@echo ""

# Project structure
structure:
	@echo "Project Structure:"
	@echo "=================="
	@find . -type d -name ".git" -prune -o -type d -print | sort
	@echo ""
	@echo "Key Directories:"
	@echo "  src/       - Source code (R scripts, C++)"
	@echo "  data/      - Data files and datasets"
	@echo "  docs/      - Documentation"
	@echo "  examples/  - Usage examples and tutorials"
	@echo "  tests/     - Test files"
	@echo "  output/    - Generated results and plots"

# Build C++ code
build:
	@echo "Compiling C++ code..."
	@Rscript -e "Rcpp::sourceCpp('src/em_algorithm.cpp'); cat('C++ compilation successful!\n')"

# Install R dependencies
install:
	@echo "Installing R dependencies..."
	@Rscript -e "install.packages(c('Rcpp', 'ggplot2', 'dplyr', 'tidyr', 'viridis', 'gganimate', 'reticulate'), repos='https://cran.rstudio.com')"
	@echo "Dependencies installed!"

# Run basic tests
test:
	@echo "Running basic functionality tests..."
	@Rscript -e "source('src/twoDyeDataAnalysis.R'); cat('Basic analysis test completed!\n')"

# Run examples
examples:
	@echo "Running example analyses..."
	@Rscript src/simulation2k.R

# Clean generated files
clean:
	@echo "Cleaning generated files..."
	@rm -rf output/*.RData
	@rm -rf output/*.pdf
	@rm -rf output/*.png
	@find . -name "*.Rout" -delete
	@find . -name ".RData" -path "./data/*" -prune -o -name "*.RData" -delete
	@echo "Cleanup completed!"

# Documentation
docs:
	@echo "Documentation Structure:"
	@echo "======================="
	@echo "Main Documentation:"
	@echo "  README.md                  - Project overview and installation"
	@echo ""
	@echo "API Documentation (docs/api/):"
	@ls docs/api/*.md | sed 's|docs/api/|  - |'
	@echo ""
	@echo "Examples and Tutorials (examples/):"
	@ls examples/*.md | sed 's|examples/|  - |'
	@echo ""
	@echo "Getting Started:"
	@echo "1. Read main README.md for project overview"
	@echo "2. Check examples/usage_examples.md for tutorials"
	@echo "3. Review docs/api/ for technical details"
