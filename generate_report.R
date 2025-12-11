#!/usr/bin/env Rscript
#' =============================================================================
#' CRISPR gRNA Counting - Automated QC Report Generator
#' =============================================================================
#'
#' This script generates a comprehensive HTML QC report with a single command.
#'
#' Usage:
#'   Rscript generate_report.R <counts_file.csv> [output_name]
#'
#' Examples:
#'   Rscript generate_report.R output/grna_counts.csv
#'   Rscript generate_report.R output/grna_counts.csv my_experiment_qc
#'   Rscript generate_report.R counts.csv report --min-count 20 --max-zero 15
#'
#' Options:
#'   --min-count N      Minimum count threshold (default: 10)
#'   --max-zero N       Maximum zero-count percentage (default: 20)
#'   --max-gini N       Maximum Gini coefficient (default: 0.85)
#'   --gc-threshold N   GC bias correlation threshold (default: 0.3)
#'   --output-dir DIR   Output directory (default: same as input)
#'
#' =============================================================================

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Show help
if (length(args) == 0 || args[1] %in% c("-h", "--help")) {
  cat("
CRISPR gRNA Counting - QC Report Generator
==========================================

Usage:
  Rscript generate_report.R <counts_file.csv> [output_name] [options]

Arguments:
  counts_file.csv    Path to the gRNA counts CSV file (required)
  output_name        Name for output report (optional, default: qc_report)

Options:
  --min-count N      Minimum count threshold (default: 10)
  --max-zero N       Maximum acceptable zero-count % (default: 20)
  --max-gini N       Maximum Gini coefficient (default: 0.85)
  --gc-threshold N   GC bias threshold (default: 0.3)
  --output-dir DIR   Output directory (default: same as input file)

Examples:
  # Basic usage
  Rscript generate_report.R output/grna_counts.csv

  # Custom output name
  Rscript generate_report.R output/grna_counts.csv experiment1_qc

  # With custom parameters
  Rscript generate_report.R counts.csv report --min-count 20 --max-zero 15

Output:
  - HTML report with interactive plots
  - QC summary and recommendations
  - Flagged guides table (exportable)
  - Gene-level statistics

")
  quit(status = 0)
}

# Required packages
required_packages <- c("rmarkdown", "knitr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing required package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
  }
}

# Parse arguments
counts_file <- args[1]
output_name <- "qc_report"

# Default parameters
params <- list(
  min_count_threshold = 10,
  max_zero_pct = 20,
  max_gini = 0.85,
  gc_bias_threshold = 0.3
)

output_dir <- NULL

# Parse remaining arguments
i <- 2
while (i <= length(args)) {
  arg <- args[i]

  if (arg == "--min-count" && i < length(args)) {
    params$min_count_threshold <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (arg == "--max-zero" && i < length(args)) {
    params$max_zero_pct <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (arg == "--max-gini" && i < length(args)) {
    params$max_gini <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (arg == "--gc-threshold" && i < length(args)) {
    params$gc_bias_threshold <- as.numeric(args[i + 1])
    i <- i + 2
  } else if (arg == "--output-dir" && i < length(args)) {
    output_dir <- args[i + 1]
    i <- i + 2
  } else if (!startsWith(arg, "--")) {
    output_name <- arg
    i <- i + 1
  } else {
    i <- i + 1
  }
}

# Validate input file
if (!file.exists(counts_file)) {
  stop(paste("Error: File not found:", counts_file))
}

# Set output directory
if (is.null(output_dir)) {
  output_dir <- dirname(counts_file)
  if (output_dir == ".") output_dir <- getwd()
}

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Find the Rmd template
script_dir <- dirname(sys.frame(1)$ofile)
if (is.null(script_dir) || script_dir == "") {
  script_dir <- getwd()
}

rmd_file <- file.path(script_dir, "qc_report.Rmd")

if (!file.exists(rmd_file)) {
  # Try current directory
  rmd_file <- "qc_report.Rmd"
}

if (!file.exists(rmd_file)) {
  stop("Error: Cannot find qc_report.Rmd template file")
}

# Generate output path
output_file <- file.path(output_dir, paste0(output_name, ".html"))

# Print info
cat("\n")
cat("=================================================\n")
cat("  CRISPR gRNA Counting - QC Report Generator\n")
cat("=================================================\n\n")
cat("Input file:     ", counts_file, "\n")
cat("Output file:    ", output_file, "\n")
cat("Template:       ", rmd_file, "\n")
cat("\nParameters:\n")
cat("  Min count threshold: ", params$min_count_threshold, "\n")
cat("  Max zero-count %:    ", params$max_zero_pct, "%\n")
cat("  Max Gini coefficient:", params$max_gini, "\n")
cat("  GC bias threshold:   ", params$gc_bias_threshold, "\n")
cat("\n")

# Add counts_file to params
params$counts_file <- normalizePath(counts_file)

# Render the report
cat("Generating report...\n\n")

tryCatch({
  rmarkdown::render(
    input = rmd_file,
    output_file = basename(output_file),
    output_dir = output_dir,
    params = params,
    quiet = FALSE,
    envir = new.env()
  )

  cat("\n")
  cat("=================================================\n")
  cat("  Report generated successfully!\n")
  cat("=================================================\n")
  cat("\nOutput saved to:", output_file, "\n")
  cat("\nOpen the HTML file in a web browser to view.\n\n")

}, error = function(e) {
  cat("\nError generating report:\n")
  cat(conditionMessage(e), "\n")
  quit(status = 1)
})
