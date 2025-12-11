#!/usr/bin/env Rscript
#' =============================================================================
#' CRISPR gRNA Count QC and Visualization
#' =============================================================================
#'
#' This script provides quality control analysis and visualization for gRNA
#' counting results from CRISPR screening experiments.
#'
#' Features:
#' - Count distribution analysis (histograms, density plots, boxplots)
#' - Outlier detection (IQR method, MAD method)
#' - Low-count guide flagging
#' - GC content bias visualization
#' - Guide-level QC warnings for downstream DEG analysis
#' - Gene-level summary statistics
#'
#' Usage:
#'   Rscript qc_visualization.R <counts_file.csv> [output_dir]
#'
#' Or source in R:
#'   source("qc_visualization.R")
#'   results <- run_qc_analysis("output/grna_counts.csv")
#' =============================================================================

# Load required libraries (install if needed)
required_packages <- c("ggplot2", "dplyr", "tidyr", "scales", "gridExtra")

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste("Installing package:", pkg))
    install.packages(pkg, repos = "https://cloud.r-project.org/")
  }
  library(pkg, character.only = TRUE)
}

# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

#' Calculate GC content of a DNA sequence
#' @param seq DNA sequence string
#' @return GC content as proportion (0-1)
calculate_gc <- function(seq) {
  seq <- toupper(seq)
  gc_count <- nchar(gsub("[^GC]", "", seq))
  total <- nchar(gsub("[^ATGC]", "", seq))
  if (total == 0) return(NA)
  return(gc_count / total)
}

#' Detect outliers using IQR method
#' @param x numeric vector
#' @param k multiplier for IQR (default 1.5)
#' @return logical vector indicating outliers
detect_outliers_iqr <- function(x, k = 1.5) {
  q1 <- quantile(x, 0.25, na.rm = TRUE)
  q3 <- quantile(x, 0.75, na.rm = TRUE)
  iqr <- q3 - q1
  lower <- q1 - k * iqr
  upper <- q3 + k * iqr
  return(x < lower | x > upper)
}

#' Detect outliers using MAD (Median Absolute Deviation) method
#' @param x numeric vector
#' @param threshold MAD threshold (default 3)
#' @return logical vector indicating outliers
detect_outliers_mad <- function(x, threshold = 3) {
  med <- median(x, na.rm = TRUE)
  mad_val <- mad(x, na.rm = TRUE)
  if (mad_val == 0) return(rep(FALSE, length(x)))
  return(abs(x - med) / mad_val > threshold)
}

#' Calculate Gini coefficient
#' @param x numeric vector of counts
#' @return Gini coefficient (0-1)
calculate_gini <- function(x) {
  x <- sort(x[!is.na(x)])
  n <- length(x)
  if (n == 0 || sum(x) == 0) return(0)
  cumx <- cumsum(x)
  gini <- (2 * sum((1:n) * x) - (n + 1) * sum(x)) / (n * sum(x))
  return(max(0, min(1, gini)))
}

# =============================================================================
# QC FLAG FUNCTIONS
# =============================================================================

#' Generate QC flags for each guide
#' @param counts_df data frame with count data
#' @return data frame with QC flags added
add_qc_flags <- function(counts_df) {

  # Ensure we have count column
  if (!"count" %in% names(counts_df)) {
    stop("counts_df must have a 'count' column")
  }

  counts <- counts_df$count
  log_counts <- log10(counts + 1)

  # Calculate thresholds
  median_count <- median(counts, na.rm = TRUE)
  mean_count <- mean(counts, na.rm = TRUE)
  q10 <- quantile(counts, 0.10, na.rm = TRUE)
  q90 <- quantile(counts, 0.90, na.rm = TRUE)

  # Initialize flags
  counts_df <- counts_df %>%
    mutate(
      # Zero count flag
      flag_zero_count = count == 0,

      # Very low count (bottom 10% or < 10 reads)
      flag_low_count = count < max(q10, 10),

      # Very high count (potential outlier)
      flag_high_outlier = detect_outliers_iqr(counts, k = 3) & counts > median_count,

      # MAD-based outlier
      flag_mad_outlier = detect_outliers_mad(log_counts, threshold = 3),

      # Below median (might be depleted)
      flag_below_median = count < median_count,

      # Extreme values (top/bottom 1%)
      flag_extreme = count < quantile(counts, 0.01) | count > quantile(counts, 0.99)
    )

  # Add GC content if sequence column exists
  if ("sequence" %in% names(counts_df)) {
    counts_df$gc_content <- sapply(counts_df$sequence, calculate_gc)

    # Flag extreme GC content
    counts_df <- counts_df %>%
      mutate(
        flag_extreme_gc = gc_content < 0.25 | gc_content > 0.75
      )
  }

  # Create overall warning flag
  counts_df <- counts_df %>%
    mutate(
      qc_warning = flag_zero_count | flag_low_count | flag_high_outlier | flag_extreme,
      qc_status = case_when(
        flag_zero_count ~ "FAIL: Zero count",
        flag_low_count ~ "WARN: Low count",
        flag_high_outlier ~ "WARN: High outlier",
        flag_extreme ~ "WARN: Extreme value",
        TRUE ~ "PASS"
      )
    )

  return(counts_df)
}

# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

#' Create count distribution histogram
#' @param counts_df data frame with count data
#' @return ggplot object
plot_count_distribution <- function(counts_df) {

  p <- ggplot(counts_df, aes(x = count)) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    geom_vline(aes(xintercept = median(count)), color = "red", linetype = "dashed", size = 1) +
    geom_vline(aes(xintercept = mean(count)), color = "orange", linetype = "dashed", size = 1) +
    scale_x_continuous(labels = scales::comma) +
    labs(
      title = "gRNA Count Distribution",
      subtitle = paste("Red = Median (", round(median(counts_df$count)),
                      "), Orange = Mean (", round(mean(counts_df$count)), ")", sep = ""),
      x = "Read Count",
      y = "Number of Guides"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  return(p)
}

#' Create log-transformed count distribution
#' @param counts_df data frame with count data
#' @return ggplot object
plot_log_distribution <- function(counts_df) {

  counts_df$log_count <- log10(counts_df$count + 1)

  p <- ggplot(counts_df, aes(x = log_count)) +
    geom_histogram(bins = 50, fill = "darkgreen", color = "white", alpha = 0.7) +
    geom_density(aes(y = ..count.. * 0.1), color = "black", size = 1) +
    labs(
      title = "Log10 Count Distribution",
      subtitle = "Log10(count + 1) transformation",
      x = "Log10(Count + 1)",
      y = "Number of Guides"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  return(p)
}

#' Create boxplot by gene (top genes)
#' @param counts_df data frame with count data
#' @param top_n number of top genes to show
#' @return ggplot object
plot_gene_boxplot <- function(counts_df, top_n = 20) {

  if (!"gene_name" %in% names(counts_df)) {
    message("No gene_name column found, skipping gene boxplot")
    return(NULL)
  }

  # Get top genes by total count
  top_genes <- counts_df %>%
    group_by(gene_name) %>%
    summarize(total = sum(count), .groups = "drop") %>%
    arrange(desc(total)) %>%
    head(top_n) %>%
    pull(gene_name)

  plot_data <- counts_df %>%
    filter(gene_name %in% top_genes) %>%
    mutate(gene_name = factor(gene_name, levels = top_genes))

  p <- ggplot(plot_data, aes(x = gene_name, y = count, fill = gene_name)) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.5) +
    scale_y_log10(labels = scales::comma) +
    labs(
      title = paste("Count Distribution by Gene (Top", top_n, ")"),
      x = "Gene",
      y = "Count (log scale)"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14)
    )

  return(p)
}

#' Create QC status summary plot
#' @param counts_df data frame with QC flags
#' @return ggplot object
plot_qc_summary <- function(counts_df) {

  if (!"qc_status" %in% names(counts_df)) {
    counts_df <- add_qc_flags(counts_df)
  }

  qc_summary <- counts_df %>%
    count(qc_status) %>%
    mutate(
      pct = n / sum(n) * 100,
      label = paste0(n, " (", round(pct, 1), "%)")
    )

  # Color mapping
  colors <- c(
    "PASS" = "forestgreen",
    "WARN: Low count" = "orange",
    "WARN: High outlier" = "darkorange",
    "WARN: Extreme value" = "coral",
    "FAIL: Zero count" = "red"
  )

  p <- ggplot(qc_summary, aes(x = reorder(qc_status, -n), y = n, fill = qc_status)) +
    geom_bar(stat = "identity", alpha = 0.8) +
    geom_text(aes(label = label), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = colors) +
    labs(
      title = "Guide QC Status Summary",
      subtitle = "Flags for downstream DEG analysis",
      x = "QC Status",
      y = "Number of Guides"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1),
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    ) +
    ylim(0, max(qc_summary$n) * 1.15)

  return(p)
}

#' Create GC content vs count scatter plot
#' @param counts_df data frame with count and sequence data
#' @return ggplot object
plot_gc_bias <- function(counts_df) {

  if (!"gc_content" %in% names(counts_df)) {
    if ("sequence" %in% names(counts_df)) {
      counts_df$gc_content <- sapply(counts_df$sequence, calculate_gc)
    } else {
      message("No sequence column found, skipping GC bias plot")
      return(NULL)
    }
  }

  # Calculate correlation
  cor_val <- cor(counts_df$gc_content, log10(counts_df$count + 1),
                 use = "complete.obs", method = "spearman")

  p <- ggplot(counts_df, aes(x = gc_content, y = count + 1)) +
    geom_point(alpha = 0.4, color = "steelblue") +
    geom_smooth(method = "loess", color = "red", se = TRUE, alpha = 0.2) +
    scale_y_log10(labels = scales::comma) +
    scale_x_continuous(labels = scales::percent) +
    labs(
      title = "GC Content vs Count",
      subtitle = paste("Spearman correlation:", round(cor_val, 3)),
      x = "GC Content",
      y = "Count (log scale)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  return(p)
}

#' Create cumulative distribution plot
#' @param counts_df data frame with count data
#' @return ggplot object
plot_cumulative <- function(counts_df) {

  counts_df <- counts_df %>%
    arrange(desc(count)) %>%
    mutate(
      rank = row_number(),
      cumsum = cumsum(count),
      cum_pct = cumsum / sum(count) * 100
    )

  p <- ggplot(counts_df, aes(x = rank, y = cum_pct)) +
    geom_line(color = "steelblue", size = 1) +
    geom_hline(yintercept = c(50, 80, 95), linetype = "dashed", color = "gray50") +
    labs(
      title = "Cumulative Count Distribution",
      subtitle = "Shows count concentration (Lorenz curve)",
      x = "Guide Rank (sorted by count)",
      y = "Cumulative % of Total Counts"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  return(p)
}

#' Create outlier visualization
#' @param counts_df data frame with count data
#' @return ggplot object
plot_outliers <- function(counts_df) {

  if (!"flag_high_outlier" %in% names(counts_df)) {
    counts_df <- add_qc_flags(counts_df)
  }

  counts_df <- counts_df %>%
    mutate(
      outlier_type = case_when(
        flag_zero_count ~ "Zero",
        flag_low_count ~ "Low",
        flag_high_outlier ~ "High",
        TRUE ~ "Normal"
      ),
      log_count = log10(count + 1)
    )

  colors <- c("Normal" = "gray70", "Zero" = "red", "Low" = "orange", "High" = "purple")

  p <- ggplot(counts_df, aes(x = seq_along(count), y = log_count, color = outlier_type)) +
    geom_point(alpha = 0.6, size = 1.5) +
    scale_color_manual(values = colors) +
    labs(
      title = "Guide Counts with Outliers Highlighted",
      subtitle = "Red=Zero, Orange=Low, Purple=High outliers",
      x = "Guide Index",
      y = "Log10(Count + 1)",
      color = "Status"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 10, color = "gray40")
    )

  return(p)
}

# =============================================================================
# MAIN QC ANALYSIS FUNCTION
# =============================================================================

#' Run complete QC analysis on count data
#' @param counts_file path to counts CSV file
#' @param output_dir directory for output files (default: same as input)
#' @return list with QC results and flagged data
run_qc_analysis <- function(counts_file, output_dir = NULL) {

  message("=== CRISPR gRNA Count QC Analysis ===\n")

  # Read data
  message("Reading count data...")
  counts_df <- read.csv(counts_file, stringsAsFactors = FALSE)
  message(paste("  Loaded", nrow(counts_df), "guides"))

  # Set output directory
  if (is.null(output_dir)) {
    output_dir <- dirname(counts_file)
  }
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Add QC flags
  message("\nRunning QC checks...")
  counts_df <- add_qc_flags(counts_df)

  # Calculate summary statistics
  message("\nCalculating statistics...")
  stats <- list(
    total_guides = nrow(counts_df),
    total_reads = sum(counts_df$count),
    mean_count = mean(counts_df$count),
    median_count = median(counts_df$count),
    sd_count = sd(counts_df$count),
    cv = sd(counts_df$count) / mean(counts_df$count),
    gini = calculate_gini(counts_df$count),
    zero_count_guides = sum(counts_df$flag_zero_count),
    low_count_guides = sum(counts_df$flag_low_count),
    high_outliers = sum(counts_df$flag_high_outlier),
    qc_pass = sum(counts_df$qc_status == "PASS"),
    qc_warn = sum(grepl("WARN", counts_df$qc_status)),
    qc_fail = sum(grepl("FAIL", counts_df$qc_status))
  )

  # Print summary
  message("\n=== QC SUMMARY ===")
  message(paste("Total guides:", stats$total_guides))
  message(paste("Total reads:", format(stats$total_reads, big.mark = ",")))
  message(paste("Mean count:", round(stats$mean_count, 1)))
  message(paste("Median count:", round(stats$median_count, 1)))
  message(paste("CV (coefficient of variation):", round(stats$cv, 2)))
  message(paste("Gini coefficient:", round(stats$gini, 3)))
  message(paste("\nQC Results:"))
  message(paste("  PASS:", stats$qc_pass, "(", round(stats$qc_pass/stats$total_guides*100, 1), "%)"))
  message(paste("  WARN:", stats$qc_warn, "(", round(stats$qc_warn/stats$total_guides*100, 1), "%)"))
  message(paste("  FAIL:", stats$qc_fail, "(", round(stats$qc_fail/stats$total_guides*100, 1), "%)"))

  # Generate warnings for DEG analysis
  message("\n=== WARNINGS FOR DOWNSTREAM DEG ANALYSIS ===")

  if (stats$zero_count_guides > 0) {
    message(paste("! WARNING:", stats$zero_count_guides, "guides have zero counts - will cause issues in DEG"))
  }

  if (stats$gini > 0.8) {
    message("! WARNING: High Gini coefficient (>0.8) - very unequal count distribution")
    message("           Consider checking for technical issues or strong selection")
  }

  if (stats$cv > 2) {
    message("! WARNING: High coefficient of variation (>2) - high variability in counts")
  }

  low_pct <- stats$low_count_guides / stats$total_guides * 100
  if (low_pct > 20) {
    message(paste("! WARNING:", round(low_pct, 1), "% of guides have low counts"))
    message("           Consider filtering or using specialized low-count methods")
  }

  # Create visualizations
  message("\nGenerating visualizations...")

  plots <- list()
  plots$count_dist <- plot_count_distribution(counts_df)
  plots$log_dist <- plot_log_distribution(counts_df)
  plots$qc_summary <- plot_qc_summary(counts_df)
  plots$cumulative <- plot_cumulative(counts_df)
  plots$outliers <- plot_outliers(counts_df)

  if ("gene_name" %in% names(counts_df)) {
    plots$gene_boxplot <- plot_gene_boxplot(counts_df)
  }

  if ("sequence" %in% names(counts_df)) {
    plots$gc_bias <- plot_gc_bias(counts_df)
  }

  # Save combined plot
  message("Saving plots...")

  # Main QC report plot
  plot_list <- plots[!sapply(plots, is.null)]
  n_plots <- length(plot_list)

  if (n_plots >= 4) {
    combined <- gridExtra::arrangeGrob(grobs = plot_list[1:4], ncol = 2)
    ggsave(
      file.path(output_dir, "qc_report_main.png"),
      combined,
      width = 14,
      height = 12,
      dpi = 150
    )
  }

  # Save individual plots
  for (name in names(plots)) {
    if (!is.null(plots[[name]])) {
      ggsave(
        file.path(output_dir, paste0("qc_", name, ".png")),
        plots[[name]],
        width = 8,
        height = 6,
        dpi = 150
      )
    }
  }

  # Save flagged guides for review
  warning_guides <- counts_df %>%
    filter(qc_warning) %>%
    arrange(qc_status, count)

  write.csv(
    warning_guides,
    file.path(output_dir, "qc_warning_guides.csv"),
    row.names = FALSE
  )
  message(paste("  Saved", nrow(warning_guides), "warning guides to qc_warning_guides.csv"))

  # Save full annotated data
  write.csv(
    counts_df,
    file.path(output_dir, "counts_with_qc_flags.csv"),
    row.names = FALSE
  )
  message("  Saved annotated counts to counts_with_qc_flags.csv")

  # Save summary statistics
  stats_df <- data.frame(
    metric = names(stats),
    value = unlist(stats)
  )
  write.csv(stats_df, file.path(output_dir, "qc_statistics.csv"), row.names = FALSE)

  message("\n=== QC ANALYSIS COMPLETE ===")
  message(paste("Output saved to:", output_dir))

  return(list(
    data = counts_df,
    stats = stats,
    plots = plots,
    warnings = warning_guides
  ))
}

# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

if (!interactive()) {
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    message("Usage: Rscript qc_visualization.R <counts_file.csv> [output_dir]")
    message("\nExample:")
    message("  Rscript qc_visualization.R output/grna_counts.csv output/qc/")
    quit(status = 1)
  }

  counts_file <- args[1]
  output_dir <- if (length(args) >= 2) args[2] else NULL

  if (!file.exists(counts_file)) {
    stop(paste("File not found:", counts_file))
  }

  results <- run_qc_analysis(counts_file, output_dir)
}
