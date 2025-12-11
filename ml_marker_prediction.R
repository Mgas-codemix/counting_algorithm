#!/usr/bin/env Rscript
# ==============================================================================
# ML-based Gene Marker Prediction for CRISPR Screen QC
# ==============================================================================
# This script uses caret to predict which genes are good markers for
# further validation based on gRNA count data and QC metrics.
#
# Models used:
# 1. Random Forest (ranger)
# 2. Logistic Regression (glmnet)
# 3. Gradient Boosting Machine (gbm) - similar to XGBoost
#
# Usage: Rscript ml_marker_prediction.R <input_counts.csv> <output_dir>
# ==============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(caret)
  library(randomForest)
  library(ranger)
  library(gbm)
  library(glmnet)
  library(e1071)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readr)
  library(plotly)
  library(htmlwidgets)
})

# Set seed for reproducibility
set.seed(42)

# ==============================================================================
# Custom theme for plots
# ==============================================================================
theme_clean <- function() {
  theme_minimal() +
    theme(
      plot.background = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA),
      panel.grid.major = element_line(color = "gray90", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.3),
      axis.text = element_text(color = "black", size = 10),
      axis.title = element_text(color = "black", size = 11, face = "bold"),
      plot.title = element_text(color = "black", size = 14, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(color = "gray30", size = 10, hjust = 0.5),
      legend.background = element_rect(fill = "white", color = NA),
      legend.key = element_rect(fill = "white", color = NA),
      strip.background = element_rect(fill = "gray95", color = "black"),
      strip.text = element_text(color = "black", face = "bold")
    )
}

# ==============================================================================
# Data Loading and Feature Engineering
# ==============================================================================
prepare_ml_data <- function(counts_file) {
  cat("Loading data from:", counts_file, "\n")

  # Load count data
  counts_df <- read_csv(counts_file, show_col_types = FALSE)

  # Ensure required columns exist
  required_cols <- c("guide_id", "gene_name", "sequence", "count")
  if (!all(required_cols %in% colnames(counts_df))) {
    stop("Missing required columns. Need: ", paste(required_cols, collapse = ", "))
  }

  cat("Loaded", nrow(counts_df), "guides for", length(unique(counts_df$gene_name)), "genes\n")

  # Calculate per-guide features
  counts_df <- counts_df %>%
    mutate(
      # Log transform count
      log_count = log10(count + 1),

      # GC content
      gc_content = sapply(sequence, function(s) {
        chars <- strsplit(toupper(s), "")[[1]]
        sum(chars %in% c("G", "C")) / length(chars)
      }),

      # Sequence length
      seq_length = nchar(sequence),

      # Nucleotide frequencies
      freq_A = sapply(sequence, function(s) {
        chars <- strsplit(toupper(s), "")[[1]]
        sum(chars == "A") / length(chars)
      }),
      freq_T = sapply(sequence, function(s) {
        chars <- strsplit(toupper(s), "")[[1]]
        sum(chars == "T") / length(chars)
      }),
      freq_G = sapply(sequence, function(s) {
        chars <- strsplit(toupper(s), "")[[1]]
        sum(chars == "G") / length(chars)
      }),
      freq_C = sapply(sequence, function(s) {
        chars <- strsplit(toupper(s), "")[[1]]
        sum(chars == "C") / length(chars)
      })
    )

  # Calculate gene-level statistics
  gene_stats <- counts_df %>%
    group_by(gene_name) %>%
    summarise(
      n_guides = n(),
      mean_count = mean(count, na.rm = TRUE),
      median_count = median(count, na.rm = TRUE),
      sd_count = sd(count, na.rm = TRUE),
      cv_count = ifelse(mean_count > 0, sd_count / mean_count, NA),
      min_count = min(count, na.rm = TRUE),
      max_count = max(count, na.rm = TRUE),
      range_count = max_count - min_count,
      mean_gc = mean(gc_content, na.rm = TRUE),
      sd_gc = sd(gc_content, na.rm = TRUE),
      zero_count_guides = sum(count == 0),
      low_count_guides = sum(count < 10),
      .groups = "drop"
    ) %>%
    mutate(
      # Log transform gene-level stats
      log_mean_count = log10(mean_count + 1),
      log_median_count = log10(median_count + 1),

      # Consistency score (lower CV = more consistent)
      consistency_score = 1 / (cv_count + 0.1),

      # Replace NA/Inf values
      across(where(is.numeric), ~ifelse(is.na(.) | is.infinite(.), 0, .))
    )

  # Calculate global statistics for outlier detection
  global_median <- median(counts_df$count, na.rm = TRUE)
  global_mad <- mad(counts_df$count, na.rm = TRUE)
  q1 <- quantile(counts_df$count, 0.25, na.rm = TRUE)
  q3 <- quantile(counts_df$count, 0.75, na.rm = TRUE)
  iqr <- q3 - q1

  # Add outlier flags to gene stats
  gene_stats <- gene_stats %>%
    mutate(
      # Flag genes with any zero counts
      has_zero_counts = zero_count_guides > 0,

      # Flag genes with all low counts
      all_low_counts = low_count_guides == n_guides,

      # Flag genes with extreme CV (high variability)
      high_variability = cv_count > 1.5,

      # Flag genes with extreme GC content
      extreme_gc = mean_gc < 0.3 | mean_gc > 0.7
    )

  return(list(
    guide_data = counts_df,
    gene_data = gene_stats,
    global_stats = list(
      median = global_median,
      mad = global_mad,
      q1 = q1,
      q3 = q3,
      iqr = iqr
    )
  ))
}

# ==============================================================================
# Define Target Variable: Good Marker Gene
# ==============================================================================
define_target <- function(gene_data) {
  # A "good marker" gene has:
  # 1. Moderate to high counts (detectable signal)
  # 2. Low variability across guides (consistent)
  # 3. No zero count guides
  # 4. Reasonable GC content
  # 5. Not an extreme outlier

  gene_data <- gene_data %>%
    mutate(
      # Score components (0-1 scale)
      score_count = pmin(1, log_mean_count / max(log_mean_count, na.rm = TRUE)),
      score_consistency = pmin(1, consistency_score / max(consistency_score, na.rm = TRUE)),
      score_no_zeros = ifelse(has_zero_counts, 0, 1),
      score_gc = ifelse(extreme_gc, 0.5, 1),
      score_variability = ifelse(high_variability, 0.5, 1),

      # Combined quality score
      quality_score = (score_count * 0.25 +
                        score_consistency * 0.25 +
                        score_no_zeros * 0.2 +
                        score_gc * 0.15 +
                        score_variability * 0.15),

      # Binary classification: good marker (top 50% by quality score)
      # Exclude non-targeting controls
      is_good_marker = ifelse(
        grepl("Non-targeting|NTC", gene_name, ignore.case = TRUE),
        FALSE,
        quality_score >= median(quality_score[!grepl("Non-targeting|NTC", gene_name, ignore.case = TRUE)], na.rm = TRUE)
      ),

      # Convert to factor for classification
      marker_class = factor(
        ifelse(is_good_marker, "Good", "Poor"),
        levels = c("Poor", "Good")
      )
    )

  cat("\nTarget variable distribution:\n")
  print(table(gene_data$marker_class))

  return(gene_data)
}

# ==============================================================================
# Train ML Models using caret
# ==============================================================================
train_models <- function(gene_data, output_dir) {
  cat("\n========================================\n")
  cat("Training Machine Learning Models\n")
  cat("========================================\n")

  # Prepare features for modeling
  model_data <- gene_data %>%
    filter(!is.na(marker_class)) %>%
    select(
      gene_name,
      marker_class,
      mean_count, log_mean_count, median_count, log_median_count,
      sd_count, cv_count, min_count, max_count, range_count,
      mean_gc, sd_gc, n_guides,
      consistency_score, quality_score
    ) %>%
    # Remove any remaining NA/Inf
    mutate(across(where(is.numeric), ~ifelse(is.na(.) | is.infinite(.), 0, .)))

  # Split data
  set.seed(42)
  train_idx <- createDataPartition(model_data$marker_class, p = 0.75, list = FALSE)
  train_data <- model_data[train_idx, ]
  test_data <- model_data[-train_idx, ]

  cat("\nTraining set:", nrow(train_data), "genes\n")
  cat("Test set:", nrow(test_data), "genes\n")

  # Define cross-validation
  ctrl <- trainControl(
    method = "cv",
    number = 5,
    classProbs = TRUE,
    summaryFunction = twoClassSummary,
    savePredictions = TRUE
  )

  # Features for modeling (exclude gene_name and marker_class)
  feature_cols <- c("mean_count", "log_mean_count", "median_count", "log_median_count",
                    "sd_count", "cv_count", "min_count", "max_count", "range_count",
                    "mean_gc", "sd_gc", "n_guides", "consistency_score", "quality_score")

  train_features <- train_data[, feature_cols]
  train_labels <- train_data$marker_class
  test_features <- test_data[, feature_cols]
  test_labels <- test_data$marker_class

  # ==============================================================================
  # Model 1: Logistic Regression (glmnet)
  # ==============================================================================
  cat("\n--- Training Logistic Regression ---\n")

  lr_model <- train(
    x = train_features,
    y = train_labels,
    method = "glmnet",
    trControl = ctrl,
    metric = "ROC",
    tuneGrid = expand.grid(
      alpha = c(0, 0.5, 1),  # Ridge, Elastic Net, Lasso
      lambda = c(0.001, 0.01, 0.1)
    ),
    preProcess = c("center", "scale")
  )

  cat("Logistic Regression CV Results:\n")
  print(lr_model$results[which.max(lr_model$results$ROC), ])

  # ==============================================================================
  # Model 2: Random Forest (ranger)
  # ==============================================================================
  cat("\n--- Training Random Forest ---\n")

  rf_model <- train(
    x = train_features,
    y = train_labels,
    method = "ranger",
    trControl = ctrl,
    metric = "ROC",
    tuneGrid = expand.grid(
      mtry = c(2, 4, 6),
      splitrule = "gini",
      min.node.size = c(1, 3, 5)
    ),
    importance = "impurity",
    preProcess = c("center", "scale")
  )

  cat("Random Forest CV Results:\n")
  print(rf_model$results[which.max(rf_model$results$ROC), ])

  # ==============================================================================
  # Model 3: Gradient Boosting Machine (gbm)
  # ==============================================================================
  cat("\n--- Training Gradient Boosting Machine ---\n")

  gbm_model <- train(
    x = train_features,
    y = train_labels,
    method = "gbm",
    trControl = ctrl,
    metric = "ROC",
    tuneGrid = expand.grid(
      n.trees = c(100, 200),
      interaction.depth = c(2, 4),
      shrinkage = 0.1,
      n.minobsinnode = 5
    ),
    preProcess = c("center", "scale"),
    verbose = FALSE
  )

  cat("GBM CV Results:\n")
  print(gbm_model$results[which.max(gbm_model$results$ROC), ])

  # ==============================================================================
  # Compare Models
  # ==============================================================================
  cat("\n========================================\n")
  cat("Model Comparison\n")
  cat("========================================\n")

  # Get best results for each model
  lr_best <- lr_model$results[which.max(lr_model$results$ROC), ]
  rf_best <- rf_model$results[which.max(rf_model$results$ROC), ]
  gbm_best <- gbm_model$results[which.max(gbm_model$results$ROC), ]

  comparison <- data.frame(
    model = c("Logistic Regression", "Random Forest", "GBM"),
    ROC_AUC = c(lr_best$ROC, rf_best$ROC, gbm_best$ROC),
    Sensitivity = c(lr_best$Sens, rf_best$Sens, gbm_best$Sens),
    Specificity = c(lr_best$Spec, rf_best$Spec, gbm_best$Spec)
  )

  print(comparison)

  # Save comparison
  write_csv(comparison, file.path(output_dir, "ml_model_comparison.csv"))

  # ==============================================================================
  # Test Set Predictions
  # ==============================================================================
  cat("\n--- Test Set Performance ---\n")

  # Predictions on test set
  lr_pred <- predict(lr_model, test_features, type = "prob")
  rf_pred <- predict(rf_model, test_features, type = "prob")
  gbm_pred <- predict(gbm_model, test_features, type = "prob")

  lr_class <- predict(lr_model, test_features)
  rf_class <- predict(rf_model, test_features)
  gbm_class <- predict(gbm_model, test_features)

  # Confusion matrices
  cat("\nLogistic Regression Confusion Matrix:\n")
  print(confusionMatrix(lr_class, test_labels))

  cat("\nRandom Forest Confusion Matrix:\n")
  print(confusionMatrix(rf_class, test_labels))

  cat("\nGBM Confusion Matrix:\n")
  print(confusionMatrix(gbm_class, test_labels))

  # ==============================================================================
  # Feature Importance (from Random Forest)
  # ==============================================================================
  cat("\n--- Feature Importance (Random Forest) ---\n")

  # Extract variable importance from ranger model
  rf_final <- rf_model$finalModel
  importance_df <- data.frame(
    feature = names(rf_final$variable.importance),
    importance = rf_final$variable.importance
  ) %>%
    arrange(desc(importance))

  print(importance_df)
  write_csv(importance_df, file.path(output_dir, "ml_feature_importance.csv"))

  # ==============================================================================
  # Predict on All Genes
  # ==============================================================================
  cat("\n--- Generating Predictions for All Genes ---\n")

  all_features <- model_data[, feature_cols]

  all_predictions <- model_data %>%
    mutate(
      lr_prob = predict(lr_model, all_features, type = "prob")$Good,
      rf_prob = predict(rf_model, all_features, type = "prob")$Good,
      gbm_prob = predict(gbm_model, all_features, type = "prob")$Good,

      # Ensemble prediction (average of all models)
      ensemble_prob = (lr_prob + rf_prob + gbm_prob) / 3,

      # Final recommendation
      recommendation = case_when(
        ensemble_prob >= 0.7 ~ "Highly Recommended",
        ensemble_prob >= 0.5 ~ "Recommended",
        ensemble_prob >= 0.3 ~ "Consider",
        TRUE ~ "Not Recommended"
      )
    ) %>%
    arrange(desc(ensemble_prob))

  # Save predictions
  predictions_output <- all_predictions %>%
    select(gene_name, quality_score, lr_prob, rf_prob, gbm_prob,
           ensemble_prob, recommendation, marker_class)

  write_csv(predictions_output, file.path(output_dir, "ml_gene_predictions.csv"))

  cat("\nTop 10 Recommended Genes:\n")
  print(predictions_output %>% head(10))

  return(list(
    models = list(lr = lr_model, rf = rf_model, gbm = gbm_model),
    comparison = comparison,
    importance = importance_df,
    predictions = all_predictions,
    test_data = test_data
  ))
}

# ==============================================================================
# Generate Visualizations
# ==============================================================================
generate_ml_plots <- function(results, gene_data, output_dir) {
  cat("\n========================================\n")
  cat("Generating ML Visualizations\n")
  cat("========================================\n")

  predictions <- results$predictions
  importance <- results$importance
  comparison <- results$comparison

  # ---------------------------------------------------------------------------
  # Plot 1: Model Comparison (ROC AUC)
  # ---------------------------------------------------------------------------
  comparison_long <- comparison %>%
    pivot_longer(-model, names_to = "metric", values_to = "value")

  p1 <- ggplot(comparison_long, aes(x = model, y = value, fill = model)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~metric, scales = "free_y") +
    scale_fill_manual(values = c(
      "Logistic Regression" = "#3498db",
      "Random Forest" = "#2ecc71",
      "GBM" = "#e74c3c"
    )) +
    labs(
      title = "Model Performance Comparison",
      subtitle = "5-fold Cross-validation Results",
      x = "", y = "Score"
    ) +
    theme_clean() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")

  ggsave(file.path(output_dir, "ml_model_comparison.png"), p1,
         width = 10, height = 6, dpi = 150, bg = "white")

  # Interactive version
  p1_plotly <- ggplotly(p1) %>%
    layout(paper_bgcolor = "white", plot_bgcolor = "white")
  saveWidget(p1_plotly, file.path(output_dir, "ml_model_comparison.html"),
             selfcontained = TRUE)

  # ---------------------------------------------------------------------------
  # Plot 2: Feature Importance
  # ---------------------------------------------------------------------------
  p2 <- ggplot(importance %>% head(10),
               aes(x = reorder(feature, importance), y = importance, fill = importance)) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_gradient(low = "#3498db", high = "#e74c3c") +
    labs(
      title = "Feature Importance (Random Forest)",
      subtitle = "Top 10 Most Important Features",
      x = "", y = "Importance Score"
    ) +
    theme_clean() +
    theme(legend.position = "none")

  ggsave(file.path(output_dir, "ml_feature_importance.png"), p2,
         width = 8, height = 6, dpi = 150, bg = "white")

  # Interactive version
  p2_data <- importance %>% head(10)
  p2_data$hover_text <- paste0(
    "<b>Feature:</b> ", p2_data$feature,
    "<br><b>Importance:</b> ", round(p2_data$importance, 2)
  )

  p2_plotly <- plot_ly(p2_data,
                       y = ~reorder(feature, importance),
                       x = ~importance,
                       type = "bar",
                       orientation = "h",
                       text = ~hover_text,
                       hoverinfo = "text",
                       marker = list(
                         color = ~importance,
                         colorscale = list(c(0, "#3498db"), c(1, "#e74c3c"))
                       )) %>%
    layout(
      title = list(text = "Feature Importance (Random Forest)"),
      xaxis = list(title = "Importance Score", showline = TRUE, linecolor = "black"),
      yaxis = list(title = "", showline = TRUE, linecolor = "black"),
      paper_bgcolor = "white",
      plot_bgcolor = "white"
    )

  saveWidget(p2_plotly, file.path(output_dir, "ml_feature_importance.html"),
             selfcontained = TRUE)

  # ---------------------------------------------------------------------------
  # Plot 3: Gene Prediction Probabilities
  # ---------------------------------------------------------------------------
  predictions$hover_text <- paste0(
    "<b>Gene:</b> ", predictions$gene_name,
    "<br><b>LR Prob:</b> ", round(predictions$lr_prob, 3),
    "<br><b>RF Prob:</b> ", round(predictions$rf_prob, 3),
    "<br><b>GBM Prob:</b> ", round(predictions$gbm_prob, 3),
    "<br><b>Ensemble:</b> ", round(predictions$ensemble_prob, 3),
    "<br><b>Recommendation:</b> ", predictions$recommendation
  )

  # Scatter plot of RF vs GBM probabilities
  p3 <- plot_ly(predictions,
                x = ~rf_prob,
                y = ~gbm_prob,
                type = "scatter",
                mode = "markers",
                color = ~recommendation,
                colors = c(
                  "Highly Recommended" = "#27ae60",
                  "Recommended" = "#3498db",
                  "Consider" = "#f39c12",
                  "Not Recommended" = "#e74c3c"
                ),
                text = ~hover_text,
                hoverinfo = "text",
                marker = list(size = 10, opacity = 0.7)) %>%
    layout(
      title = list(text = "Model Agreement: Random Forest vs GBM"),
      xaxis = list(title = "Random Forest Probability",
                   showline = TRUE, linecolor = "black", linewidth = 2,
                   range = c(0, 1)),
      yaxis = list(title = "GBM Probability",
                   showline = TRUE, linecolor = "black", linewidth = 2,
                   range = c(0, 1)),
      paper_bgcolor = "white",
      plot_bgcolor = "white",
      shapes = list(
        list(type = "line", x0 = 0, x1 = 1, y0 = 0, y1 = 1,
             line = list(color = "gray", dash = "dot"))
      )
    )

  saveWidget(p3, file.path(output_dir, "ml_model_agreement.html"),
             selfcontained = TRUE)

  # ---------------------------------------------------------------------------
  # Plot 4: Gene Ranking by Ensemble Score
  # ---------------------------------------------------------------------------
  top_genes <- predictions %>%
    filter(!grepl("Non-targeting|NTC", gene_name, ignore.case = TRUE)) %>%
    head(20)

  top_genes$hover_text <- paste0(
    "<b>Gene:</b> ", top_genes$gene_name,
    "<br><b>Ensemble Score:</b> ", round(top_genes$ensemble_prob, 3),
    "<br><b>Quality Score:</b> ", round(top_genes$quality_score, 3),
    "<br><b>Recommendation:</b> ", top_genes$recommendation
  )

  p4 <- plot_ly(top_genes,
                y = ~reorder(gene_name, ensemble_prob),
                x = ~ensemble_prob,
                type = "bar",
                orientation = "h",
                color = ~recommendation,
                colors = c(
                  "Highly Recommended" = "#27ae60",
                  "Recommended" = "#3498db",
                  "Consider" = "#f39c12",
                  "Not Recommended" = "#e74c3c"
                ),
                text = ~hover_text,
                hoverinfo = "text") %>%
    layout(
      title = list(text = "Top 20 Genes for Validation"),
      xaxis = list(title = "Ensemble Probability (Good Marker)",
                   showline = TRUE, linecolor = "black", linewidth = 2,
                   range = c(0, 1)),
      yaxis = list(title = "", showline = TRUE, linecolor = "black", linewidth = 2),
      paper_bgcolor = "white",
      plot_bgcolor = "white",
      showlegend = TRUE
    )

  saveWidget(p4, file.path(output_dir, "ml_top_genes.html"),
             selfcontained = TRUE)

  # ---------------------------------------------------------------------------
  # Plot 5: Ensemble Score Distribution
  # ---------------------------------------------------------------------------
  p5 <- plot_ly(predictions,
                x = ~ensemble_prob,
                type = "histogram",
                nbinsx = 30,
                marker = list(
                  color = "#3498db",
                  line = list(color = "black", width = 1)
                )) %>%
    layout(
      title = list(text = "Distribution of Ensemble Prediction Scores"),
      xaxis = list(title = "Ensemble Probability (Good Marker)",
                   showline = TRUE, linecolor = "black", linewidth = 2),
      yaxis = list(title = "Count",
                   showline = TRUE, linecolor = "black", linewidth = 2),
      paper_bgcolor = "white",
      plot_bgcolor = "white",
      shapes = list(
        list(type = "line", x0 = 0.5, x1 = 0.5, y0 = 0, y1 = 1, yref = "paper",
             line = list(color = "red", dash = "dash", width = 2))
      )
    )

  saveWidget(p5, file.path(output_dir, "ml_score_distribution.html"),
             selfcontained = TRUE)

  # ---------------------------------------------------------------------------
  # Plot 6: Quality Score vs Ensemble Prediction
  # ---------------------------------------------------------------------------
  p6 <- plot_ly(predictions,
                x = ~quality_score,
                y = ~ensemble_prob,
                type = "scatter",
                mode = "markers",
                color = ~recommendation,
                colors = c(
                  "Highly Recommended" = "#27ae60",
                  "Recommended" = "#3498db",
                  "Consider" = "#f39c12",
                  "Not Recommended" = "#e74c3c"
                ),
                text = ~hover_text,
                hoverinfo = "text",
                marker = list(size = 10, opacity = 0.7)) %>%
    layout(
      title = list(text = "Quality Score vs ML Prediction"),
      xaxis = list(title = "Quality Score (Rule-based)",
                   showline = TRUE, linecolor = "black", linewidth = 2),
      yaxis = list(title = "Ensemble Probability (ML)",
                   showline = TRUE, linecolor = "black", linewidth = 2),
      paper_bgcolor = "white",
      plot_bgcolor = "white"
    )

  saveWidget(p6, file.path(output_dir, "ml_quality_vs_prediction.html"),
             selfcontained = TRUE)

  cat("All visualizations saved to:", output_dir, "\n")
}

# ==============================================================================
# Generate Summary Report
# ==============================================================================
generate_summary <- function(results, gene_data, output_dir) {
  predictions <- results$predictions

  summary_text <- c(
    "# ML Marker Prediction Summary",
    "",
    "## Model Performance (5-fold CV)",
    "",
    knitr::kable(results$comparison, format = "markdown"),
    "",
    "## Top 10 Feature Importance",
    "",
    knitr::kable(results$importance %>% head(10), format = "markdown"),
    "",
    "## Gene Recommendations Summary",
    "",
    paste("- Highly Recommended:", sum(predictions$recommendation == "Highly Recommended")),
    paste("- Recommended:", sum(predictions$recommendation == "Recommended")),
    paste("- Consider:", sum(predictions$recommendation == "Consider")),
    paste("- Not Recommended:", sum(predictions$recommendation == "Not Recommended")),
    "",
    "## Top 15 Genes for Validation",
    "",
    knitr::kable(
      predictions %>%
        filter(!grepl("Non-targeting|NTC", gene_name, ignore.case = TRUE)) %>%
        head(15) %>%
        select(gene_name, ensemble_prob, recommendation),
      format = "markdown"
    )
  )

  writeLines(summary_text, file.path(output_dir, "ml_summary.md"))
  cat("\nSummary report saved to:", file.path(output_dir, "ml_summary.md"), "\n")
}

# ==============================================================================
# Main Function
# ==============================================================================
main <- function() {
  # Parse command line arguments
  args <- commandArgs(trailingOnly = TRUE)

  if (length(args) < 1) {
    input_file <- "output/grna_counts.csv"
    output_dir <- "output"
  } else {
    input_file <- args[1]
    output_dir <- ifelse(length(args) >= 2, args[2], "output")
  }

  # Create output directory if needed
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  cat("========================================\n")
  cat("ML-based Gene Marker Prediction\n")
  cat("========================================\n")
  cat("Input file:", input_file, "\n")
  cat("Output directory:", output_dir, "\n")

  # Step 1: Prepare data
  cat("\n--- Step 1: Preparing Data ---\n")
  data <- prepare_ml_data(input_file)

  # Step 2: Define target variable
  cat("\n--- Step 2: Defining Target Variable ---\n")
  gene_data <- define_target(data$gene_data)

  # Step 3: Train models
  cat("\n--- Step 3: Training ML Models ---\n")
  results <- train_models(gene_data, output_dir)

  # Step 4: Generate visualizations
  cat("\n--- Step 4: Generating Visualizations ---\n")
  generate_ml_plots(results, gene_data, output_dir)

  # Step 5: Generate summary
  cat("\n--- Step 5: Generating Summary Report ---\n")
  generate_summary(results, gene_data, output_dir)

  cat("\n========================================\n")
  cat("ML Analysis Complete!\n")
  cat("========================================\n")
  cat("\nOutput files:\n")
  cat("- ml_model_comparison.csv: Cross-validation metrics\n")
  cat("- ml_feature_importance.csv: Feature importance scores\n")
  cat("- ml_gene_predictions.csv: All gene predictions\n")
  cat("- ml_summary.md: Summary report\n")
  cat("- Interactive HTML plots for exploration\n")
}

# Run main function
main()
