# ML-Based Gene Marker Prediction

This document provides a step-by-step explanation of the machine learning pipeline used to predict which genes are good markers for further validation in CRISPR screens.

## Overview

The `ml_marker_prediction.R` script uses three different machine learning models to classify genes as "Good" or "Poor" markers based on their gRNA count data and derived quality metrics.

## Table of Contents

1. [Input Requirements](#input-requirements)
2. [Data Preparation](#data-preparation)
3. [Feature Engineering](#feature-engineering)
4. [Target Variable Definition](#target-variable-definition)
5. [Model 1: Logistic Regression](#model-1-logistic-regression)
6. [Model 2: Random Forest](#model-2-random-forest)
7. [Model 3: Gradient Boosting Machine](#model-3-gradient-boosting-machine)
8. [Ensemble Prediction](#ensemble-prediction)
9. [Output Files](#output-files)
10. [Usage](#usage)

---

## Input Requirements

### Input File Format

The script expects a CSV file with the following required columns:

| Column | Type | Description |
|--------|------|-------------|
| `guide_id` | string | Unique identifier for each gRNA |
| `gene_name` | string | Target gene name |
| `sequence` | string | gRNA sequence (e.g., "ACGTACGTACGTACGTACGT") |
| `count` | integer | Raw read count for the gRNA |

### Example Input

```csv
guide_id,gene_name,sequence,count
guide_001,TP53,ACGTACGTACGTACGTACGT,1523
guide_002,TP53,TGCATGCATGCATGCATGCA,1456
guide_003,BRCA1,GCTAGCTAGCTAGCTAGCTA,892
```

---

## Data Preparation

### Step 1: Load and Validate Data

```r
counts_df <- read_csv(counts_file)
```

The script:
1. Loads the CSV file
2. Validates that all required columns are present
3. Reports the number of guides and genes loaded

---

## Feature Engineering

### Step 2: Calculate Per-Guide Features

For each gRNA, the following features are calculated:

| Feature | Formula | Description |
|---------|---------|-------------|
| `log_count` | log10(count + 1) | Log-transformed count for normalization |
| `gc_content` | (G + C) / length | Proportion of G and C nucleotides |
| `seq_length` | nchar(sequence) | Length of the gRNA sequence |
| `freq_A` | A / length | Frequency of adenine |
| `freq_T` | T / length | Frequency of thymine |
| `freq_G` | G / length | Frequency of guanine |
| `freq_C` | C / length | Frequency of cytosine |

### Step 3: Calculate Gene-Level Statistics

Guides are aggregated by gene to create gene-level features:

| Feature | Formula | Description |
|---------|---------|-------------|
| `n_guides` | count of guides | Number of gRNAs per gene |
| `mean_count` | mean(count) | Average count across all guides |
| `median_count` | median(count) | Median count across all guides |
| `sd_count` | sd(count) | Standard deviation of counts |
| `cv_count` | sd / mean | Coefficient of variation (consistency measure) |
| `min_count` | min(count) | Minimum count |
| `max_count` | max(count) | Maximum count |
| `range_count` | max - min | Range of counts |
| `mean_gc` | mean(gc_content) | Average GC content |
| `sd_gc` | sd(gc_content) | Variability in GC content |
| `zero_count_guides` | sum(count == 0) | Number of guides with zero counts |
| `low_count_guides` | sum(count < 10) | Number of guides with low counts |
| `log_mean_count` | log10(mean + 1) | Log-transformed mean count |
| `log_median_count` | log10(median + 1) | Log-transformed median count |
| `consistency_score` | 1 / (cv + 0.1) | Higher score = more consistent |

### Step 4: Create Quality Flags

Binary flags are created to identify potential issues:

| Flag | Condition | Meaning |
|------|-----------|---------|
| `has_zero_counts` | zero_count_guides > 0 | Gene has guides with no reads |
| `all_low_counts` | low_count_guides == n_guides | All guides have low counts |
| `high_variability` | cv_count > 1.5 | High variability across guides |
| `extreme_gc` | mean_gc < 0.3 or > 0.7 | GC content outside optimal range |

---

## Target Variable Definition

### Step 5: Calculate Quality Score

A composite quality score (0-1) is calculated for each gene:

```
quality_score = (score_count × 0.25) +
                (score_consistency × 0.25) +
                (score_no_zeros × 0.20) +
                (score_gc × 0.15) +
                (score_variability × 0.15)
```

| Component | Weight | Description |
|-----------|--------|-------------|
| `score_count` | 25% | Normalized log mean count |
| `score_consistency` | 25% | Normalized consistency score |
| `score_no_zeros` | 20% | 1 if no zeros, 0 otherwise |
| `score_gc` | 15% | 1 if normal GC, 0.5 if extreme |
| `score_variability` | 15% | 1 if low CV, 0.5 if high |

### Step 6: Define Binary Target

Genes are classified as:
- **Good**: quality_score >= median (top 50%)
- **Poor**: quality_score < median (bottom 50%)

Non-targeting controls are always classified as "Poor".

---

## Model 1: Logistic Regression

### Algorithm: glmnet (Elastic Net Regularization)

Logistic regression models the probability of a gene being a "Good" marker using a linear combination of features.

### How It Works

1. **Linear Combination**: Calculates weighted sum of features
   ```
   z = β₀ + β₁×feature₁ + β₂×feature₂ + ... + βₙ×featureₙ
   ```

2. **Sigmoid Function**: Converts to probability
   ```
   P(Good) = 1 / (1 + e^(-z))
   ```

3. **Regularization**: Prevents overfitting using elastic net
   - **L1 (Lasso)**: Encourages sparse solutions (some coefficients = 0)
   - **L2 (Ridge)**: Shrinks coefficients toward zero
   - **Alpha parameter**: Controls L1/L2 balance (0=Ridge, 1=Lasso, 0.5=Elastic Net)

### Hyperparameters Tuned

| Parameter | Values Tested | Description |
|-----------|---------------|-------------|
| `alpha` | 0, 0.5, 1 | Regularization type |
| `lambda` | 0.001, 0.01, 0.1 | Regularization strength |

### Strengths
- Interpretable coefficients
- Fast training
- Works well with fewer samples
- Provides probability estimates

### Limitations
- Assumes linear relationship between features and log-odds
- May underperform with complex feature interactions

---

## Model 2: Random Forest

### Algorithm: ranger (Fast Random Forest Implementation)

Random Forest is an ensemble method that builds multiple decision trees and combines their predictions.

### How It Works

1. **Bootstrap Sampling**: Create N random samples from training data (with replacement)

2. **Tree Building**: For each sample:
   - At each node, randomly select `mtry` features
   - Find the best split among selected features
   - Continue until minimum node size reached

3. **Prediction**:
   - Each tree votes for a class
   - Final prediction = majority vote
   - Probability = proportion of trees voting "Good"

```
                    [Data]
                       |
        ┌──────────────┼──────────────┐
        v              v              v
     [Tree 1]       [Tree 2]   ... [Tree 500]
        |              |              |
     [Vote]         [Vote]         [Vote]
        └──────────────┼──────────────┘
                       v
               [Majority Vote]
```

### Hyperparameters Tuned

| Parameter | Values Tested | Description |
|-----------|---------------|-------------|
| `mtry` | 2, 4, 6 | Number of features to consider at each split |
| `splitrule` | gini | Splitting criterion |
| `min.node.size` | 1, 3, 5 | Minimum samples in terminal node |

### Feature Importance

Importance is calculated using **Gini impurity decrease**:
- Measures how much each feature reduces classification error
- Higher values = more important feature

### Strengths
- Handles non-linear relationships
- Robust to outliers
- Provides feature importance
- Less prone to overfitting

### Limitations
- Less interpretable than logistic regression
- Can be slow with very large datasets

---

## Model 3: Gradient Boosting Machine

### Algorithm: gbm (Gradient Boosting)

GBM builds trees sequentially, with each tree correcting errors from previous trees.

### How It Works

1. **Initialize**: Start with a simple prediction (e.g., mean)

2. **Iterative Boosting**:
   ```
   For each iteration t:
     1. Calculate residuals (errors) from current model
     2. Fit a new tree to predict residuals
     3. Add new tree to model (with shrinkage)

   Final Model = Tree₁ + η×Tree₂ + η×Tree₃ + ... + η×Treeₙ
   ```

3. **Shrinkage (Learning Rate)**: Each tree's contribution is scaled by η (0.1)

```
Iteration 1:  [Initial] ──────────────────> [Model₁]
                                               |
Iteration 2:  [Residuals₁] ──> [Tree₂] ──> [Model₂]
                                               |
Iteration 3:  [Residuals₂] ──> [Tree₃] ──> [Model₃]
                                               |
              ...                              ...
```

### Hyperparameters Tuned

| Parameter | Values Tested | Description |
|-----------|---------------|-------------|
| `n.trees` | 100, 200 | Number of boosting iterations |
| `interaction.depth` | 2, 4 | Maximum tree depth |
| `shrinkage` | 0.1 | Learning rate |
| `n.minobsinnode` | 5 | Minimum samples in terminal node |

### Strengths
- Often achieves best predictive performance
- Handles complex feature interactions
- Can model non-linear relationships

### Limitations
- Prone to overfitting if not tuned properly
- Slower training than Random Forest
- Less interpretable

---

## Ensemble Prediction

### Step 7: Combine Model Predictions

The final prediction combines all three models:

```
ensemble_prob = (lr_prob + rf_prob + gbm_prob) / 3
```

### Recommendation Categories

| Category | Threshold | Interpretation |
|----------|-----------|----------------|
| Highly Recommended | ≥ 0.70 | Strong candidate for validation |
| Recommended | 0.50 - 0.69 | Good candidate |
| Consider | 0.30 - 0.49 | May be worth investigating |
| Not Recommended | < 0.30 | Poor marker characteristics |

---

## Output Files

### CSV Files

| File | Description |
|------|-------------|
| `ml_model_comparison.csv` | Cross-validation metrics for all models |
| `ml_feature_importance.csv` | Feature importance scores from Random Forest |
| `ml_gene_predictions.csv` | Predictions for all genes |

### ml_gene_predictions.csv Columns

| Column | Description |
|--------|-------------|
| `gene_name` | Gene identifier |
| `quality_score` | Rule-based quality score (0-1) |
| `lr_prob` | Logistic Regression probability |
| `rf_prob` | Random Forest probability |
| `gbm_prob` | GBM probability |
| `ensemble_prob` | Average of all three models |
| `recommendation` | Categorical recommendation |
| `marker_class` | Actual class (Good/Poor) |

### Interactive HTML Visualizations

| File | Description |
|------|-------------|
| `ml_model_comparison.html` | Bar chart comparing model performance |
| `ml_feature_importance.html` | Horizontal bar chart of feature importance |
| `ml_model_agreement.html` | Scatter plot: RF vs GBM predictions |
| `ml_top_genes.html` | Ranked bar chart of top genes |
| `ml_score_distribution.html` | Histogram of ensemble scores |
| `ml_quality_vs_prediction.html` | Quality score vs ML prediction |

### Summary Report

| File | Description |
|------|-------------|
| `ml_summary.md` | Markdown summary with tables |

---

## Usage

### Basic Usage

```bash
Rscript ml_marker_prediction.R <input_file> <output_dir>
```

### Examples

```bash
# Using default paths (output/grna_counts.csv -> output/)
Rscript ml_marker_prediction.R

# Specify input file
Rscript ml_marker_prediction.R my_counts.csv

# Specify input and output
Rscript ml_marker_prediction.R my_counts.csv results/
```

### Required R Packages

| Package | Purpose |
|---------|---------|
| caret | ML framework and cross-validation |
| ranger | Fast Random Forest |
| gbm | Gradient Boosting Machine |
| glmnet | Regularized Logistic Regression |
| dplyr | Data manipulation |
| tidyr | Data reshaping |
| ggplot2 | Static plots |
| plotly | Interactive plots |
| htmlwidgets | Save HTML widgets |
| readr | CSV reading/writing |

### Installation (Ubuntu/Debian)

```bash
sudo apt-get install r-cran-caret r-cran-ranger r-cran-gbm r-cran-glmnet \
                     r-cran-dplyr r-cran-tidyr r-cran-ggplot2 r-cran-plotly \
                     r-cran-htmlwidgets r-cran-readr
```

---

## Interpreting Results

### Model Performance Metrics

| Metric | Description | Good Value |
|--------|-------------|------------|
| ROC AUC | Area under ROC curve | > 0.8 |
| Sensitivity | True positive rate | > 0.8 |
| Specificity | True negative rate | > 0.8 |

### Feature Importance Interpretation

Higher importance means the feature is more useful for distinguishing Good from Poor markers:

1. **quality_score**: Composite quality metric (most important)
2. **log_median_count**: Count levels matter
3. **consistency_score**: Consistent guides = better marker
4. **cv_count**: Lower variability = better marker

### Recommendations for Validation

1. Start with "Highly Recommended" genes
2. Check ensemble probability and individual model agreement
3. Review raw count data for top candidates
4. Consider biological relevance alongside ML predictions

---

## Workflow Diagram

```
┌─────────────────────────────────────────────────────────────────┐
│                         INPUT                                    │
│                    grna_counts.csv                               │
│            (guide_id, gene_name, sequence, count)                │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    FEATURE ENGINEERING                           │
│  ┌─────────────────┐    ┌──────────────────────────────────┐    │
│  │ Per-Guide       │    │ Gene-Level Aggregation            │    │
│  │ - GC content    │───>│ - mean/median/sd counts           │    │
│  │ - Nucleotide %  │    │ - CV, consistency score           │    │
│  │ - Log count     │    │ - Quality flags                   │    │
│  └─────────────────┘    └──────────────────────────────────┘    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    TARGET DEFINITION                             │
│           Quality Score → Binary Classification                  │
│                    (Good vs Poor markers)                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    MODEL TRAINING                                │
│  ┌─────────────┐  ┌─────────────┐  ┌─────────────┐              │
│  │  Logistic   │  │   Random    │  │    GBM      │              │
│  │ Regression  │  │   Forest    │  │  Boosting   │              │
│  │  (glmnet)   │  │  (ranger)   │  │             │              │
│  └──────┬──────┘  └──────┬──────┘  └──────┬──────┘              │
│         │                │                │                      │
│         └────────────────┼────────────────┘                      │
│                          ▼                                       │
│                   5-Fold Cross-Validation                        │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    ENSEMBLE PREDICTION                           │
│         ensemble_prob = (LR + RF + GBM) / 3                      │
│                          │                                       │
│         ┌────────────────┼────────────────┐                      │
│         ▼                ▼                ▼                      │
│   Highly Recommended  Recommended    Not Recommended             │
│      (≥0.70)         (0.50-0.69)       (<0.30)                   │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                         OUTPUT                                   │
│  ┌─────────────────┐  ┌─────────────────┐  ┌─────────────────┐  │
│  │ CSV Files       │  │ HTML Plots      │  │ Summary         │  │
│  │ - predictions   │  │ - interactive   │  │ - ml_summary.md │  │
│  │ - importance    │  │ - hover info    │  │                 │  │
│  │ - comparison    │  │ - gene names    │  │                 │  │
│  └─────────────────┘  └─────────────────┘  └─────────────────┘  │
└─────────────────────────────────────────────────────────────────┘
```
