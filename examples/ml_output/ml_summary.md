# ML Marker Prediction Summary

## Model Performance (5-fold CV)

|model               | ROC_AUC| Sensitivity| Specificity|
|:-------------------|-------:|-----------:|-----------:|
|Logistic Regression |       1|        0.95|           1|
|Random Forest       |       1|        1.00|           1|
|GBM                 |       1|        1.00|           1|

## Top 10 Feature Importance

|                  |feature           | importance|
|:-----------------|:-----------------|----------:|
|quality_score     |quality_score     |  7.7885829|
|log_median_count  |log_median_count  |  2.1654307|
|consistency_score |consistency_score |  1.6870468|
|median_count      |median_count      |  1.6506505|
|cv_count          |cv_count          |  1.3127389|
|min_count         |min_count         |  0.7669364|
|range_count       |range_count       |  0.7144984|
|mean_gc           |mean_gc           |  0.6910176|
|log_mean_count    |log_mean_count    |  0.6127663|
|mean_count        |mean_count        |  0.5366914|

## Gene Recommendations Summary

- Highly Recommended: 24
- Recommended: 1
- Consider: 0
- Not Recommended: 26

## Top 15 Genes for Validation

|gene_name | ensemble_prob|recommendation     |
|:---------|-------------:|:------------------|
|DNMT3A    |     0.9972561|Highly Recommended |
|NF1       |     0.9670700|Highly Recommended |
|ARID1A    |     0.9574801|Highly Recommended |
|NPM1      |     0.9512786|Highly Recommended |
|HRAS      |     0.9478378|Highly Recommended |
|KIT       |     0.9476186|Highly Recommended |
|FBXW7     |     0.9327709|Highly Recommended |
|ZNRF3     |     0.9095106|Highly Recommended |
|TP53      |     0.9036213|Highly Recommended |
|SETD2     |     0.8862570|Highly Recommended |
|IDH1      |     0.8861161|Highly Recommended |
|BRCA1     |     0.8819793|Highly Recommended |
|ASXL1     |     0.8813215|Highly Recommended |
|MSH2      |     0.8767343|Highly Recommended |
|EZH2      |     0.8735201|Highly Recommended |
