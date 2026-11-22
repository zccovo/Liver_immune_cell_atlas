# Liver Immune Cell Atlas

This repository provides R scripts and differential-expression results for our single-cell RNA sequencing (scRNA-seq) study of hepatic immune cells in dairy cows with clinical ketosis.

# Analysis code

The R scripts are organized by main figure and its associated supplementary figures. Each script contains the analysis and visualization code used to generate the corresponding results.

# Differential-expression results

The Excel workbook Supplementary_DE_results.xlsx contains 19 worksheets covering three types of analyses:

FindAllMarkers: Six worksheets containing marker-analysis results for major immune-cell populations and their subclusters, including raw and Bonferroni-adjusted P values.

FindMarkers: Twelve worksheets containing differential-expression results between the healthy and ketosis groups within individual cell types or subsets, including raw and Bonferroni-adjusted P values.

DESeq2: One worksheet containing cow-level pseudobulk differential-expression results for the overall hepatic CD45⁺ immune-cell compartment, including raw and Benjamini–Hochberg-adjusted P values.

Each worksheet begins with a brief description and includes gene identifiers, statistical results, and other analysis-specific information. The workbook retains the supplied analysis outputs, including entries that do not meet the adjusted P-value threshold. Significant genes should be identified using the expression, fold-change, and adjusted P-value criteria described in the manuscript.

# Contact

For questions about the code or analyses, please contact jluzhaocc@126.com.
