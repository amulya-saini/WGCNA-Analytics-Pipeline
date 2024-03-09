Weighted Gene Co-expression Network Analysis (WGCNA) Pipeline
Overview:
This GitHub repository contains a comprehensive WGCNA pipeline for analyzing gene expression data and extracting meaningful biological insights. The pipeline covers data loading, pre-processing, soft threshold selection, hierarchical clustering, module identification, and Module-Trait Relationship analysis.

Key Features:
Data Loading and Pre-processing:

Load gene expression and traits data from CSV files.
Handle missing values, set appropriate row/column names, and perform necessary data transformations.
Soft Threshold Selection:

Choose an optimal soft-thresholding power using network topology analysis.
Visualize the scale-free topology model fit and mean connectivity for different soft-thresholding powers.
Gene Clustering:

Perform hierarchical clustering on the dissimilarity matrix derived from the Topological Overlap Matrix (TOM).
Identify gene modules using dynamic tree cut approach and visualize them on a dendrogram.
Module Eigengenes:

Calculate module eigengenes representing overall gene expression patterns within each module.
Cluster module eigengenes and visualize the dynamic relationships between them.
Merging Modules:

Automatically merge gene modules based on a dissimilarity threshold, providing a hierarchical structure of gene expression data.
Module-Trait Relationship Analysis:

Assess the correlation between module eigengenes and external traits.
Generate a heatmap illustrating the relationships between gene modules and traits.
