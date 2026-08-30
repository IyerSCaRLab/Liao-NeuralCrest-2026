# Liao-NeuralCrest-2026
Code repository accompanying the paper by Liao et al. (2026)

## Overview
This repository contains all code used in bulk RNA sequencing analysis in Liao et al. (2026) after processing raw data through the nf-core/rnaseq pipeline. 

## Data Availability
Data associated with this paper is located in GEO.

## Code Usage
1. DataProcessing.R performs QC and generates PCA plots and Spearman's Correlation Heatmaps. 
2. HOXGenePlotting.R generates heatmaps using the VST-transformed expression values of selected HOX genes. 
3. MarkerGenePlotting.R generates heatmaps using the VST-transformed expression values of selected marker genes. 

## Contact
For questions about this code or bioinformatics analysis used in this publication, please contact: 
- Zhixin Liao (Primary Author): zhixin.liao@tufts.edu
- Julian Chow (Code Author): cchow01@tufts.edu
