# This file takes the counts file and processes the data for and plots marker gene heatmaps. 
# This specific file performs the analysis on Schwann Cells, and was modified to analyze each cell type. 

# Load libraries
library(rlang)
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(EnhancedVolcano)
library(clusterProfiler)
library(org.Hs.eg.db)
library(dplyr)
library(tidyverse)
library(tidyr)
library(tibble)
library(RColorBrewer)
library(ggrepel)
library(grid)

# File management
counts_file <- "../Combined_scaled_counts_v3.tsv"
metadata_file <- "../sample_metadata_v3.tsv"
out_dir <- "../Scaled_RNAAnalysisOutput/"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
counts_raw <- read.table(
  counts_file,
  header = TRUE,
  sep    = "\t"
)
metadata <- read.table(
  metadata_file,
  header = TRUE,
  sep ="\t"
)

# Preprocessing of counts matrix to prepare for DEseq2
# remove duplicate genes
rownames(counts_raw) <- make.unique(counts_raw$gene_name)

# Drop annotation columns; keep only the numeric count columns
counts_raw$gene_id   <- NULL
counts_raw$gene_name <- NULL

# Round fractional Salmon counts to integers (required by DESeq2)
counts <- round(counts_raw)

# Extract only the Schwann Cells 
schwann_samples <- metadata$cell_type == "Schwanns"
counts <- counts[, schwann_samples]
metadata <- metadata[schwann_samples, , drop = FALSE]

# Build DeSeqDataSet 
dds <- DESeqDataSetFromMatrix(
  countData = counts,
  colData   = metadata,
  design    = ~ cell_publication
)

# set this study as the reference
dds$cell_source <- relevel(dds$cell_publication, ref = "This Study")
dds <- estimateSizeFactors(dds)

# No removal of lowly expressed genes to obtain full gene panel
# VST normalization
vsd <- vst(dds, blind = TRUE)

# Extract Marker Genes
marker_genes <- read.csv("../MarkerGenes_v1.csv",header = TRUE)

neural_crest_genes <- marker_genes$Neural.Crest
schwann_genes <- marker_genes$Schwann

# Clean genes
neural_crest_genes <- neural_crest_genes[
  !is.na(neural_crest_genes) &
    trimws(neural_crest_genes) != "" &
    !grepl("\\s", neural_crest_genes)
]

schwann_genes <- schwann_genes[
  !is.na(schwann_genes) &
    trimws(schwann_genes) != "" &
    !grepl("\\s", schwann_genes)
]

selected_genes <- c(neural_crest_genes, schwann_genes)

# Extract VST-normalized expression for these genes
marker_intersects <- selected_genes[selected_genes %in% rownames(assay(vsd))]
mat <- assay(vsd)[marker_intersects, ]

# Gap positioning
neural_crest_count <- sum(marker_intersects %in% neural_crest_genes)
gap_position <- neural_crest_count

# Order by cell_source
metadata$cell_source <- factor(
  metadata$cell_source,
  levels = c("Other", "Iyer", "Primary")
)

col_order <- order(metadata$cell_source)
mat <- mat[, col_order]

# Z-score each gene (row) so we see relative patterns across samples
mat_scaled <- t(scale(t(mat)))

# Clip extreme values to prevent outliers from dominating the color scale
mat_scaled[mat_scaled >  2.5] <-  2.5
mat_scaled[mat_scaled < -2.5] <- -2.5

# Ordering mat_scaled and adding colours 
sample_order <- c("PrimarySchwanns1", "PrimarySchwanns2",
                  "TNC_SC_D64_1003", "TNC_SC_D64_1015", "SC_old1", "SC_old2",
                  "SC_young1", "SC_young2")
mat_scaled <- mat_scaled[, sample_order]

annotation_df <- data.frame(
  sample = factor(colnames(mat_scaled), levels = colnames(mat_scaled))
)
rownames(annotation_df) <- colnames(mat_scaled)

sample_colors <- list(
  sample = c(
  "PrimarySchwanns1" = "#6C8EA6", "PrimarySchwanns2" = "#6C8EA6", 
  "TNC_SC_D64_1003" = "#D67634", "TNC_SC_D64_1015" = "#D67634", 
  "SC_old1" = "#6B9080", "SC_old2" = "#6B9080", "SC_young1" = "#A8CEBD",
  "SC_young2" = "#A8CEBD"
  ))

# Draw Heatmap
pheatmap(
  mat_scaled, # mat for VSD and mat_scaled for z-scored
  main = "",
  border_color = "grey30",
  gaps_row = gap_position,
  fontsize = 24,
  show_colnames = FALSE,
  show_rownames = TRUE,
  cluster_cols      = FALSE,  
  cluster_rows      = FALSE,
  annotation_col = annotation_df,
  annotation_colors = sample_colors,
  annotation_names_col = FALSE,
  annotation_legend = FALSE
)


