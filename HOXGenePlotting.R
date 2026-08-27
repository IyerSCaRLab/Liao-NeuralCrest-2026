# This file takes the counts file and processes the data for and plots HOX gene heatmaps. 
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

# set our cells as the reference
dds$cell_publication <- relevel(dds$cell_publication, ref = "This Study")
dds <- estimateSizeFactors(dds)

# No removal of lowly expressed genes to obtain full HOX panel
# VST normalization
vsd <- vst(dds, blind = TRUE)

# Extract HOX Genes
hox_genes <- counts[grep("^HOX", rownames(counts), ignore.case = TRUE), ]
hox_coding <- rownames(hox_genes[grep("^HOX[ABCD][0-9]+$", rownames(hox_genes)), ])

# Extract VST-normalized expression for these genes
hox_interects <- intersect(hox_coding, rownames(assay(vsd)))
mat <- assay(vsd)[hox_interects, ]

# Order HOX genes correctly
HOXcluster <- sub("^HOX([A-D]).*", "\\1", rownames(mat))
HOXnumber <- as.numeric(sub("^HOX[A-D](\\d+)$", "\\1", rownames(mat)))
HOXorder <- order(HOXnumber, HOXcluster)
mat <- mat[HOXorder,]

# Select only representative samples
plot_samples <- c(
  "TNC_SC_D64_1003",
  "TNC_SC_D64_1015",       
  "PrimarySchwanns1",
  "SC_young1",
  "SC_old2"
)

# Subset metadata and matrix to only these samples, and enforce this order
metadata_sub <- metadata[metadata$sample %in% plot_samples, ]
metadata_sub$sample <- factor(metadata_sub$sample, levels = plot_samples)
metadata_sub <- metadata_sub[order(metadata_sub$sample), ]
mat_sub <- mat[, as.character(metadata_sub$sample), drop = FALSE]

# Build annotation data frame with a column name that matches anno_colors
anno_col <- data.frame(
  sample = factor(metadata_sub$sample, levels = plot_samples),
  row.names = colnames(mat_sub),
  check.names = FALSE
)

anno_colors <- list(
  "sample" = c(
    "TNC_SC_D64_1003" = "#D67634",
    "TNC_SC_D64_1015" = "#D67634",
    "PrimarySchwanns1" = "#6C8EA6",
    "SC_young1" = "#A8CEBD",
    "SC_old2" = "#6B9080",
    
  )
)

t_mat = t(mat_sub)
gap_positions <- c(2, 3)

# Draw Heatmap
pheatmap(
  t_mat,
  main = "",
  border_color = "grey30",
  fontsize_col = 8,
  show_colnames = TRUE,
  angle_col = 90,
  show_rownames = FALSE,
  cluster_cols      = FALSE,  
  cluster_rows      = FALSE,
  gaps_row = gap_positions,
  annotation_row = anno_col,
  annotation_colors = anno_colors,
  annotation_names_row = FALSE,
  annotation_legend = FALSE
)

# Forcing bolding of labels
grid.force()
grid.ls()
grid.gedit("GRID.text.186", gp = gpar(fontface = "bold"), grep = TRUE)
grid.gedit("GRID.text.189", gp = gpar(fontface = "bold"), grep = TRUE)
