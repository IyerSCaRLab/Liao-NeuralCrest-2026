# This code takes the counts file and processes the data for each cell type, 
# performing QC on the data, and generating plots for the PCA and Spearman's Correlation.
# This specific file performs the analysis on Schwann Cells, and was modified to analyze 
# each cell type. 

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

# set Iyer cells as the reference
dds$cell_publication <- relevel(dds$cell_publication, ref = "This Study")

## Filter lowly expressed genes
if (!is.null(sizeFactors(dds))) {
  cat("DESeq2 has already been fitted — skipping filter step.\n")
  cat("Genes in dds:", nrow(dds), "\n")
} else {
  n_before <- nrow(dds)
  keep     <- rowSums(counts(dds)) >= 10
  dds      <- dds[keep, ]
  cat("Genes before filtering:", n_before, "\n")
  cat("Genes after filtering:  ", nrow(dds), "\n")
  cat("Genes removed:          ", n_before - nrow(dds), "\n")
}


# VST normalization
vsd <- vst(dds, blind = TRUE)

# Plot Library sizes
lib_sizes <- data.frame(
  sample    = colnames(dds),
  lib_size  = colSums(counts(dds)),
  cell_publication = dds$cell_publication
)
ggplot(lib_sizes, aes(x = sample, y = lib_size / 1e6, fill = cell_publication)) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  geom_hline(
    yintercept = mean(lib_sizes$lib_size) / 1e6,
    linetype   = "dashed", color = "red", linewidth = 0.8
  ) +
  scale_fill_brewer(
    palette = "Set2",
  ) +
  labs(
    title    = "Library sizes per sample",
    subtitle = "Red dashed line = average library size",
    x        = "Sample",
    y        = "Million mapped reads",
    fill     = "Day"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# VST-normalized count distributions
sample_info <- as.data.frame(colData(dds)) 

# Reshape VST matrix to long format for ggplot
vsd_long <- assay(vsd) %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols      = -gene,
    names_to  = "sample",
    values_to = "vst_count"
  ) %>%
  left_join(sample_info, by = "sample") %>%
  # Order samples by cell source
  mutate(sample = factor(sample, levels = unique(sample[order(cell_publication)])))

ggplot(vsd_long, aes(x = sample, y = vst_count, fill = cell_publication)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.8) +
  geom_boxplot(width = 0.08, fill = "white",
               outlier.shape = NA, coef = 0) +
  scale_fill_brewer(
    palette = "Set2",
  ) +
  labs(
    title    = "VST-normalized count distributions",
    subtitle = "All samples should show similar shape and median",
    x        = "Sample",
    y        = "VST-normalized expression",
    fill     = "Cell Type"
  ) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

# PCA Plot of cell source 
# Extract PCA coordinates from VST-normalized data
pca_data_cellpub <- plotPCA(vsd, intgroup = c("cell_publication"), returnData = TRUE)

# Variance explained by each PC
pct_var_cellpub <- round(100 * attr(pca_data_cellpub, "percentVar"), 1)

# Create factoring for legend order
pca_data_cellpub$cell_publication <- factor(
  pca_data_cellpub$cell_publication,
  levels = c(
    "This Study",
    "Majd et al. LP (2023)",
    "Majd et al. HP (2023)",
    "Primary Human Schwann Cells (Majd 2023)"
  )
)

ggplot(pca_data_cellpub, aes(x = PC1, y = PC2, color = cell_publication, label = sample)) +
  geom_point(size = 5, alpha = 0.85) +
  scale_color_manual(
    values = c(
    "Majd et al. HP (2023)" = "#6B9080",
    "Majd et al. LP (2023)" = "#A8CEBD",
    "This Study"= "#D67634",
    "Primary Human Schwann Cells (Majd 2023)" = "#6C8EA6"), 
    labels = c(
      "Majd et al. HP (2023)" = "Majd et al. HP (2023)",
      "Majd et al. LP (2023)" = "Majd et al. LP (2023)",
      "This Study"= "This Study",
      "Primary Human Schwann Cells (Majd 2023)" = "Primary Human Schwann Cells (Majd 2023)"
    ) 
    ) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 1.2)) +
  labs(
    x        = paste0("PC1: ", pct_var_cellpub[1], "% variance"),
    y        = paste0("PC2: ", pct_var_cellpub[2], "% variance"),
    color    = NULL
  ) +
  theme_bw(base_family = "Arial") + 
  theme(
    panel.grid = element_blank(),                    
    axis.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 20),
    legend.key.height = unit(1.15, "lines"),
    legend.background = element_rect(fill = "white", color = NA),
    legend.text = element_text(size = 15),
    legend.position   = c(0.99, 0.99), # x,y position of legend
    legend.justification = c("right", "top")
  )


# Plotting the Correlation Matrix
vsd_mat <- assay(vsd)
correlation_mat <- cor(vsd_mat, method = "spearman")
cor_dist <- as.dist(1-correlation_mat)

sample_order <- c("PrimarySchwanns1", "PrimarySchwanns2",
                  "TNC_SC_D64_1003", "TNC_SC_D64_1015", "SC_old1", "SC_old2",
                  "SC_young1", "SC_young2")
correlation_mat <- correlation_mat[sample_order, sample_order]

heat_colors <- colorRampPalette(brewer.pal(9, "Blues")[1:7])(255)

annotation_df <- data.frame(
  sample = factor(colnames(correlation_mat), levels = colnames(correlation_mat))
)
rownames(annotation_df) <- colnames(correlation_mat)
sample_colors <- list(
  sample = c(
  "PrimarySchwanns1" = "#6C8EA6", "PrimarySchwanns2" = "#6C8EA6", 
  "TNC_SC_D64_1003" = "#D67634", "TNC_SC_D64_1015" = "#D67634", 
  "SC_old1" = "#6B9080", "SC_old2" = "#6B9080", "SC_young1" = "#A8CEBD",
  "SC_young2" = "#A8CEBD"
  ))

pheatmap(
  correlation_mat,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  annotation_col = annotation_df,
  annotation_row = annotation_df,
  annotation_colors = sample_colors,
  annotation_names_row = FALSE,
  annotation_names_col = FALSE,
  annotation_legend = FALSE,
  show_rownames = FALSE,
  show_colnames = FALSE,
  col = heat_colors,
  fontsize = 12,
  display_numbers = TRUE,
  fontsize_number = 12,
  number_color = "black",
)


# Manually calculating the PCA to find the drivers of genes 
# Get the top variable genes, run PCA directly
rv <- rowVars(assay(vsd))
select <- order(rv, decreasing = TRUE)[seq_len(min(500, length(rv)))]
pca_full <- prcomp(t(assay(vsd)[select, ]), scale. = FALSE)

# Loadings for each gene on each PC1/PC2
loadings <- pca_full$rotation

# Display top genes driving PC1/PC2
n_top <- 15
top_pc1_genes <- names(sort(abs(loadings[, "PC1"]), decreasing = TRUE))[1:n_top]
top_pc2_genes <- names(sort(abs(loadings[, "PC2"]), decreasing = TRUE))[1:n_top]
top_genes <- union(top_pc1_genes, top_pc2_genes)

arrow_df <- data.frame(
  gene = top_genes,
  PC1  = loadings[top_genes, "PC1"],
  PC2  = loadings[top_genes, "PC2"]
)

ggplot(arrow_df) +
  geom_segment(
    aes(x = 0, y = 0, xend = PC1, yend = PC2),
    arrow = arrow(length = unit(0.2, "cm")),
    color = "grey30",
    linewidth = 0.5
  ) +
  geom_text_repel(
    aes(x = PC1, y = PC2, label = gene),
    size = 3.2,
    color = "black",
    max.overlaps = 20
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  scale_x_continuous(expand = expansion(mult = 0.15)) +
  scale_y_continuous(expand = expansion(mult = 0.15)) +
  labs(
    x = paste0("PC1 loading (", pct_var_cellpub[1], "% variance)"),
    y = paste0("PC2 loading (", pct_var_cellpub[2], "% variance)"),
    title = "Top Gene Drivers of PC1 and PC2"
  ) +
  theme_bw(base_family = "Arial") +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(size = 16, face = "bold"),
    axis.title = element_text(size = 14, face = "bold"),
    axis.text  = element_text(size = 12)
  )






