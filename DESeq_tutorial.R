# install.packages("BiocManager")
# BiocManager::install("DESeq2")
# install.packages("dplyr")
# install.packages("readr")
# install.packages("tidyr")
# install.packages("purrr")
# install.packages("pheatmap")

library(DESeq2)
library(dplyr)
library(readr)
library(tidyr)
library(purrr)
library(EnhancedVolcano)
library(pheatmap)

# Define folder containing the count files
alignment_dir <- "/mnt/raid1/despvoulg/Desktop/rna_seq_tutorial/5_Alignment_outputs"

# List all count files
count_files <- list.files(path = alignment_dir, pattern = "*ReadsPerGene.out.tab", full.names = TRUE)

# Initialize an empty list to store the data frames
counts_list <- list()

# Loop though each file and read the data
for (file in count_files) {
  col_count <- 2
  sample_name <- gsub("_ReadsPerGene.out.tab", "", basename(file))

# Skip the 4 first lines and read the data
  temp_df <- read.table(file, skip = 4, header = FALSE, sep = "\t") %>% 
    select(V1, col_count)
  
  colnames(temp_df) <- c("GeneID", sample_name)
  
  counts_list[[sample_name]] <- temp_df
}

# Combine all data frames into one
count_matrix <- purrr::reduce(counts_list, full_join, by = "GeneID")

count_matrix <- count_matrix %>%
                tibble::column_to_rownames(var = "GeneID")

head(count_matrix)

# Organize the data into a DESeq2 compatible format
colData <- data.frame(
  samples = colnames(count_matrix),
  condition = c("LumCtr", "LumPreg", "LumLac", "BasCtr", "BasPreg", "BasLac")
)

# Match the sample names with the count matrix
rownames(colData) <- colData$sample

head(colData)

# same order in ColData and count_matrix
all(rownames(colData) == colnames(count_matrix))

# Create the DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = colData,
                              design = ~ condition)

dds <- DESeq(dds)

# Results
resultsNames(dds)

# Compare conditions

# Effect of pregnancy in Luminal cells
res_LumPreg_vs_LumCtr <- results(dds, contrast = c("condition", "LumPreg", "LumCtr"))

# Effect of lactation in Luminal cells
res_LumLac_vs_LumCtr <- results(dds, contrast = c("condition", "LumLac", "LumCtr"))

# Effect of lactation relative to pregnancy in Luminal cells
res_LumLac_vs_LumPreg <- results(dds, contrast = c("condition", "LumLac", "LumPreg"))

# For the Basal cells
res_BasPreg_vs_BasCtr <- results(dds, contrast = c("condition", "BasPreg", "BasCtr"))
res_BasLac_vs_BasCtr <- results(dds, contrast = c("condition", "BasLac", "BasCtr"))
res_BasLac_vs_BasPreg <- results(dds, contrast = c("condition", "BasLac", "BasPreg"))

results_list <- list(
  "LumPreg_vs_LumCtr" = res_LumPreg_vs_LumCtr,
  "LumLac_vs_LumCtr" = res_LumLac_vs_LumCtr,
  "LumLac_vs_LumPreg" = res_LumLac_vs_LumPreg,
  "BasPreg_vs_BasCtr" = res_BasPreg_vs_BasCtr,
  "BasLac_vs_BasCtr" = res_BasLac_vs_BasCtr,
  "BasLac_vs_BasPreg" = res_BasLac_vs_BasPreg
)

# Loop through the list and create an MA plot for each
for (name in names(results_list)) {
  res <- results_list[[name]]
  png(paste0("MA_plot_", name, ".png")) 
  plotMA(res, main = paste("MA Plot:", name))
  dev.off()
}

for (name in names(results_list)) {
  res <- results_list[[name]]
  if (nrow(res) > 0) {
    png(paste0("Volcano_plot_", name, ".png"), width = 800, height = 800)
    p <- EnhancedVolcano(res,
                         lab = rownames(res),
                         x = 'log2FoldChange',
                         y = 'pvalue',
                         pCutoff = 0.05,
                         FCcutoff = 1.5,
                         title = paste("Volcano Plot:", name))
    
    # Explicitly print the plot to the active graphics device (the PNG file)
    print(p)
    dev.off()
    
    cat(paste0("Created Volcano plot for: ", name, "\n"))
  } else {
    cat(paste0("No data to plot for: ", name, ". Skipping.\n"))
  }
}

# Get significant genes (padj < 0.05, abs(log2FC) > 1) for each comparison
sig_Lum_Preg_vs_Ctr <- rownames(subset(res_LumPreg_vs_LumCtr, padj < 0.05 & abs(log2FoldChange) > 1))
sig_Lum_Lac_vs_Ctr <- rownames(subset(res_LumLac_vs_LumCtr, padj < 0.05 & abs(log2FoldChange) > 1))
sig_Lum_Lac_vs_Preg <- rownames(subset(res_LumLac_vs_LumPreg, padj < 0.05 & abs(log2FoldChange) > 1))

# Combine into a list
gene_lists <- list(
  "LumPreg_vs_LumCtr" = sig_Lum_Preg_vs_Ctr,
  "LumLac_vs_LumCtr" = sig_Lum_Lac_vs_Ctr,
  "LumLac_vs_LumPreg" = sig_Lum_Lac_vs_Preg
)

# Plot using UpSetR
library(UpSetR)
upset(fromList(gene_lists), order.by = "freq")

# Combine all significant genes from all comparisons
all_sig_genes <- unique(c(
  rownames(subset(res_LumPreg_vs_LumCtr, padj < 0.05)),
  rownames(subset(res_LumLac_vs_LumCtr, padj < 0.05)),
  rownames(subset(res_BasPreg_vs_BasCtr, padj < 0.05)),
  rownames(subset(res_BasLac_vs_BasCtr, padj < 0.05))
))

top_genes <- all_sig_genes

# Get the normalized counts for these genes
vsd <- vst(dds, blind = FALSE)
mat <- assay(vsd)[top_genes, ]
mat_scaled <- t(scale(t(mat)))

# Use the pheatmap function
pheatmap::pheatmap(mat_scaled,
                   show_rownames = FALSE,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   main = "Heatmap of Significant Genes across all Conditions")