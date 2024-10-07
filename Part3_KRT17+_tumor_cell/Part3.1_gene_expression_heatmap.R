library(Seurat)
library(pheatmap)
library(ggplot2)
library(MuDataSeurat)
library(ggheatmap)
library(tidyr)
library(Matrix)
library(matrixStats)

seurat_obj <- ReadH5AD('./allsample.h5ad')


genes <- c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT17")
seurat_obj <- NormalizeData(seurat_obj)
data_normalized <- GetAssayData(seurat_obj, assay = "RNA", slot = "data")
batch_info <- seurat_obj@meta.data$batch
merged_cluster_info <- seurat_obj@meta.data$merged_cluster

batch_cluster_means <- data.frame()

# Calculate the mean gene expression for each sample and merged cluster
batch_order <- c("HBCP1", "HBCP3A", "HBCP3B", "HBCP3C", "HBCP5A", "HBCP5B", "HBCP6", "HBCP7", "HBCP8", "HBCP9", "HBCP10", "HBCP11", "HBCP12", "HBCP13", "HBCP14", "HBCP15", "HBCP17", "HBCP18", "All")
cluster_order <- c("11", "2", "10", "0", "3", "14")

for (batch in batch_order) {
  for (cluster in cluster_order) {
    indices <- which(batch_info == batch & merged_cluster_info == cluster)
    if(length(indices) > 0){
      data_subset <- data_normalized[genes, indices, drop = FALSE]
      data_subset_dense <- as.matrix(data_subset)
      if(ncol(data_subset_dense) > 0 && nrow(data_subset_dense) > 0){
        group_means <- rowMeans(data_subset_dense)
        group_means_df <- data.frame(t(group_means))
        colnames(group_means_df) <- genes
        group_means_df$batch <- batch
        group_means_df$merged_cluster <- cluster
        batch_cluster_means <- rbind(batch_cluster_means, group_means_df)
      }
    }
  }
}

# Calculate the mean gene expression of each merged cluster for all samples
for (cluster in cluster_order) {
  indices <- which(merged_cluster_info == cluster)
  if(length(indices) > 0){
    data_subset <- data_normalized[genes, indices, drop = FALSE]
    data_subset_dense <- as.matrix(data_subset)
    if(ncol(data_subset_dense) > 0 && nrow(data_subset_dense) > 0){
      group_means <- rowMeans(data_subset_dense)
      group_means_df <- data.frame(t(group_means))
      colnames(group_means_df) <- genes
      group_means_df$batch <- "All"
      group_means_df$merged_cluster <- cluster
      batch_cluster_means <- rbind(batch_cluster_means, group_means_df)
    }
  }
}


batch_colors <- c("All" = "#FFE4B5", "HBCP1" = "#1f77b4", "HBCP3A" = "#ff7f0e", "HBCP3B" = "#2ca02c", "HBCP3C" = "#d62728",
                  "HBCP5A" = "#9467bd", "HBCP5B" = "#8c564b", "HBCP6" = "#e377c2", "HBCP7" = "#7f7f7f",
                  "HBCP8" = "#bcbd22", "HBCP9" = "#17becf", "HBCP10" = "#ffbb78", "HBCP11" = "#ff9896",
                  "HBCP12" = "#f7b6d2", "HBCP13" = "#c5b0d5", "HBCP14" = "#c49c94", "HBCP15" = "#e377c2",
                  "HBCP17" = "#f5b0ae", "HBCP18" = "#20B2AA")
merged_cluster_colors <- c("0" = "#FF8C00", "2" = "#683ec2", "3" = "#ede737", "10" ="#6bd155", "11" = "#d34a76", "14" = "#00FFFF")
batch_colors <- batch_colors[batch_order]
merged_cluster_colors <- merged_cluster_colors[cluster_order]
annotation_col <- data.frame(
  SD = factor(batch_cluster_means$merged_cluster, levels = cluster_order),
  Sample = factor(batch_cluster_means$batch, levels = batch_order)
)
batch_cluster_means <- batch_cluster_means[order(annotation_col$Sample, annotation_col$SD), ]
annotation_col <- annotation_col[order(annotation_col$Sample, annotation_col$SD), ]
annotation_colors <- list(
  Sample = batch_colors, 
  SD = merged_cluster_colors
)
batch_cluster_means_transposed <- t(batch_cluster_means[, -c(ncol(batch_cluster_means)-1, ncol(batch_cluster_means))])
color_palette <- colorRampPalette(c("#92c5de", "#f7f7f7", "#b2182b"))(100)


pdf("Mean_heatmap_normalized.pdf", width = 20, height = 16)
pheatmap(
  batch_cluster_means_transposed, 
  cluster_rows = FALSE, 
  cluster_cols = FALSE, 
  show_rownames = TRUE, 
  show_colnames = FALSE,
  annotation_col = annotation_col, 
  annotation_colors = annotation_colors, 
  main = "Gene Expression Heatmap (Mean)",
  fontsize_row = 20,
  fontsize_col = 15,
  annotation_legend = TRUE, 
  gaps_col = cumsum(table(annotation_col$Sample)), 
  cellwidth = 10,
  cellheight = 40,
  fontsize = 20, 
  annotation_names_col = TRUE,
  color = color_palette
)
dev.off()