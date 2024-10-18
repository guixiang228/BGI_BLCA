## Functional gene set score of macrophages
library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(MuDataSeurat)
library(reticulate)
library(Matrix)
library(reshape2)
library(ComplexHeatmap)
library(pheatmap)
library(rstatix)
seurat_obj <- ReadH5AD('macro_bigcell.h5ad')


Inflammation_promoting <- list(c("CCL5", "CD19", "CD8B", "CXCL10", "CXCL13", "CXCL9", "GNLY", "GZMB", "IFNG", "IL12A", "IL12B", "IRF1", "PRF1", "STAT1", "TBX21"))
Parainflammation <- list(c("CXCL10", "PLAT", "CCND1", "LGMN", "PLAUR", "AIM2", "MMP7", "ICAM1", "MX2", "CXCL9", "ANXA1", "TLR2", "PLA2G2D", "ITGA2", "MX1", "HMOX1", "CD276", "TIRAP", "IL33", "PTGES", "TNFRSF12A", "SCARB1", "CD14", "BLNK", "IFIT3", "RETNLB", "IFIT2", "ISG15", "OAS2", "REL", "OAS3", "CD44", "PPARG", "BST2", "OAS1", "NOX1", "PLA2G2A", "IFIT1", "IFITM3", "IL1RN"))
HLA <- list(c("HLA-E", "HLA-DPB2", "HLA-C", "HLA-J", "HLA-DQB1", "HLA-DQB2", "HLA-DQA2", "HLA-DQA1", "HLA-A", "HLA-DMA", "HLA-DOB", "HLA-DRB1", "HLA-H", "HLA-B", "HLA-DRB5", "HLA-DOA", "HLA-DPB1", "HLA-DRA", "HLA-DRB6", "HLA-L", "HLA-F", "HLA-G", "HLA-DMB", "HLA-DPA1"))
Type_I_IFN_Reponse <- list(c("DDX4", "IFIT1", "IFIT2", "IFIT3", "IRF7", "ISG20", "MX1", "MX2", "RSAD2", "TNFSF10"))
Type_II_IFN_Reponse <- list(c("GPR146", "SELP", "AHR"))

type_list <- c(Inflammation_promoting, Parainflammation,  HLA, Type_I_IFN_Reponse,Type_II_IFN_Reponse)
names(type_list) <- c('Inflammation_promoting', 'Parainflammation',  'HLA', 'Type_I_IFN_Reponse','Type_II_IFN_Reponse')


seurat_obj <- CreateSeuratObject(counts = seurat_obj@assays$RNA@counts, meta.data = seurat_obj@meta.data)

seurat_obj <- AddModuleScore(seurat_obj, features = Inflammation_promoting, ctrl = 100, name = "Inflammation_promoting_score")
seurat_obj <- AddModuleScore(seurat_obj, features = Parainflammation, ctrl = 100, name = "Parainflammation_score")
seurat_obj <- AddModuleScore(seurat_obj, features = HLA, ctrl = 100, name = "HLA_score")
seurat_obj <- AddModuleScore(seurat_obj, features = Type_I_IFN_Reponse, ctrl = 100, name = "Type_I_IFN_Reponse_score")
seurat_obj <- AddModuleScore(seurat_obj, features = Type_II_IFN_Reponse, ctrl = 100, name = "Type_II_IFN_Reponse_score")


colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Inflammation_promoting_score1")] <- "Inflammation_promoting"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Parainflammation_score1")] <- "Parainflammation"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "HLA_score1")] <- "HLA"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Type_I_IFN_Reponse_score1")] <- "Type_I_IFN_Reponse"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Type_II_IFN_Reponse_score1")] <- "Type_II_IFN_Reponse"


seurat_obj@meta.data <- seurat_obj@meta.data %>%
  filter(spatial_domain %in% c("SD0", "SD1", "SD2", "SD3", "SD4", "SD5", "SD6", "SD7", "SD8", "SD9", "SD10", "SD11", "SD12", "SD13", "SD14"))
seurat_obj@meta.data$spatial_domain <- factor(seurat_obj@meta.data$spatial_domain, levels = paste0("SD", 0:14))
score_list <- c('Inflammation_promoting', 'Parainflammation',  'HLA', 'Type_I_IFN_Reponse','Type_II_IFN_Reponse')
mycol <- c("#FF8C00", "#d5492f", "#683ec2", "#ede737", "#6d76cb",
           "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
           "#6bd155", "#d34a76", "#c9d754", "#d04dc7", "#00FFFF")

    
for (score in score_list) {
    pdf(paste0('Mac_', '_', score, 'domain_boxplot.pdf'), width = 10, height = 4)

    y_values <- seurat_obj@meta.data[[score]]
    lower_limit <- quantile(y_values, 0.05, na.rm = TRUE)
    upper_limit <- quantile(y_values, 0.95, na.rm = TRUE)

    p <- ggplot(seurat_obj@meta.data, aes(x = spatial_domain, y = !!rlang::sym(score), fill = spatial_domain)) +
        geom_boxplot(color = 'black', width = 0.6, size = 0.4, outlier.shape = NA) +
        geom_jitter(aes(color = spatial_domain), width = 0.05, size = 0.1, show.legend = FALSE) +
        scale_fill_manual(values = mycol) +
        scale_color_manual(values = mycol) +
        scale_y_continuous(limits = c(lower_limit, upper_limit)) +
        theme_classic() +
        theme(
            axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
            plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "npc"),
            axis.title.x = element_blank()
        )

    print(p)
    dev.off()
}


features <- c('Inflammation_promoting', 'Parainflammation',  'HLA', 'Type_I_IFN_Reponse','Type_II_IFN_Reponse')
all_p_values <- list()
for (feature in features) {
  if (feature %in% colnames(seurat_obj@meta.data)) {
    x <- seurat_obj@meta.data
    x$spatial_domain <- as.factor(x$spatial_domain)

    # t test
    formula <- as.formula(paste(feature, "~", "spatial_domain"))
    stat.test <- t_test(x, formula)

    group1 <- stat.test$group1
    group2 <- stat.test$group2
    p_values <- stat.test$p.adj
    significance <- ifelse(p_values < 0.001, "****",
                           ifelse(p_values < 0.01, "***",
                                  ifelse(p_values < 0.05, "**",
                                         ifelse(p_values < 0.1, "*", "ns"))))
    p_value_df <- data.frame(Feature = feature, Group1 = group1, Group2 = group2, Significance = significance)
    all_p_values[[feature]] <- p_value_df
  }
}

combined_p_values <- do.call(rbind, all_p_values)
write.csv(combined_p_values, "p_values_significance.csv", row.names = FALSE)
