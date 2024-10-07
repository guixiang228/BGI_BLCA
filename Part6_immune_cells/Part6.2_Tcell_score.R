## Functional gene set score of T cell
library(tidyverse)
library(Seurat)
library(ggplot2)
library(dplyr)
library(MuDataSeurat)
library(reticulate)
library(Matrix)
library(reshape2)
library(ComplexHeatmap)
library(cowplot)
library(rstatix)
seurat_obj <- ReadH5AD('Tcell_bigcell.h5ad')

CoStimulatory <- list(c("CD27", "CD28", "CD40LG", "ICOS", "SLAMF1", "TNFRSF14", "TNFRSF18", "TNFRSF4", "TNFRSF9", "ICOSLG", "TNFRSF8", "CD226", "TNFRSF25"))
Cytotoxic <- list(c("GNLY", "GZMA", "GZMB", "GZMK", "IFNG", "NKG7", "PRF1", "CST7", "TNFSF10", "CCL4", "CCL3", "FASLG", "CD44"))
Exhaustion <- list(c("BTLA", "CD276", "CTLA4", "ENTPD1", "HAVCR2", "IDO1", "KLRC1", "LAG3", "LAYN", "LGALS9", "LILRB2", "LILRB4", "CD274", "PDCD1LG2", "TDO2", "TIGIT", "VSIR", "PDCD1", "CXCL13"))
Inhibitory <- list(c("BTLA", "VSIR", "CD160", "CD244", "CD274", "CTLA4", "HAVCR2", "LAG3", "LAIR1", "TIGIT"))

type_list <- c(CoStimulatory, Cytotoxic, Exhaustion, Inhibitory)
names(type_list) <- c('CoStimulatory', 'Cytotoxic', 'Exhaustion', 'Inhibitory')


seurat_obj <- CreateSeuratObject(counts = seurat_obj@assays$RNA@counts, meta.data = seurat_obj@meta.data)
seurat_obj <- AddModuleScore(seurat_obj,
                            features = CoStimulatory,
                            ctrl = 100,
                            name = paste0("CoStimulatory_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Cytotoxic,
                            ctrl = 100,
                            name = paste0("Cytotoxic_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Exhaustion,
                            ctrl = 100,
                            name = paste0("Exhaustion_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Inhibitory,
                            ctrl = 100,
                            name = paste0("Inhibitory_score"))

colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "CoStimulatory_score1")] <- "CoStimulatory"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Cytotoxic_score1")] <- "Cytotoxic"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Exhaustion_score1")] <- "Exhaustion"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Inhibitory_score1")] <- "Inhibitory"


seurat_obj@meta.data <- seurat_obj@meta.data %>%
  filter(spatial_domain %in% c("SD0", "SD1", "SD2", "SD3", "SD4", "SD5", "SD6", "SD7", "SD8", "SD9", "SD10", "SD11", "SD12", "SD13", "SD14"))
seurat_obj@meta.data$spatial_domain <- as.factor(seurat_obj@meta.data$spatial_domain)
score_list <- c("CoStimulatory", "Cytotoxic", "Exhaustion", "Inhibitory")
seurat_obj@meta.data$spatial_domain <- factor(seurat_obj@meta.data$spatial_domain, levels = paste0("SD", 0:14))
mycol <- c("#FF8C00", "#d5492f", "#683ec2", "#ede737", "#6d76cb",
           "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
           "#6bd155", "#d34a76", "#c9d754", "#d04dc7", "#00FFFF")
  
for (score in score_list) {
    pdf(paste0('T_', '_', score, 'sixdomain_boxplot.pdf'), width = 10, height = 4)
    
    y_values <- seurat_obj@meta.data[[score]]
    lower_limit <- quantile(y_values, 0.05, na.rm = TRUE)
    upper_limit <- quantile(y_values, 0.95, na.rm = TRUE)

    p <- ggplot(seurat_obj@meta.data, aes(x = spatial_domain, y = !!rlang::sym(score), fill = spatial_domain)) +
        geom_boxplot(width = 0.6, size = 0.4, outlier.shape = NA) +
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