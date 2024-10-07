# calculate proportion of immune cells in domains of all samples
library(tidyverse)
library(Seurat)
library(MuDataSeurat)
library(SeuratDisk)
library(ggplot2)
library(cowplot)
library(dplyr)
library(rstatix)

workdir <- "./domain_immune_stat/"
setwd(workdir)
seurat_obj <- ReadH5AD('ALL_immune_subtype.h5ad')

seurat_obj@meta.data$spatial_domain <- as.factor(seurat_obj@meta.data$spatial_domain)
selected_subtypes <- c("Tcell") 
# selected_subtypes <- c("Macrophage") 
# selected_subtypes <- c("Bcell","Plasmocyte")

# Calculate the percentage of selected subtypes in each spatial domain and sample
percentage_by_domain_batch <- seurat_obj@meta.data %>%
  filter(annotated_cluster %in% selected_subtypes) %>%
  group_by(spatial_domain, batch) %>%
  summarize(percent = n() / nrow(seurat_obj@meta.data[seurat_obj@meta.data$spatial_domain == unique(spatial_domain) & seurat_obj@meta.data$batch == unique(batch), ]) * 100, .groups = 'drop')

seurat_obj@meta.data$spatial_domain <- factor(seurat_obj@meta.data$spatial_domain, levels = paste0("SD", 0:14))

mycol <- c("#FF8C00", "#d5492f", "#683ec2", "#ede737", "#6d76cb",
           "#ce9d3f", "#81357a", "#d3c3a4", "#3c2f5a", "#b96f49",
           "#6bd155", "#d34a76", "#c9d754", "#d04dc7", "#00FFFF")

seurat_obj@meta.data$batch <- factor(seurat_obj@meta.data$batch, levels = c('HBCP1','HBCP3A', 'HBCP3B', 'HBCP3C', 'HBCP5A', 
                                                                            'HBCP5B', 'HBCP6', 'HBCP7', 'HBCP8', 'HBCP9', 
                                                                            'HBCP10', 'HBCP11', 'HBCP12', 'HBCP13', 'HBCP14', 
                                                                            'HBCP15', 'HBCP17', 'HBCP18'))  
batch_color <- c("#ed1299", "#246b93", "#cc8e12", "#d561dd", "#4aef7b", 
                   "#e86502", "#9ed84e", "#39ba30", "#6ad157", "#8249aa", 
                   "#99db27", "#e07233", "#ff523f", "#ce2523", "#f7aa5d", 
                   "#cebb10", "#931635", "#373bbf")


pdf('Tcell_domain_batch_violin.pdf', width = 10, height = 4)
p <- ggplot(percentage_by_domain_batch, aes(x = spatial_domain, y = percent, fill = spatial_domain)) +
    geom_violin(alpha = 0.8, scale = 'width', linewidth = 0.3) +  
    geom_jitter(aes(color = batch), width = 0.2, alpha = 0.5, size = 1) +
    geom_boxplot(color = 'lightgrey', width = 0.4, size = 0.4, fill = NA, outlier.shape = NA) +
    labs(x = "Spatial Domain", y = "Proportion of Annotated Cluster") +
    scale_fill_manual(values = mycol) +
    scale_color_manual(values = batch_color) +
    theme_classic() +
    theme(
        axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "npc"),
        axis.title.x = element_blank()
    )
print(p)
dev.off()