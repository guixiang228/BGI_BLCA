#### box plot of gene set score among regions
library(Seurat)
library(dplyr)
library(cowplot)
library(patchwork)
library(ggplot2)
library(future)
options(bitmapType='cairo')
library(igraph)
library(Matrix)
library(MuDataSeurat)
library(reticulate)

sc <- ReadH5AD("Region_tumor_bigcell.h5ad")
str(sc)
seurat_obj <- sc


Basal <- list(c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C"))
Luminal <- list(c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2"))
P53 <- list(c("ACTG2", "CNN1", "MYH11", "MFAP4", "PGM5", "FLNC", "ACTC1", "DES", "PCP4"))
Squamous <- list(c("DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3", "S100A7", "S100A8"))
Neuroendocrine <- list(c("CHGA", "CHGB", "SCG2", "ENO2", "SYP", "NCAM1"))
Cancerstem <- list(c("CD44", "KRT5", "RPSA", "ALDH1A1I"))
EMT <- list(c("ZEB1", "ZEB2", "VIM", "SNAIL", "TWIST1", "FOXC2", "CDH2"))
Claudinlow <- list(c("CLDN3", "CLDN7", "CLDN4", "CDH1", "VIM", "SNAI2", "TWIST1", "ZEB1", "ZEB2"))

type_list <- c(Basal, Luminal, P53, Squamous, Neuroendocrine, Cancerstem, EMT, Claudinlow)
names(type_list) <- c("Basal", "Luminal", "P53", "Squamous", "Neuroendocrine", "Cancerstem", "EMT", "Claudinlow")


seurat_obj <- AddModuleScore(seurat_obj,
                            features = Basal,
                            ctrl = 100,
                            name = paste0("Basal_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Luminal,
                            ctrl = 100,
                            name = paste0("Luminal_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = P53,
                            ctrl = 100,
                            name = paste0("P53_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Squamous,
                            ctrl = 100,
                            name = paste0("Squamous_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Neuroendocrine,
                            ctrl = 100,
                            name = paste0("Neuroendocrine_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Cancerstem,
                            ctrl = 100,
                            name = paste0("Cancerstem_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = EMT,
                            ctrl = 100,
                            name = paste0("EMT_score"))

seurat_obj <- AddModuleScore(seurat_obj,
                            features = Claudinlow,
                            ctrl = 100,
                            name = paste0("Claudinlow_score"))

colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Basal_score1")] <- "Basal"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Luminal_score1")] <- "Luminal"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "P53_score1")] <- "P53"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Squamous_score1")] <- "Squamous"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Neuroendocrine_score1")] <- "Neuroendocrine"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Cancerstem_score1")] <- "Cancerstem"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "EMT_score1")] <- "EMT"
colnames(seurat_obj@meta.data)[which(colnames(seurat_obj@meta.data) == "Claudinlow_score1")] <- "Claudinlow"


library(Seurat)
library(ggplot2)
library(rstatix)
library(ggpubr)


mycol <- c("#48C0AA", "#456990", "#EF767A")

# Basal score
pdf("Basal_score.pdf",height=6,width=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Basal <= 15)
p <- ggplot(data = x, aes(x = merged_cluster, y = Basal)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-10,20) +
  labs(fill = "Region", y='Basal score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank()) 
stat.test <- x %>% t_test(Basal ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Basal)
y_positions <- c(11.5+ 5, 11.5+ 7, 11.5+ 5)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# Luminal score
pdf("Luminal_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Luminal <= 27)
p <- ggplot(data = x, aes(x = merged_cluster, y = Luminal)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +
  ylim(-2,33) +
  labs(fill = "Region", y='Luminal score') +  
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(Luminal ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Luminal)
y_positions <- c(25 + 5, 25 + 7, 25 + 5)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# P53 score
pdf("P53_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$P53 <= 3)
p <- ggplot(data = x, aes(x = merged_cluster, y = P53)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
              width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-2,4) +
  labs(fill = "Region", y='P53 score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(P53 ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$P53)
y_positions <- c(2.3 + 1, 2.3 + 1.35, 2.3 + 1)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# Squamous score
pdf("Squamous_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Squamous <= 7)
p <- ggplot(data = x, aes(x = merged_cluster, y = Squamous)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-5,9) +
  labs(fill = "Region", y='Squamous score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(Squamous ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Squamous)
y_positions <- c(6 + 2, 6 + 2.9, 6 + 2)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# Neuroendocrine score
pdf("Neuroendocrine_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Neuroendocrine <= 0.8)
p <- ggplot(data = x, aes(x = merged_cluster, y = Neuroendocrine)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-0.5,1.1) +
  labs(fill = "Region", y='Neuroendocrine score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(Neuroendocrine ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Neuroendocrine)
y_positions <- c(0.7 + 0.2, 0.7 + 0.3, 0.7 + 0.2)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

## Cancerstem score
pdf("Cancerstem_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Cancerstem <= 25)
p <- ggplot(data = x, aes(x = merged_cluster, y = Cancerstem)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-18,31) +
  labs(fill = "Region", y='Cancer-stem score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(Cancerstem ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Cancerstem)
y_positions <- c(22 + 6, 22 + 9, 22 + 6)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# EMT score
pdf("EMT_score.pdf", width=6,height=6)
p <- ggplot(data = seurat_obj@meta.data, aes(x = merged_cluster, y = EMT)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-1.5,4) +
  labs(fill = "Region", y='EMT score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- seurat_obj@meta.data %>% t_test(EMT ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$EMT)
y_positions <- c(2.5 + 1, 2.5 + 1.3, 2.5 + 1)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

# Claudinlow score
pdf("Claudinlow_score.pdf",width=6,height=6)
x <- filter(seurat_obj@meta.data, seurat_obj@meta.data$Claudinlow <= 10)
p <- ggplot(data = x, aes(x = merged_cluster, y = Claudinlow)) +
  geom_violin(aes(fill = merged_cluster),
              color = 'grey', alpha = 0.8,
              scale = 'width',
              linewidth = 0.3) +  
  geom_boxplot(color = 'lightgrey',  
               width = 0.4,  
               size = 0.4,
               fill = NA,
               outlier.shape = NA) +
  scale_fill_manual(values = mycol) +  
  ylim(-3,12) +
  labs(fill = "Region", y='Claudin-low score') +
  theme_classic() +
  theme(plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "npc"),
        axis.title.x = element_blank())
stat.test <- x %>% t_test(Claudinlow ~ merged_cluster)
str(stat.test)
max_y <- max(seurat_obj@meta.data$Claudinlow)
y_positions <- c(10 + 1, 10 + 2, 10 + 1)
p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", y.position = y_positions)
print(p)
dev.off()

