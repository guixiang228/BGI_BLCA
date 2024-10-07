#### box plot of gene set score among domains
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

paths <- c("HBCP3A", "HBCP3B", "HBCP3C","HBCP11")

for (path in paths) {
  full_path <- paste0("./domain_tumor_DEG/",path,"_bigcell.h5ad")
  sc <- ReadH5AD(full_path)
  str(sc)

 seurat_obj <- sc

  workdir <- ("./box_plot/")
  
  if (!dir.exists(workdir)) {
    dir.create(workdir)
  }
  
setwd(workdir)

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

#######
library(rstatix)
library(ggpubr)
library(gridExtra)


my_order <- c("11", "2", "10", "0", "3", "14")
mycol <- c("#d34a76", "#683ec2", "#6bd155", "#FF8C00", "#ede737", "#00FFFF")
metascape_features <- c("Basal", "Luminal", "P53", "Squamous", "Neuroendocrine", "Cancerstem", "EMT", "Claudinlow")

plots <- list()

for (feature in metascape_features) {
  if (feature %in% colnames(seurat_obj@meta.data)) {
    feature_data <- seurat_obj@meta.data[[feature]]
    min_value <- min(feature_data, na.rm = TRUE)
    max_value <- max(feature_data, na.rm = TRUE)
    x <- seurat_obj@meta.data %>% filter(.data[[feature]] <= max_value)
    x$merged_cluster <- factor(x$merged_cluster, levels = my_order)
    box_summary <- x %>% group_by(merged_cluster) %>% summarise(
      lower = quantile(.data[[feature]], 0.25),
      upper = quantile(.data[[feature]], 0.75)
    )
    x_top_bottom <- x %>%
      left_join(box_summary, by = "merged_cluster") %>%
      filter(.data[[feature]] < lower |.data[[feature]] > upper)
    ylim_range <- c(min_value - 0.1 * (max_value - min_value), max_value + 0.23 * (max_value - min_value))
    
    p <- ggplot(data = x, aes(x = merged_cluster, y =.data[[feature]])) +
      geom_boxplot(aes(fill = merged_cluster),
                   alpha = 0.65, linewidth = 0.5, width = 0.8, outliers = FALSE, staplewidth = 0.5) +  
      geom_jitter(data = x_top_bottom, aes(color = merged_cluster), width = 0.05, size = 0.1, show.legend = FALSE) +
      scale_fill_manual(values = mycol, labels = my_order) +  
      scale_color_manual(values = mycol, labels = my_order) +  
      ylim(ylim_range) +  
      labs(fill = "Region", y = paste(feature, 'score')) +
      theme_classic() +
      theme(plot.margin = unit(c(0, 0, 0, 0), "npc"),
            axis.title.x = element_blank(),
            legend.position = "right",
            legend.justification = "center",
            legend.margin = margin(0, 0, 0, 0))
    formula <- as.formula(paste(feature, "~ merged_cluster"))
    stat.test <- x %>% 
      t_test(formula, paired = FALSE) %>% 
      filter(group1 == "0" | group2 == "0") %>% 
      filter(group1 != group2)
manual_y_positions <- function(stat.test) {
  y_positions <- rep(NA, nrow(stat.test))
  
  for (i in 1:nrow(stat.test)) {
    if ((stat.test$group1[i] == "0" && stat.test$group2[i] == "11") ||
        (stat.test$group1[i] == "11" && stat.test$group2[i] == "0") ||
        (stat.test$group1[i] == "0" && stat.test$group2[i] == "14") ||
        (stat.test$group1[i] == "14" && stat.test$group2[i] == "0")) {
      y_positions[i] <- max_value + 0.03 * (max_value - min_value)  
    } else if ((stat.test$group1[i] == "0" && stat.test$group2[i] == "2") ||
               (stat.test$group1[i] == "2" && stat.test$group2[i] == "0") ||
               (stat.test$group1[i] == "0" && stat.test$group2[i] == "3") ||
               (stat.test$group1[i] == "3" && stat.test$group2[i] == "0")) {
      y_positions[i] <- max_value + 0.13 * (max_value - min_value)  
    } else if ((stat.test$group1[i] == "0" && stat.test$group2[i] == "10") ||
               (stat.test$group1[i] == "10" && stat.test$group2[i] == "0")) {
      y_positions[i] <- max_value + 0.23 * (max_value - min_value)  
    } else {
      y_positions[i] <- max_value + 0.1*(i-1) * (max_value - min_value) + 0.03 * (max_value - min_value) # For additional cases if needed
    }
  }
  
  return(y_positions)
}
    
    y_positions <- manual_y_positions(stat.test)
    
    
    p <- p + stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                                y.position = y_positions, 
                                tip.length = 0.02)
    
    p <- p + scale_x_discrete()
    
    
    plots[[feature]] <- p
  } else {
    warning(paste("Feature column", feature, "does not exist in meta.data"))
  }
}
num_groups <- length(my_order)
plot_width <- 0.5 * num_groups + 2.1 

pdf(paste(path, "_Metascape_scores.pdf"), height = 6, width = plot_width * 2.5)
grid.arrange(grobs = plots, ncol = 4, nrow = 2)
dev.off()

}