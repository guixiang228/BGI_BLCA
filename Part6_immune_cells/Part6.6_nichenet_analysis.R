### cell communication between immune cell and tumor cell
library(nichenetr)
library(tidyverse)
library(circlize)
library(Seurat)
library(Matrix)
library(Cairo)
library(MuDataSeurat)
options(bitmapType = "cairo")

sc <- ReadH5AD("bigcell.h5ad")
SD11<-subset(x = sc, subset = spatial_domain == "SD11")
SD11$SD <- 'super'
SD2<-subset(x = sc, subset = spatial_domain == "SD2")
SD2$SD <- 'invasive'
SD10<-subset(x = sc, subset = spatial_domain == "SD10")
SD10$SD <- 'invasive'
SD0<-subset(x = sc, subset = spatial_domain == "SD0")
SD0$SD <- 'invasive'

SD=merge(SD11,y=c(SD2,SD10,SD0))

#### select tumor cell
Basal <- list(c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C"))
Luminal <- list(c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2"))
SD <- AddModuleScore(SD, features = Basal, ctrl = 50, nbin=12, name = paste0("Basal_score"))
SD <- AddModuleScore(SD, features = Luminal, ctrl = 50, nbin=12, name = paste0("Luminal_score"))
##1
# remove_cells <- SD@meta.data$merged_cluster %in% c('NMIBC','MIBC') &
#                 SD@meta.data$SD == "invasive" &
#                 SD@meta.data$Basal_score1 < SD@meta.data$Luminal_score1
# cells_to_remove <- rownames(SD@meta.data)[remove_cells]
# keep_cells <- !(colnames(SD) %in% cells_to_remove)
# SD <- SD[,keep_cells]
##1.1
SD1 <- SD %>% subset(merged_cluster %in% c('NMIBC','MIBC') & SD == "invasive")
median_basal <- median(SD1@meta.data$Basal_score1)
remove_cells <- SD@meta.data$merged_cluster %in% c('NMIBC','MIBC') &
                SD@meta.data$SD == "invasive" &
                SD@meta.data$Basal_score1 < median_basal
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]
##2
remove_cells <- SD@meta.data$merged_cluster %in% c('NMIBC','MIBC') &
                SD@meta.data$SD == "super" &
                SD@meta.data$Basal_score1 > SD@meta.data$Luminal_score1
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]
##2.1
SD2 <- SD %>% subset(merged_cluster %in% c('NMIBC','MIBC') & SD == "super")
median_luminal <- median(SD2@meta.data$Luminal_score1)
remove_cells <- SD@meta.data$merged_cluster %in% c('NMIBC','MIBC') &
                SD@meta.data$SD == "super" &
                SD@meta.data$Luminal_score1 < median_luminal
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]

## filter cell types
SD = subset(SD, merged_cluster %in% c('APOE_Macrophages','B_cells','CA10_T_cells','CD1C_cDC2','CD8A_CTL','CLEC9A_cDC1','LAMP3_cDC3','LYVE1_Macrophages','MIBC',
'NCAM1_NK_cells','NLRP3_Macrophages','NMIBC','Plasma_cells','TCF7_T_cells','TNFRSF9_T_cells','TOP2A_T_cells'))
Idents(SD) = SD$merged_cluster


sender_celltypes = c('APOE_Macrophages','B_cells','CA10_T_cells','CD1C_cDC2','CD8A_CTL','CLEC9A_cDC1','LAMP3_cDC3','LYVE1_Macrophages',
'NCAM1_NK_cells','NLRP3_Macrophages','Plasma_cells','TCF7_T_cells','TNFRSF9_T_cells','TOP2A_T_cells')
receiver = c('NMIBC','MIBC')

ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network = readRDS("lr_network_human_21122021.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


rm(seurat_obj)
seurat_obj<-SD
seurat_obj<-subset(SD, merged_cluster %in% c(sender_celltypes,receiver))
Idents(seurat_obj) <- seurat_obj$merged_cluster

list_expressed_genes_receiver = receiver %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.3) # lapply to get the expressed genes of every receiver cell type separately here
expressed_genes_receiver = list_expressed_genes_receiver %>% unlist() %>% unique()

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]


list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, seurat_obj, 0.2) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

DE_table_receiver = FindMarkers(object = seurat_obj, ident.1 = receiver, ident.2 = sender_celltypes, min.pct = 0.10) %>% rownames_to_column("gene")
geneset_oi = DE_table_receiver %>% filter(p_val_adj <= 0.05 & abs(avg_log2FC) >= 0.5) %>% pull(gene)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)

ligand_activities = ligand_activities %>% arrange(-aupr_corrected) %>% mutate(rank = rank(desc(aupr_corrected)))
best_upstream_ligands = ligand_activities %>% top_n(100, aupr_corrected) %>% arrange(-aupr_corrected) %>% pull(test_ligand) %>% unique()

#select ligands with high expression in sender cells
best_upstream_ligands <- c('SPP1','SEMA4D','PLXNA1','S100A4','CD96','NEGR1','LGALS1','COMP','TIMP1','BACE2','CALM1','GZMB','LRPAP1','SIGLEC1','TGFBI','MFGE8','NAMPT',
'SEMA4C','FURIN','CADM1','WNT5B','SEMA6C','ZNRF3','CRLF1','ADGRE5','PLAU','NECTIN2','ST6GAL1','LFNG','PSAP','EPHA4','LRIG2')

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.6)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

### select target genes
interactions <- active_ligand_target_links[order_targets, ]
interaction_strength <- rowSums(interactions, na.rm = TRUE)
sorted_genes <- order(interaction_strength, decreasing = TRUE)
top_genes <- order_targets[sorted_genes[1:30]]

order_targets1 <- c('KRT17','CD44','CDH3','KRT5','KRT14','KRT16','JUN','JUNB','TCF7L2','HES1','FGFR3','KRT18','KRT19','KRT7','GATA3','KRT8','UPK2','FOXA1','ERBB3','PPARG')
order_targets1 <- unique(c(order_targets1,top_genes))
order_targets <- order_targets1[order_targets1 %in% order_targets]
vis_ligand_target <- active_ligand_target_links[order_targets,order_ligands] %>% t()
vis_ligand_target <- vis_ligand_target[,order_targets]
p_ligand_target_network <- vis_ligand_target %>%
  make_heatmap_ggplot("Prioritized ligands", "Predicted target genes", color = "purple", legend_position = "top", x_axis_position = "top", legend_title = "Regulatory potential") +
  theme(
    axis.title.x = element_text(size = 12, face = "italic"), 
    axis.title.y = element_text(size = 12, face = "italic"), 
    axis.text.x = element_text(size = 10, face = "italic"),  
    axis.text.y = element_text(size = 10),                   
    legend.title = element_text(size = 12),                  
    legend.text = element_text(size = 10)                    
  ) +
  scale_fill_gradient2(low = "whitesmoke", high = "purple", breaks = c(0, 0.07, 0.14))
#Diagram of Interaction Strength
pdf("p_ligand_target_network2_interested_ImmunetoTumor.pdf",height=5,width=4.5)
print(p_ligand_target_network)
dev.off()


###ligand expression in sender cells
seurat_obj<-SD
seurat_obj <- subset(seurat_obj, merged_cluster %in% c('APOE_Macrophages','B_cells','CA10_T_cells','CD1C_cDC2','CD8A_CTL','CLEC9A_cDC1','LAMP3_cDC3','LYVE1_Macrophages',
'NCAM1_NK_cells','NLRP3_Macrophages','Plasma_cells','TCF7_T_cells','TNFRSF9_T_cells','TOP2A_T_cells'))
Idents(seurat_obj) <- seurat_obj$merged_cluster
seurat_obj <- NormalizeData(seurat_obj)
ligands <- c('SPP1','SEMA4D','PLXNA1','S100A4','CD96','NEGR1','LGALS1','COMP','TIMP1','BACE2','CALM1','GZMB','LRPAP1','SIGLEC1','TGFBI','MFGE8','NAMPT',
'SEMA4C','FURIN','CADM1','WNT5B','SEMA6C','ZNRF3','CRLF1','ADGRE5','PLAU','NECTIN2','ST6GAL1','LFNG','PSAP','EPHA4','LRIG2')
ligands <- rev(ligands)
pdf("order_targets_interested_ligands_ImmunetoTumor.pdf", height=5.8, width=4)
p <- DotPlot(seurat_obj, features = ligands, cols = "RdYlBu", group.by = "merged_cluster") +
  coord_flip() + 
  scale_y_discrete(position = "right") + 
  RotatedAxis() +
  theme(legend.box = "vertical",
        legend.key.size = unit(0.2, "cm"),
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        axis.title.y = element_blank(),
        axis.title.x = element_blank(), 
        axis.text.x = element_text(size = 8.7,angle = 90,hjust = 0), 
        axis.text.y = element_text(size = 8.5))
print(p)
dev.off()


###target gene expression in receiver cells
seurat_obj<-SD
seurat_obj <- subset(seurat_obj, merged_cluster %in% c('NMIBC','MIBC'))
seurat_obj <- NormalizeData(seurat_obj)
desired_order <- c('SD11','SD2','SD10','SD0')
seurat_obj@meta.data$spatial_domain <- factor(seurat_obj@meta.data$spatial_domain, levels = desired_order)
order_targets <- c('CD44','JUN','JUNB','TCF7L2','HES1','GATA3','FOXA1','ERBB3','PPARG','CDKN1A','VEGFA','CCND1','FOS','TP53','ID1','ID2','BRCA1','ATF3','DUSP1','GDF15','TGIF1','THBS1','CDH1','AR','FN1','PPP1R15A',
'NFKB1','SLC2A1','LPP','JUND','SIRT1','STAT3','TNFSF10','FOSL2','CEBPB','PERP','CKS1B','HSP90AB1','CFLAR','MAP1LC3B','AURKA')
pdf("order_targets_interested_targets_ImmunetoTumor.pdf", height=2, width=8.5)
print(DotPlot(seurat_obj, features = order_targets, cols = "RdYlBu", group.by = "spatial_domain") +
      RotatedAxis() +
      theme(legend.box = "vertical",
            legend.key.size = unit(0.06, "cm"),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title.x = element_text(size = 10)) +
      xlab("Target genes in tumor") 
)
dev.off()



