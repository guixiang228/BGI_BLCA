### cell communication analysis between caf and tumor cells
library(nichenetr)
library(tidyverse)
library(circlize)
library(Seurat)
library(Matrix)
library(Cairo)
library(MuDataSeurat)
options(bitmapType = "cairo")

load('SD.RData')

#### select tumor cell
Basal <- list(c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C"))
Luminal <- list(c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2"))
SD <- AddModuleScore(SD, features = Basal, ctrl = 50, nbin=12, name = paste0("Basal_score"))
SD <- AddModuleScore(SD, features = Luminal, ctrl = 50, nbin=12, name = paste0("Luminal_score"))
##1
# remove_cells <- SD@meta.data$merged_cluster == "Metastatic_tumor" &
#                 SD@meta.data$spatial_domain == "SD14" &
#                 SD@meta.data$Basal_score1 < SD@meta.data$Luminal_score1
# cells_to_remove <- rownames(SD@meta.data)[remove_cells]
# keep_cells <- !(colnames(SD) %in% cells_to_remove)
# SD <- SD[,keep_cells]
##1.1
SD1 <- SD %>% subset(merged_cluster == "Metastatic_tumor" & spatial_domain == "SD14")
median_basal <- median(SD1@meta.data$Basal_score1)
remove_cells <- SD@meta.data$merged_cluster == "Metastatic_tumor" &
                SD@meta.data$spatial_domain == "SD14" &
                SD@meta.data$Basal_score1 < median_basal
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]
##2
remove_cells <- SD@meta.data$merged_cluster == "Metastatic_tumor" &
                SD@meta.data$spatial_domain == "SD3" &
                SD@meta.data$Basal_score1 > SD@meta.data$Luminal_score1
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]
##2.1
SD2 <- SD %>% subset(merged_cluster == "Metastatic_tumor" & spatial_domain == "SD3")
median_luminal <- median(SD2@meta.data$Luminal_score1)
remove_cells <- SD@meta.data$merged_cluster == "Metastatic_tumor" &
                SD@meta.data$spatial_domain == "SD3" &
                SD@meta.data$Luminal_score1 < median_luminal
cells_to_remove <- rownames(SD@meta.data)[remove_cells]
keep_cells <- !(colnames(SD) %in% cells_to_remove)
SD <- SD[,keep_cells]


sender_celltypes = c('ABL2+CAF','FBLN1+CAF',"SPINT2+CAF","ABCC4+myCAF","MYOCD+myCAF",'FBXO32+myCAF','IL6+iCAF','POSTN+myCAF','SLC14A1+myCAF','ZBTB16+CAF')
receiver = c('Metastatic_tumor')

ligand_target_matrix <- readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network = readRDS("lr_network_human_21122021.rds")
weighted_networks <- readRDS("weighted_networks_nsga2r_final.rds")
weighted_networks_lr = weighted_networks$lr_sig %>% inner_join(lr_network %>% distinct(from,to), by = c("from","to"))


seurat_obj<-SD
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
best_upstream_ligands <- c('CADM1','CMTM8','PDGFA','DKK3','IGF1','CX3CL1','WNT7B','POSTN','PLAU','LAMA2','FURIN','ANGPTL4','HMGB1','CALR','ANGPTL2','INHBA','BMP3','TGFA','PSEN1','DCN','WNT5A','MDK','GPI','CLDN1','CALM2','EFNA5','CCN2','IGF2','CALM1')

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 50) %>% bind_rows() %>% drop_na()
active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.6)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names()
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names()

### select target genes
interactions <- active_ligand_target_links[order_targets, ]
interaction_strength <- rowSums(interactions, na.rm = TRUE)
sorted_genes <- order(interaction_strength, decreasing = TRUE)
top_genes <- order_targets[sorted_genes[1:20]]

order_targets1 <- c('KRT17','CD44','CDH3','KRT5','KRT16','JUN','GATA3','KRT8','UPK2','FOXA1','ERBB3','PPARG','JUNB','FOS')
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
  scale_fill_gradient2(low = "whitesmoke", high = "purple", breaks = c(0, 0.05, 0.1))
#Diagram of Interaction Strength
pdf("p_ligand_target_network2_interested_CAFtoTumor.pdf",height=6,width=6)
print(p_ligand_target_network)
dev.off()


###ligand expression in sender cells
ligands <- c('CADM1','CMTM8','PDGFA','DKK3','IGF1','CX3CL1','WNT7B','POSTN','PLAU','LAMA2','FURIN','ANGPTL4','HMGB1','CALR','ANGPTL2','INHBA','BMP3','TGFA','PSEN1','DCN','WNT5A','MDK','GPI','CLDN1','CALM2','EFNA5','CCN2','IGF2','CALM1')
ligands <- rev(ligands)
seurat_obj<-SD
seurat_obj <- subset(seurat_obj, merged_cluster %in% c('ABL2+CAF','FBLN1+CAF','SPINT2+CAF','MYOCD+myCAF','FBXO32+myCAF','IL6+iCAF','POSTN+myCAF','SLC14A1+myCAF','ABCC4+myCAF'))
seurat_obj <- NormalizeData(seurat_obj)
pdf("order_targets_interested_ligands_CAFtoTumor.pdf", height=5.5, width=4)
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
        axis.text.x = element_text(size = 8.5,angle = 90,hjust = 0), 
        axis.text.y = element_text(size = 8.5))
print(p)
dev.off()

###target gene expression in receiver cells
order_targets <- c('KRT17','CD44','CDH3','KRT5','KRT16','JUN','GATA3','KRT8','UPK2','FOXA1','ERBB3','PPARG','JUNB','FOS','GDF15','SERPINE1','STAT3','CCND1','CDKN1A','CYP1A1','DDIT4','MET','BHLHE40','CDH1','TXNIP','HES1')
seurat_obj<-SD
seurat_obj <- subset(seurat_obj, merged_cluster %in% c('Metastatic_tumor'))
seurat_obj <- NormalizeData(seurat_obj)
pdf("order_targets_interested_targets_CAFtoTumor.pdf", height=2, width=5.5)
print(DotPlot(seurat_obj, features = order_targets, cols = "RdYlBu", group.by = "spatial_domain") +
      RotatedAxis() +
      theme(legend.box = "vertical",
            legend.key.size = unit(0.06, "cm"),
            legend.title = element_text(size = 5),
            legend.text = element_text(size = 5),
            axis.title.y = element_blank(),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.x = element_text(size=10)) + 
      xlab("Target genes in tumor") 
)
dev.off()


#### ligand-receptor-target network
library(DiagrammeR)
library(DiagrammeRsvg)

weighted_networks = readRDS("weighted_networks_nsga2r_final.rds")
ligand_tf_matrix = readRDS("ligand_target_matrix_nsga2r_final.rds")
lr_network = readRDS("lr_network_human_allInfo_30112033.rds")
sig_network = readRDS("signaling_network_human_21122021.rds")
gr_network = readRDS("gr_network_human_21122021.rds")


ligands_all = c("POSTN")
targets_all = c("FOS")
# extract network
active_signaling_network = get_ligand_signaling_path(
  ligand_tf_matrix = ligand_tf_matrix, 
  ligands_all = ligands_all, 
  targets_all = targets_all,
  weighted_networks = weighted_networks)
lapply(active_signaling_network, head)
## normalize weights
active_signaling_network_min_max = active_signaling_network
active_signaling_network_min_max$sig = active_signaling_network_min_max$sig %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)
active_signaling_network_min_max$gr = active_signaling_network_min_max$gr %>% 
  mutate(weight = ((weight-min(weight))/(max(weight)-min(weight))) + 0.75)

graph_min_max = diagrammer_format_signaling_graph(
  signaling_graph_list = active_signaling_network_min_max, 
  ligands_all = ligands_all, 
  targets_all = targets_all, 
  sig_color = "indianred", 
  gr_color = "steelblue")
DiagrammeR::render_graph(graph_min_max, 
                         output="graph", # graph, visNetwork
                         layout = "kk" # nicely, circle, tree, kk, and fr
)
