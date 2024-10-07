### script of monocle2 analysis for tumors in different domains
library(monocle)
library(Seurat)
library(MuDataSeurat)
library(patchwork)
library(gridExtra)
library(RColorBrewer)
library(pheatmap)

sam = 'all_down500_mocnole2_'

obj1 <- ReadH5AD('allsample_umap.h5ad')
Idents(obj1) <- 'merged_cluster'
obj1$merged_cluster <- as.factor(obj1$merged_cluster)
obj <- subset(obj1, downsample = 500)

pd <- new("AnnotatedDataFrame", data = obj@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(obj), row.names = rownames(obj)))
cds <- newCellDataSet(obj@assays$RNA@counts, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr=0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

genes <- FindAllMarkers(obj) %>% filter(p_val_adj < 0.001) 

ordergene <- unique(genes$gene)
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds, root_state = '3')

custom_colors <- c("#FFFF33", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#E41A1C")
p1 <- plot_cell_trajectory(cds, 
                           color_by = "Pseudotime", 
                           cell_size = 0.2, 
                           show_tree = F, 
                           show_branch_points = F) + 
    scale_color_gradientn(colors = c("#09030C","#492960", "#B8504E", "#E3A748", "#FCF7AF"))

p2 <- plot_cell_trajectory(cds, 
                           color_by = "merged_cluster", 
                           cell_size = 0.2, 
                           show_tree = F, 
                           show_branch_points = F, 
                           cell_name_size = 10) + 
    scale_color_manual(values = custom_colors)

p3 <- plot_cell_trajectory(cds, 
                           color_by = "State",
                           cell_size=0.2, 
                           show_tree = F, 
                           show_branch_points = F)
pdf(paste0(sam, '_monocle.pdf'), height = 7, width = 7)
print((p1+p2)/(p3+p3))
dev.off()

p4 <- plot_cell_trajectory(cds, 
                           color_by = "merged_cluster",
                           cell_size = 0.5, 
                           show_tree = F, 
                           show_branch_points = F) + 
    facet_wrap(~merged_cluster, nrow = 2) + scale_color_manual(values = custom_colors)

pdf(paste0(sam, '_monocle_per_SD.pdf'), height = 6, width = 9)
print(p4)
dev.off()

pdf(paste0(sam, '_monocle_tree_plot.pdf'), height = 8, width = 4)
plot_complex_cell_trajectory(cds, 
                             cell_size = 0.6,
                             x = 1, 
                             y = 2, 
                             color_by = "merged_cluster",
                             root_states = '3') + 
    scale_color_manual(values = custom_colors) +
    theme(legend.title = element_blank()) 
dev.off()


pheno_df <- as.data.frame(pData(cds))
table(pheno_df$State, pheno_df$merged_cluster)


experData <- ScaleData(LogNormalize(obj@assays$RNA@counts))
experGenes <- c("UCA1", "FOXA1", "ID1", "GNL3", "NEK7",
                "HMGCS2", "PRSS23", "CRACR2B", "FGF3", "GDF15",
                "TMPRSS2", "TXNRD1", "IFITM2", "GSTA1", "CSTA", 
                "KRT4", "SLPI", "KRT5", "SELENOP", "KRT15",
                "PSMG1", "GATA3", "FASN", "SEPT9", "ITM2C", # SD3
                "TIMP1", "KRT17", "CYP1A1", "DMKN", "KRT16")
picList <- list()
for (i in experGenes){
    cds[[i]] <- experData[i,]
    picList[[i]] <- plot_cell_trajectory(cds, color_by = i, cell_size = 0.2) 
    + scale_color_gradientn(colors = c("#09030C","#492960", "#B8504E", "#E3A748", "#FCF7AF"))
}
pdf(paste0(sam, 'gene_in_track.pdf'), width = 17.5, height = 21)
grid.arrange(grobs = picList, ncol = 5)
dev.off()


listEpi <- list(Basal <- list(c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C","KRT17")),
                Luminal <- list(c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2")),
                Squamous <- list(c("DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3", "S100A7", "S100A8")),
                EMT <- list(c("ZEB1", "ZEB2", "VIM", "SNAIL", "TWIST1", "FOXC2", "CDH2")))
Genesheat <- unlist(listEpi)
# whole pseudotime 
pdf('pseudotime_heatmap_monocle.pdf',width = 3, height = 5)
Pheat <- plot_pseudotime_heatmap(cds[rownames(cds) %in% Genesheat,],
                              #num_clusters = 3, 
                              cores = 4,
                              cluster_rows = F,
                              show_rownames = T,
                              hmcols = colorRampPalette(rev(brewer.pal(9, "RdBu")))(100))
print(Pheat)
dev.off()

listEpi <- list(Basal = c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C","KRT17"),
                Luminal = c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2"),
                Squamous = c("DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3", "S100A7", "S100A8"),
                EMT = c("ZEB1", "ZEB2", "VIM", "SNAIL", "TWIST1", "FOXC2", "CDH2"))
Genesheat <- do.call(rbind, lapply(names(listEpi), function(category) {
    data.frame(Gene = listEpi[[category]], Category = category, stringsAsFactors = FALSE)
}))
Genesheat <- Genesheat[Genesheat$Gene %in% c('DSC1','KRT7','DSG1','KRT1','DSG3','KRT6C',
                                             'KRT16','DSC2','KRT6A','ERBB3','UPK2','XBP1',
                                             'PPARG','KRT19','CD44','CDH3','GATA3','FOXA1',
                                             'KRT14','KRT5','KRT17','KRT18','KRT8'),]
rownames(Genesheat) <- Genesheat$Gene
colnames(Genesheat) <- c('Gene','Set')
Genesheat <- Genesheat[2]
# per branch 
mat <- plot_genes_branched_heatmap(cds = cds[rownames(cds) %in% rownames(Genesheat) ,],
                                   branch_point = 1,                                 
                                   # num_clusters = 1,                               
                                   cores = 4,
                                   cluster_rows = F,
                                   branch_labels = c("Cell fate 1", "Cell fate 2"),                                 
                                   hmcols = colorRampPalette(rev(brewer.pal(9, "RdBu")))(62),                                 
                                   show_rownames = T,
                                   return_heatmap = T
)
pdf('pseudotime_heatmap_monocle_sep_nocluster.pdf', width = 5.4, height = 4)
pheatmap(mat$heatmap_matrix[rownames(Genesheat),], 
         fontsize = 7,
         cellheight = 10,
         cellwidth = 1.2,
         useRaster = TRUE, 
         show_colnames = F,
         cluster_cols = FALSE, 
         cluster_rows = FALSE,
         annotation_row = Genesheat,
         annotation_colors = list(
             `Cell Type` = c("Pre-branch" = "lightgrey", "Cell fate 1" = "#F05662", "Cell fate 2" = '#7990C8')
         ),
         # annotation_row = Pheat$annotation_row, 
         annotation_col = mat$annotation_col,  
         color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(10000),
         gaps_col = mat$col_gap_ind)
dev.off()

save.image(file = "monocel2_downsample_500_environment.RData")







