## script for distinguishing three tumor states and following analysis
library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggpubr)
set.seed(42)


obj <- readRDS('ALL_tumor_score.RDS')
obj <- subset(obj, nFeature_RNA >1000)
expr_matrix <- obj@assays$RNA@counts

Basal_gene <- c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT17")
Luminal_gene <- c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2", "EPCAM")
obj$Basal_gene_count <- colSums(expr_matrix[Basal_gene, ])
obj$Luminal_gene_count <- colSums(expr_matrix[Luminal_gene, ])
obj <- subset(obj, Basal_gene_count > 0 & Luminal_gene_count > 0)

obj <- CreateSeuratObject(counts = obj@assays$RNA@counts, meta.data = obj@meta.data)
obj <- NormalizeData(obj)


Basal_gene <- c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C", "KRT17")
Luminal_gene <- c("CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2", "EPCAM")
gene_union <- union(Basal_gene, Luminal_gene)
for (i in gene_union){
    obj[[i]] <-  obj@assays$RNA@data[i,]
}

type_list <- list(
    Basal =list(c("CD44", "CDH3", "KRT1", "KRT14", "KRT16", "KRT5", "KRT6A", "KRT6B", "KRT6C","KRT17")),
    Luminal = list(c("EPCAM","CYP2J2", "ERBB2", "ERBB3", "FGFR3", "FOXA1", "GATA3", "GPX2", "KRT18", "KRT19", "KRT20", "KRT7", "KRT8", "PPARG", "XBP1", "UPK1A", "UPK2")),
    P53 = list(c("ACTG2", "CNN1", "MYH11", "MFAP4", "PGM5", "FLNC", "ACTC1", "DES", "PCP4")),
    Squamous = list(c("DSC1", "DSC2", "DSC3", "DSG1", "DSG2", "DSG3", "S100A7", "S100A8")),
    Neuroendocrine = list(c("CHGA", "CHGB", "SCG2", "ENO2", "SYP", "NCAM1")),
    Cancerstem = list(c("CD44", "KRT5", "RPSA", "ALDH1A1I")),
    EMT = list(c("ZEB1", "ZEB2", "VIM", "SNAIL", "TWIST1", "FOXC2", "CDH2")),
    Claudinlow = list(c("CLDN3", "CLDN7", "CLDN4", "CDH1", "VIM", "SNAI2", "TWIST1", "ZEB1", "ZEB2")))

for (type in names(type_list)) {
    score_name <- paste0(type, "_score")
    obj <- AddModuleScore(obj, features = type_list[[type]], ctrl = 100, name = score_name)
    colnames(obj@meta.data) <- sub("score1$", "score", colnames(obj@meta.data))
}

meta <- obj@meta.data[,50:ncol(obj@meta.data)]
meta$max <-  pmax(meta$Basal_score, meta$Luminal_score)
meta$subs <- meta$Basal_score - meta$Luminal_score
meta$distance <- (meta$subs)^2 + (meta$max)^2
meta$sum <- meta$Basal_score + meta$Luminal_score
meta$cellID <- rownames(meta)

pdf('EPCAM_KRT17_new.pdf')
for (i in colnames(meta)[1:35]){
    meta <- meta[order(meta[[i]]), ]
    print(ggscatter(meta, x = 'subs', y = 'max', color = i,size = 1)+ scale_color_gradientn(colours = c("#92c5de", "#f7f7f7", "#b2182b")))
}
dev.off()


basal_id <- meta %>%
    filter(subs > 0.3) %>%
    filter(Basal_score  > 0.6) %>%
    filter(Luminal_score  < 0.3) %>%
    arrange(desc(Basal_score)) %>%
    slice(1:1000) %>%
    #arrange(desc(distance)) %>%
    #slice(1:1000) %>%
    pull(cellID)

luminal_id <- meta %>% 
    filter(subs < -0.5)  %>%
    filter(Luminal_score > 0.6) %>%
    filter(Basal_score  < 0.3) %>%
    arrange(desc(Luminal_score)) %>%
    slice(1:2000) %>%
    arrange(desc(distance)) %>%
    slice(1:1000) %>%
    pull(cellID)

Intermediate_ID <- meta %>% 
    #filter(abs(subs) < 0.3)  %>%
    filter(Basal_score > 0.5)  %>%
    filter(Luminal_score > 0.5) %>%
    #filter(sum > 0.3)  %>%
    arrange(desc(max)) %>%
    slice(1:1500) %>%
    pull(cellID)

meta <- meta %>% arrange(factor(Condition, levels = c('None', 'Basal', 'Luminal', 'Intermediate')))
    
meta$Condition <- ifelse(rownames(meta) %in% basal_id, "Basal",  
                          ifelse(rownames(meta) %in% luminal_id, "Luminal",
                                 ifelse(rownames(meta) %in% setdiff(Intermediate_ID, basal_id)[1:1000], 'Intermediate', 
                                        ifelse(meta$max < 0,'None','None'))))
pdf('EPCAM_KRT17_group.pdf')
ggscatter(meta, x = 'subs', y = 'max', color = 'Condition',size = 0.5) + scale_color_manual(values = c('Basal' = '#EA8379', 'Luminal' = '#B395BD', 'None' = 'lightgrey', 'Intermediate' = '#7DAEE0'))
dev.off()
table(meta$Condition)

### 3000 rds processing
meta <- meta[rownames(obj@meta.data),]
obj$Condition <- meta$Condition
Idents(obj) <- 'Condition'
obj1 <- subset(obj, idents = 'None', invert = T)
dim(obj1);saveRDS(obj1,'3000_new.rds')

obj1@meta.data <- obj1@meta.data[,c(1:49,85)]
obj1 <- ScaleData(obj1)
obj1@assays$RNA@var.features <- gene_union
obj1 <- RunPCA(obj1)
obj1 <- RunUMAP(obj1, dims = 1:5)
obj1 <- RunTSNE(obj1, dims = 1:5)

pdf('3000_new_umap.pdf')
DimPlot(obj1, group.by = 'Condition', reduction = 'umap')
dev.off()

pdf('3000_new_umap_gene.pdf', width = 20, height = 20)
FeaturePlot(obj1, features = gene_union, ncol = 5, reduction = 'umap', slot = 'data', order = T) * scale_color_gradientn(colours = c("#92c5de", "#f7f7f7", "#b2182b"))
dev.off()


markers <- FindAllMarkers(obj1, only.pos = T, logfc.threshold = 0.5)
tf <- read.table('hs_hgnc_curated_tfs.txt')
markers$TF <- rownames(markers) %in% tf$V1
openxlsx::write.xlsx(markers, '3000_cells_markers.xlsx')

#### cytotrace 
expression_data <- obj1@assays$RNA@counts
annotation <- obj1[['Condition']]
library(CytoTRACE2) #loading
cytotrace2_result <- cytotrace2(expression_data, species = 'human', ncores = 1)
plots <- plotData(cytotrace2_result = cytotrace2_result, 
                  annotation = annotation,
                  expression_data = expression_data
)
obj1$cytotrace <- cytotrace2_result$CytoTRACE2_Score
plots$CytoTRACE2_UMAP
plots$Phenotype_UMAP
pdf('cytotrace2.pdf', width = 4, height = 4)
plots$CytoTRACE2_Boxplot_byPheno
FeaturePlot(obj1, features = 'cytotrace', reduction = 'umap')* scale_color_gradientn(colours = c("#92c5de", "#f7f7f7", "#b2182b"))
dev.off()



####### monocle2 
library(monocle)
pd <- new("AnnotatedDataFrame", data = obj1@meta.data)
fd <- new("AnnotatedDataFrame", data = data.frame(gene_short_name = rownames(obj1), row.names = rownames(obj1)))
cds <- newCellDataSet(obj1@assays$RNA@counts, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr=0.1)
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))

markers <- FindAllMarkers(obj1, only.pos = T, logfc.threshold = 0.5)
genes <- markers %>%
    filter(p_val_adj < 0.001) %>%
    arrange(p_val_adj) %>%
    dplyr::slice(1:500)

ordergene <- unique(genes$gene)
cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components = 2, method = "DDRTree")
cds <- orderCells(cds)
cds <- orderCells(cds, root_state = '5')

custom_colors <- c("#377EB8", "#4DAF4A", "#FF7F00")
p1 <- plot_cell_trajectory(cds, 
                           color_by = "Pseudotime", 
                           cell_size = 0.2, 
                           show_tree = F, 
                           show_branch_points = F) + 
    scale_color_gradientn(colors = c("#09030C","#492960", "#B8504E", "#E3A748", "#FCF7AF"))

p2 <- plot_cell_trajectory(cds, 
                           color_by = "Condition", 
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
sam = '3000'
pdf(paste0(sam, '_monocle.pdf'), height = 8, width = 8)
print((p1+p2)/(p3+p3))
dev.off()


my_color_palette <- colorRampPalette(c("#92c5de", "#f7f7f7", "#b2182b"))(100)
genes <- rev(c("UPK1A", "CEBPA", "FXYD3", "IFI27", "EPCAM", "HLA-A", "HLA-B", "KRT18", "FOXA1","GATA3", "CEBPB", "SOX4", "KRT20","ERBB3", "PPARG","JUN", "KRT17", "JUND", "JUNB", "CEBPD", "KRT5", "FOS", "FOSB", "CD44", "RUNX1", "KRT14", "TCF7L2", "KRT16", "CUX1", "TP63","HES2","FOSL1"))
pdf('select_gene_expression_track_monocle2.pdf', width = 6, height = 6)
plot_pseudotime_heatmap(cds[genes,], 
                        #num_clusters = 4, 
                        cluster_rows = F,
                        cores = 1, 
                        return_heatmap = T,
                        show_rownames = TRUE,
                        hmcols = my_color_palette)
dev.off()





##   find degs markers local script
DEG_wilcox2 <- function(samplesCluster, mat){
    library(future.apply)
    
    DEG <- list()
    if(length(samplesCluster) < 3){
        n <- 1
    }else{
        n <- length(samplesCluster)
    }
    
    options(future.globals.maxSize= 8912896000)
    plan(multisession, workers=2)
    
    for(i in 1:n){
        selectSamples <- samplesCluster[[i]]
        otherSamples <- colnames(mat)[!(colnames(mat)%in%selectSamples)]
        
        DEG[[i]] <- matrix(0, nrow = nrow(mat), ncol = 6)
        rownames(DEG[[i]]) <- rownames(mat)
        colnames(DEG[[i]]) <- c("Stat", "Pval", "MeanSelect", "MeanOther", "FC", "FDR")
        
        wiltest <- future_lapply(seq(nrow(mat)), function(x){
            wilcox.test(mat[x, selectSamples], mat[x, otherSamples], alternative = "two.sided")
        })
        
        DEG[[i]][,1] <- unlist(lapply(wiltest, function(x) x$statistic))
        DEG[[i]][,2] <- unlist(lapply(wiltest, function(x) x$p.value))
        DEG[[i]][,3] <- rowMeans(mat[, selectSamples])
        DEG[[i]][,4] <- rowMeans(mat[, otherSamples])
        DEG[[i]][,5] <- log2(DEG[[i]][,3]/DEG[[i]][,4])
        
        
        DEG[[i]][,6] <- p.adjust(p = DEG[[i]][,"Pval"], method = "fdr")
        DEG[[i]] <- DEG[[i]][order(DEG[[i]][,5], decreasing = TRUE),]
        
        rm(selectSamples, otherSamples)
    }
    
    if(length(samplesCluster) < 3){
        names(DEG) <- names(samplesCluster)[1]
    }else{
        names(DEG) <- names(samplesCluster)
    }
    
    return(DEG)
}
sc <- readRDS('3000_new.rds')
for (i in unique(sc$Condition)) {
    sc <- readRDS('3000_new.rds')
    sc$Condition[sc$Condition != i] <- paste0(i,'_refer')
    samplesCluster = list()
    for (subclo in sort(unique(sc$Condition))) {
        samplesCluster[[subclo]] = rownames(sc@meta.data)[sc@meta.data$Condition == subclo]
    }
    
    #options(future.globals.maxSize= 8912896000)
    degs_raw = DEG_wilcox2(samplesCluster, as.matrix(sc@assays$RNA@data))
    degs_raw = data.frame(degs_raw)
    degs_raw$gene <- rownames(degs_raw)
    degs = subset(degs_raw, abs(degs_raw[[paste0(i,'.FC')]]) > 1 & degs_raw[[paste0(i,'.Pval')]] < 0.01 & degs_raw[[paste0(i,'.MeanSelect')]])
    degs[paste0(i,'_TF')] <- rownames(degs) %in% tf$V1
    degs <- degs[order(-degs[[paste0(i,'.FC')]]), ][1:200,]
    openxlsx::write.xlsx(degs, paste0(i, '_degs.xlsx'))
}







