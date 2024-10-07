### infer the trajectory path of myocyte and fibroblast subtype by monocle2
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle)
library(MuDataSeurat)

x <- load('merge_CAF_5Myo.RData')
seurat_comb_singlet <- get(x)
Idents(seurat_comb_singlet) = seurat_comb_singlet$celltype

seurat_comb_singlet@meta.data$idents <- Idents(seurat_comb_singlet)
current.cluster.ids <- c('7','FBXO32+CAF','Muscle_like_CAF','Myocyte','POSTN+CAF','SLC14A1+CAF')
new.cluster.ids <- c('ABCC4+myCAF','FBXO32+myCAF','MYOCD+myCAF','Myocyte','POSTN+myCAF','SLC14A1+myCAF')
seurat_comb_singlet@meta.data$celltype = plyr::mapvalues(x = seurat_comb_singlet@meta.data[,"idents"], from = current.cluster.ids, to = new.cluster.ids)
Idents(seurat_comb_singlet) <- seurat_comb_singlet@meta.data$celltype
str(seurat_comb_singlet)

expr_matrix <- as(as.matrix(seurat_comb_singlet@assays$RNA@counts),"sparseMatrix")
p_data <- seurat_comb_singlet@meta.data
f_data <- data.frame(gene_short_name=row.names(seurat_comb_singlet), row.names = row.names(seurat_comb_singlet))

pd <- new("AnnotatedDataFrame", data= p_data)
fd <- new("AnnotatedDataFrame", data= f_data)
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
cds <- detectGenes(cds, min_expr=0.1)
print(head(fData(cds)))
save(cds, file = "cds.RData")
load('cds.RData')

deg.cluster <- FindAllMarkers(seurat_comb_singlet)
expressed_genes <- subset(deg.cluster, p_val_adj<0.001)$gene
cds <- setOrderingFilter(cds, expressed_genes)
save(cds, file = "cds_OF.RData")

diff <- differentialGeneTest(cds[expressed_genes,], fullModelFormulaStr="~celltype", core=1)
deg <- subset(diff, qval<0.01)
deg <- deg[order(deg$qval,decreasing = F),]
write.table(deg, file="train.cluster.DEG.xls",col.names=T,row.names=F,sep="\t",quote=F)
ordergene <- rownames(deg)

cds <- setOrderingFilter(cds, ordergene)
cds <- reduceDimension(cds, max_components=2, method="DDRTree")
save(cds, file = "cds_RD.RData")
cds <- orderCells(cds)
save(cds, file = "cds_order.RData")

pdf(paste0(workdir,"monocle.pseudotime.pdf"),width=7,height=7)
plot_cell_trajectory(cds, color_by="Pseudotime", cell_size=0.5, show_backbone=TRUE)
dev.off()
custom_colors <- c(
  "ABCC4+myCAF" = "#377eb8",
  "POSTN+myCAF" = "#ff7f00",
  "SLC14A1+myCAF" = "#984ea3",
  "FBXO32+myCAF" = "#4daf4a",
  "MYOCD+myCAF" = "#f0da3a",
  "Myocyte" = "#e41a1c") 
pdf(paste0(workdir,"monocle.celltype.pdf"),width=7,height=7)
plot_cell_trajectory(cds, 
                          color_by = "celltype", 
                          cell_size = 1,
                          show_backbone = TRUE, 
                          cell_name_size = 10, 
                          cell_link_size = 0.5) +
  scale_color_manual(values = custom_colors) +
  theme(legend.text = element_text(size = 12),
        legend.title = element_text(size = 14)) + 
  guides(color = guide_legend(override.aes = list(size = 3))) 
dev.off()
pdf(paste0(workdir,"monocle.state.pdf"),width=7,height=7)
plot_cell_trajectory(cds, color_by="State", cell_size=0.5, show_backbone=TRUE)
dev.off()

