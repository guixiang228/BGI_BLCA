### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(future)
library(Cairo)
library(harmony)
options(bitmapType = "cairo")

dim.usage <- 30
workdir <- "./intersectHVG/"
setwd(workdir)

x <- load('fibroblast.RData')
seurat_comb_singlet <- get(x)

#### use common genes among samples
geneset_test <- c()
for (i in c("HBCP1","HBCP2","HBCP3A","HBCP3B","HBCP3C","HBCP4A","HBCP4B","HBCP5A","HBCP5B","HBCP6","HBCP7","HBCP8","HBCP9","HBCP10","HBCP11","HBCP12","HBCP13","HBCP14","HBCP15","HBCP16","HBCP17","HBCP18","HBCP19","HBCP20","HBCP21","HBCP22")) {
    sc_test <- subset(seurat_comb_singlet, Patients == i)
    sc_test <- NormalizeData(sc_test)
    sc_test <- FindVariableFeatures(sc_test, selection.method = "vst", nfeatures = 2000)
    geneset_test <- c(geneset_test, VariableFeatures(sc_test))
}
geneset <- c()
for (i in geneset_test) {
    if (table(geneset_test)[i] >= 15) {
        geneset <- c(geneset, i)
    }
}
geneset <- unique(geneset)


seurat_comb_singlet <- RunPCA(seurat_comb_singlet, features = geneset)
seurat_comb_singlet <- RunUMAP(seurat_comb_singlet, dims = 1:dim.usage, reduction.use='pca')

# #### cluster all the cells  ####
seurat_comb_singlet <- FindNeighbors(seurat_comb_singlet, reduction.use='pca')
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.2)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.4)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.5)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.6)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 0.8)
seurat_comb_singlet <- FindClusters(seurat_comb_singlet, resolution = 1)

save(seurat_comb_singlet, file = paste0(workdir, "_cluster.RData"))



x <- load(file = paste0(workdir, "_cluster.RData"))
seurat_obj <- get(x)

replacement_rules <- c(
  '0' = 'ABL2+CAF',
  '1' = 'SPINT2+CAF',
  '2' = 'FBXO32+myCAF',
  '3' = 'ZBTB16+CAF',
  '4' = 'SLC14A1+myCAF',
  '5' = 'ZBTB16+CAF',
  '6' = 'NAV3+CAF',
  '7' = 'ABCC4+myCAF',
  '8' = 'FBLN1+CAF',
  '9' = 'ITGA8+myCAF',
  '10' = 'DEPTOR+CAF',
  '11' = 'IL6+iCAF',
  '12' = 'MYL9+myCAF',
  '13' = 'ABL2+CAF',
  '14' = 'MYOCD+myCAF',
  '15' = 'POSTN+myCAF'
)

seurat_obj@meta.data$fib_type <- seurat_obj@meta.data$RNA_snn_res.1
seurat_obj@meta.data$fib_type <- ifelse(seurat_obj@meta.data$RNA_snn_res.1 %in% names(replacement_rules),
                                        replacement_rules[seurat_obj@meta.data$RNA_snn_res.1],
                                        seurat_obj@meta.data$fib_type)
save(seurat_obj, file="fib_type.RData")

colors <- c('ABL2+CAF' = "#FFE4B5", 'SPINT2+CAF' = "#1f77b4", 'FBXO32+myCAF' = "#ff7f0e", 'ZBTB16+CAF' = "#2ca02c",
            'SLC14A1+myCAF' = "#d62728", 'NAV3+CAF' = "#9467bd", 'ABCC4+myCAF' = "#8c564b", 'FBLN1+CAF' = "#e377c2",
            'ITGA8+myCAF' = "#e86502", 'DEPTOR+CAF' = "#bcbd22", 'IL6+iCAF' = "#17becf", 'MYL9+myCAF' = "#ffbb78",
            'MYOCD+myCAF' = "#f7b6d2", 'POSTN+myCAF' = "#c5b0d5")

p <- DimPlot(seurat_obj, reduction = "umap", group.by = "fib_type", label = TRUE, pt.size = 1) + 
  scale_color_manual(values = colors)
ggsave("fib_type_cluster.pdf", plot = p, width = 10, height = 8)



