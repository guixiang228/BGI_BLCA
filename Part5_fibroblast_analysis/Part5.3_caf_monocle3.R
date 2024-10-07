### infer the fibroblast trajectory path by monocle3
suppressPackageStartupMessages({
    library(monocle3)
    library(Seurat)
    library(MuDataSeurat)
    library(patchwork)
    library(SeuratWrappers)
})
rm(list=ls())
sam = 'seurat_wrapper_'

load('fib_type.RData')
Idents(seurat_obj)  <- 'RNA_snn_res.1'
new_cluster_names <- c(
    '0'  = 'ABL2+CAF',
    '1'  = 'SPINT2+CAF',
    '2'  = 'FBXO32+myCAF',
    '3'  = 'ZBTB16+CAF',
    '4'  = 'SLC14A1+myCAF',
    '5'  = 'ZBTB16+CAF',
    '6'  = 'NAV3+CAF',
    '7'  = 'ABCC4+myCAF',
    '8'  = 'FBLN1+CAF',
    '9'  = 'ITGA8+myCAF',
    '10' = 'DEPTOR+CAF',
    '11' = 'IL6+iCAF',
    '12' = 'MYL9+myCAF',
    '13' = 'ABL2+CAF',
    '14' = 'MYOCD+myCAF',
    '15' = 'POSTN+myCAF'
)
seurat_obj <- RenameIdents(seurat_obj, new_cluster_names)
seurat_obj$fib_type <- as.factor(seurat_obj@active.ident)

cds <- as.cell_data_set(seurat_obj, assay = 'RNA', reductions = 'umap', group.by = 'fib_type')
fData(cds)$gene_short_name <- rownames(fData(cds))
recreate.partitions <- c(rep(1, length(cds@colData@rownames)))
names(recreate.partitions) <- cds@colData@rownames
recreate.partitions <- as.factor(recreate.partitions)

cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partitions
cds@clusters@listData[["UMAP"]][["clusters"]] <- seurat_obj$fib_type
cds@int_colData@listData[["reducedDims"]]@listData[["UMAP"]] <- seurat_obj@reductions$umap@cell.embeddings


cds <- learn_graph(cds, use_partition = F)
earliest_principal_node <- function(cds, time_bin='ITGA8+myCAF'){
    cell_ids <- which(colData(cds)[, "fib_type"] == time_bin)
    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
    root_pr_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(table(closest_vertex[cell_ids,]))))]
    root_pr_nodes
}
cds <- order_cells(cds, reduction_method = "UMAP", root_pr_nodes = earliest_principal_node(cds))


pdf(paste0(sam,'monocle3.pdf'),width = 6, height = 5)
plot_cells(cds, label_groups_by_cluster=T, color_cells_by = "fib_type", cell_size = 0.3, group_label_size = 6, 
           label_cell_groups = F, label_leaves = FALSE, label_branch_points = FALSE) + 
    scale_color_manual(values = c('ABL2+CAF'    = "#FFE4B5", 'SPINT2+CAF' = "#1f77b4", 'FBXO32+myCAF' = "#ff7f0e", 'ZBTB16+CAF'= "#2ca02c",
                                  'SLC14A1+myCAF' = "#d62728", 'NAV3+CAF'   = "#9467bd", 'ABCC4+myCAF'  = "#8c564b", 'FBLN1+CAF' = "#e377c2",
                                  'ITGA8+myCAF'  = "#e86502", 'DEPTOR+CAF' = "#bcbd22", 'IL6+iCAF'   = "#17becf", 'MYL9+myCAF' = "#ffbb78",
                                  'ZBTB16+CAF'   = "#ff9896", 'MYOCD+myCAF' = "#f7b6d2", 'POSTN+myCAF' = "#c5b0d5"))
plot_cells(cds, color_cells_by = "pseudotime", label_groups_by_cluster = T, label_branch_points = F, label_roots = T, label_leaves = F)
dev.off()


DimPlot(obj, reduction = "umap", label = TRUE)
