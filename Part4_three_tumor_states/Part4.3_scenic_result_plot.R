### visualization of scenic result
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)
library(Seurat)

loom = open_loom('out_SCENIC.loom')

regulons_incidMat<-get_regulons(loom,column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
regulons<-regulonsToGeneLists(regulons_incidMat)
regulonAUC<-get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAucThresholds<-get_regulon_thresholds(loom)
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])

embeddings<-get_embeddings(loom)
close_loom(loom)

sub_regulonAUC = regulonAUC[,match(colnames(obj),colnames(regulonAUC))]


cellClusters <- data.frame(row.names = colnames(obj), 
                           seurat_clusters = as.character(obj$condition))

cellTypes <- data.frame(row.names = colnames(obj), 
                        celltype = obj$condition)

save(sub_regulonAUC, cellTypes, cellClusters, obj,
     file = 'for_rss_and_visual.Rdata')


selectedResolution <- "celltype"
cellsPerGroup <- split(rownames(cellTypes), cellTypes[, selectedResolution])

# Remove extended regulatory factors
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),] 

# Calculate the average expression of each group
regulonActivity_byGroup <- sapply(cellsPerGroup, function(cells) 
    rowMeans(getAUC(sub_regulonAUC)[,cells]))
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup), center = TRUE, scale = TRUE)) 
regulonActivity_byGroup_Scaled <- na.omit(regulonActivity_byGroup_Scaled)

pdf('SCENIC_3000.pdf', width = 4, height = 20)
Heatmap(
    regulonActivity_byGroup_Scaled,
    name                            = "z-score",
    col                             = c("#92c5de", "#f7f7f7", "#b2182b"),
    show_row_names                  = TRUE,
    show_column_names               = TRUE,
    row_names_gp                    = gpar(fontsize = 6), 
    clustering_method_rows          = "ward.D2",
    clustering_method_columns       = "ward.D2",
    row_title_rot                   = 0, 
    cluster_rows                    = TRUE,
    cluster_row_slices              = FALSE,
    cluster_columns                 = FALSE
)
dev.off()


# Calculate Regulatory Factor Specific Score (RSS)
rss <- calcRSS(AUC = getAUC(sub_regulonAUC), cellAnnotation = cellTypes[colnames(sub_regulonAUC), selectedResolution])
rss <- na.omit(rss)

pdf('rssPlot_z=1.pdf',height=20, width = 4)
rss <- rss[,c(2,3,1)]
print(plotRSS(rss,zThreshold = 1,cluster_columns = FALSE,order_rows = TRUE,thr=0.1,varName = "condition",col.low = '#92c5de',col.mid = '#f7f7f7',col.high = '#b2182b'))
dev.off()

# selection process and plot
genes <- c("STAT5A", "STAT2", "RUNX3", "TCF7L2", "HES2", "RUNX1", "FOSL1", "TP63", "ETV6", "TFAP2A", "CUX1", "FOS", "FOSB", "HMGA2", "JUNB", "IRF6", "JUN", "JUND", "CEBPD", "CEBPB", "HMGA1", "PPARG", "IRF1","STAT1", "SOX4", "CEBPA", "FOXA1")
gene1 <- paste0(genes,'(+)')
rss_sel <- rss[rownames(rss) %in% gene1,]
df1 <- plotRSS(rss,zThreshold = 1,cluster_columns = FALSE,order_rows = TRUE,thr=0.1,varName = "condition",col.low = '#92c5de',col.mid = '#f7f7f7',col.high = '#b2182b')
df <- df1$df[df1$df$Topic %in% gene1, ]
df$condition <- as.character(df$condition)
df$condition[df$condition == "Basel"] <- "Basal"
df$condition <- as.factor(df$condition)
levels(df$condition)[levels(df$condition) == "Basel"] <- "Basal"
table(df$condition)

pdf('rssPlot_selection.pdf',height=6, width = 4)
ggplot(df, aes(x = condition, y = Topic)) +
    geom_point(aes(size = RSS, color = Z)) +
    scale_size_continuous(range = c(1,4)) +
    scale_color_gradientn(colors = c('#92c5de', '#f7f7f7', '#b2182b'))+
    labs(x = "Group", y = "TF", size = "RSS", color = "Z") +
    theme_bw() +
    guides(
        size = guide_legend(order = 1), 
        color = guide_colorbar(order = 2)
    )
dev.off()
