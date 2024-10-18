### plot gene expression or score distribution of mouse spatial data
suppressPackageStartupMessages({
    library(Seurat)
    library(MuDataSeurat)
    library(ggplot2)
    library(gridExtra)
}
)
setwd('~/Desktop/')
PlotList <- list()
name_list <- list(
    c("6w", 0.9),
    c("8w", 0.7),
    c("12w", 1),
    c("24w", 0.8)
)


genes <- c('Carmn', 'Lmod1', 'Myocd','Spp1','Apoe','Postn','Krt17','Epcam')


#load data and add x,y 
for (sam in 1:length(name_list)){
    i = name_list[[sam]][1]
    
    obj <- ReadH5AD(paste0('~/Desktop/BCa/Mouse_st/BinSize/sc_ad_mouse_', i, '_anno.h5ad'))
    obj@assays$Spatial <- obj@assays$RNA
    coord.df = data.frame(
        x = obj@reductions$spatial@cell.embeddings[, 'spatial_2'],
        y = obj@reductions$spatial@cell.embeddings[, 'spatial_1'],
        stringsAsFactors = FALSE
    )
    rownames(coord.df) = rownames(obj@reductions$spatial@cell.embeddings)
    colnames(coord.df) <- c('x', 'y')
    coord.df <- coord.df[rownames(obj@meta.data),]
    obj@images$image =  new(Class = 'SlideSeq', assay = "Spatial",key = "image_", coordinates = coord.df)
    
    DefaultAssay(obj) <- 'Spatial'
    obj <- NormalizeData(obj, assay = 'Spatial')
    #cols = c("#92c5de", "#f7f7f7", "#b2182b")
    cols =  c('#FFeFe7','#F7AA77','#D8525F','#992C79')
    
    
    plot_data <- as.data.frame(t(obj@assays$Spatial@data[genes,]))
    
    plot_data$X <- obj$X
    plot_data$Y <- -obj$Y
    
    dot_size = as.numeric(name_list[[sam]][2])
    
    for (j in seq_along(colnames(plot_data)[1:8])) {
        gene <- colnames(plot_data)[j]
        if (gene == 'Apoe'){
            PlotList[[paste0(i,'_',j)]] <- ggplot(plot_data, aes(x = X, y = Y)) +  
                geom_point(aes_string(fill = gene), size = dot_size, stroke = NA, shape = 21) +
                scale_fill_gradientn(colours = cols, limits = c(0, 5),  oob = scales::squish) + 
                labs(title = paste0('BCa_mouse_', i)) +
                theme_void() + 
                theme(
                 #  panel.background = element_rect(fill = "#E8E8E8", colour = NA),
                 #  plot.background = element_rect(fill = "#E8E8E8", colour = NA), 
                    legend.position = "right",
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    aspect.ratio = (max(obj$Y)-min(obj$Y))/(max(obj$X)-min(obj$X))
                )
        } else {
            PlotList[[paste0(i,'_',j)]] <- ggplot(plot_data, aes(x = X, y = Y)) + 
                geom_point(aes_string(fill = gene), size = dot_size, stroke = NA, shape = 21) +
                scale_fill_gradientn(colours = cols, limits = c(0, 3),  oob = scales::squish) +
                labs(title = paste0('BCa_mouse_', i)) +
                theme_void() + 
                theme( 
                  #  panel.background = element_rect(fill = "#E8E8E8", colour = NA),
                  #  plot.background = element_rect(fill = "#E8E8E8", colour = NA),
                    legend.position = "right",
                    panel.grid = element_blank(),
                    plot.title = element_text(hjust = 0.5),
                    aspect.ratio = (max(obj$Y)-min(obj$Y))/(max(obj$X)-min(obj$X))
                )
        }
    }
}

pdf('mouse_SD.pdf',  width = 40, height = 20)
cowplot::plot_grid(plotlist = PlotList, ncol = 8)  
dev.off()

type_list <- list(
    SD11 = list(c("Acer2", "Bhmt", "Foxa1", "Cyp4f8", "Upk1a", "Ina", "Uca1", "Fam13a", "Tox3", "Limch1", "Tmem45b", "Ctse", "Inpp4b", "B2m", "Fabp5", "Pla2g2a", "Gpr160", "Bambi", "Ndufc2-Kctd14", "St3gal4", "Slitrk6", "Vegfa", "Gnl3", "Pnck", "Scd", "Egln3", "Abcd3", "Nek7", "Upk2", "Srebf1", "Pof1b", "Dapk1", "Idh1", "Ccser1", "Rnf128", "Hpgd", "Grhl1", "Fer1l4", "Snx31", "Cyp4f22")),
    SD2  = list(c("Krt17", "Vim", "Gdf15", "Hmgcs2", "Crisp3", "Cracr2b", "Prss23", "Krt4", "Crip1", "Myeov", "Krt5", "Bcam", "Ccn2")),
    SD10 = list(c("Krt6a", "Myeov", "S100a9", "Nr4a1", "Egr1", "Fosb", "Sod3", "Ifitm2", "Maml2", "Cxcl14", "Cd74", "Rgs5", "Cd93", "Krt5", "Slpi", "Stab1", "Vim", "A2m", "Sparc", "Cav1", "Cavin1", "Thy1", "Tgm2", "Plat", "Crip1", "Pmepa1", "Timp2", "Tymp", "Lgals1", "Ifi6", "Msn")),
    SD0  = list(c("Timp1", "Krt4", "S100a9", "Krt5", "Krt15", "Slpi", "Dusp1", "Selenop", "Cxcl14", "Fkbp5", "Ccn2", "Cd74", "A2m", "Vim", "Ifitm2", "Tgfbi", "Hla-b", "Sparc", "Hla-dqb1", "Timp2", "Rgs5", "Ctsa", "Ctsd", "Rnase1", "Lgals1", "Cavin1", "Apoe", "Rhoc", "Tgm2", "Crip1", "Tmsb10", "Ndufs6")),
    SD3  = list(c("Emb", "Itm2c", "Siae", "Bhmt", "Degs1", "Gata3", "Snhg25", "Snhg19", "Evl", "Fam3b", "Acox3", "Ints1", "Psmg1", "Casp8ap2", "Mrpl35", "Bri3bp", "Rabgap1l", "Trit1", "Dctn1", "Sept9", "Rabl2b", "Fasn", "Wdr5", "Ndufa10")),
    SD14 = list(c("Cyp1b1", "Cyp1a1", "Apoe", "Fn1", "Fabp4", "Aebp1", "Cav1", "Sparc", "Lbh", "A2m", "Polr2j4", "Pkp1", "Timp1", "Upk3bl1", "Dmkn", "Itm2c", "Pygl", "Timp3", "Krt17", "Spock1", "Actn1", "Itgav", "Mgll", "Cpt1a", "Man1a1", "Csrp1", "Arl2", "Ano1", "Krt16", "Snhg25", "Krt5", "Pacs1", "Tnni2", "Thbs2", "Plec", "Cyb561a3", "Pdxk", "Hla-b", "Cd59", "Slc29a1", "Tymp", "Pmepa1", "Slc7a5", "Eno1", "S100a9")),
    Basal = list(c("Cd44", "Cdh3", "Krt1", "Krt14", "Krt16", "Krt5", "Krt6a", "Krt6b", "Krt6c")),
    Luminal = list(c("Cyp2j2", "Erbb2", "Erbb3", "Fgfr3", "Foxa1", "Gata3", "Gpx2", "Krt18", "Krt19", "Krt20", "Krt7", "Krt8", "Pparg", "Xbp1", "Upk1a", "Upk2"))
)
for (type in names(type_list)) {
    score_name <- paste0(type, "_score")
    obj <- AddModuleScore(obj, features = type_list[[type]], ctrl = 100, name = score_name)
    colnames(obj@meta.data) <- sub("score1$", "score", colnames(obj@meta.data))
}


