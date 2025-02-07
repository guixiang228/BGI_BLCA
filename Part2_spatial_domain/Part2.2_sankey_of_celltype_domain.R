#### Sankey plot to show the proportion of cell types in each domain
library(tidyverse)
library(Seurat)
library(ggplot2)
library(MuDataSeurat)
library(networkD3)
library(dplyr)
library(purrr)
library(future)
library(igraph)
library(htmlwidgets)
library(anndataR)
library(SeuratObject)
library(pandoc)
library(knitr)
library(rmarkdown)

seurat_obj <- ReadH5AD('allsample_domain.h5ad')
workdir <- "./sankey/"
setwd(workdir)

seurat_obj@meta.data$annotated_cluster <- as.factor(seurat_obj@meta.data$annotated_cluster)
meta_data <- seurat_obj@meta.data


meta_data <- meta_data %>%
  mutate(spatial_domain = as.character(spatial_domain)) %>%
  filter(!is.na(spatial_domain))

meta_data <- meta_data %>%
  mutate(annotated_cluster = as.character(annotated_cluster)) %>%
  filter(!is.na(annotated_cluster))

# Calculate the number of different celltypes in each domain
counts <- meta_data %>%
  count(annotated_cluster, spatial_domain) %>%
  mutate(Weight = n)

sankey_data <- counts %>%
  mutate(Source = annotated_cluster,
         Target = spatial_domain,
         Weight = n) %>%
  select(Source, Target, Weight)
nodes<- data.frame(name = c(as.character(sankey_data$Source), 
                             as.character(sankey_data$Target)) %>% unique())
                          
print(2)

nodes <- data.frame(
  name = c( "B cell", "DCs", "Endothelial", "Fibroblast", "Macrophage", "Mastcell", "Myocyte", "Normal Epithelial", "Pericyte", "Plasmocyte", "T cell", "Tumor cell",
            '0', '1', '10', '11', '12', '13', '14', '2', '3', '4', '5', '6', '7', '8', '9'),
  color = c('#FFBBFF', '#81d8ae', '#607d3a', '#FFD700', '#87CEFF', '#AB82FF', '#3399FF', '#663300', '#FFD39B', '#6bd155', '#B23AEE', '#CD3278',
             '#FF8C00', '#d5492f',  '#6bd155', '#d34a76', '#c9d754', '#d04dc7', '#00FFFF','#683ec2', '#ede737', '#6d76cb', '#ce9d3f', '#81357a', '#d3c3a4', '#3c2f5a', '#b96f49'),
          stringsAsFactors = FALSE
)

sankey_data <- data.frame( Source = as.character(sapply(sankey_data$Source, function(x) if (x %in% nodes$name) x else NA())),
                           Target = as.character(sapply(sankey_data$Target, function(x) if (x %in% nodes$name) x else NA())),
                           Weight = sankey_data$Weight
                          )


sankey_data$IDsource = match(sankey_data$Source, nodes$name)-1
sankey_data$IDtarget = match(sankey_data$Target, nodes$name)-1


my_colors <- JS("d3.scaleOrdinal().range(['#FFBBFF', '#81d8ae', '#607d3a', '#FFD700', '#87CEFF', '#AB82FF', '#3399FF', '#663300', '#FFD39B', '#6bd155', '#B23AEE', '#CD3278',
                                          '#FF8C00', '#d5492f', '#6bd155', '#d34a76', '#c9d754', '#d04dc7', '#00FFFF','#683ec2', '#ede737', '#6d76cb', '#ce9d3f', '#81357a', '#d3c3a4', '#3c2f5a', '#b96f49'])")
print(nodes)
print(sankey_data)
print(nodes)
print(sankey_data)
print(my_colors)

p <- sankeyNetwork(Links = sankey_data,
                   Nodes = nodes,
                   Source = "IDsource",  
                   Target = "IDtarget",  
                   Value = "Weight",   
                   NodeID = "name",    
                   LinkGroup = 'Source',
                   nodeWidth = 25,     
                   fontSize = 10,      
                   nodePadding = 3,   
                   sinksRight = FALSE,
                   NodeGroup = "name",
                   colourScale = my_colors
                   )
saveRDS(p, file = "sankey_celltype_to_domain.RDS")
rmarkdown::find_pandoc(dir='./software/anaconda/envs/R4.2.1/bin')
saveWidget(p, "sankey_domain_to_celltype.html", selfcontained = TRUE)
