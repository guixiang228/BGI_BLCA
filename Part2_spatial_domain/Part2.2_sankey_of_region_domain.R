#### Sankey plot to show the proportion of domains in each region
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


seurat_obj <- readRDS('allsample_domain.rds')
meta_data <- seurat_obj@meta.data
meta_data <- meta_data[, !duplicated(names(meta_data))]
meta_data <- meta_data %>%
  filter(spatial_domain %in% c("11", "2", "10", "0", "3", "14"))

meta_data <- meta_data %>%
  mutate(Region = as.character(Region),
         Region = if_else(Region == "Black", NA_character_, Region),
         Region = recode(Region, `Red` = "ET", `Blue` = "IT", `Green` = "ST")) %>%
  filter(!is.na(Region))


counts <- meta_data %>%
  count(spatial_domain, Region) %>%
  mutate(Weight = n)


sankey_data <- counts %>%
  mutate(Source = spatial_domain,
         Target = Region,
         Weight = n) %>%
  select(Source, Target, Weight)


nodes<- data.frame(name = c(as.character(sankey_data$Source), 
                             as.character(sankey_data$Target)) %>% unique())


nodes <- data.frame(
  name = c("11", "2", "10", "0", "3", "14", "ET", "IT", "ST"),
  color = c("#d34a76", "#683ec2", "#6bd155", "#FF8C00", "#ede737", "#00FFFF", "#EF767A", "#456990", "#48C0AA"),
  stringsAsFactors = FALSE
)


sankey_data <- data.frame(Source = as.character(sapply(sankey_data$Source, function(x) if (x %in% nodes$name) x else NA())),
                          Target = as.character(sapply(sankey_data$Target, function(x) if (x %in% nodes$name) x else NA())),
                          Weight = sankey_data$Weight)

sankey_data$IDsource = match(sankey_data$Source, nodes$name)-1
sankey_data$IDtarget = match(sankey_data$Target, nodes$name)-1

my_colors <- JS("d3.scaleOrdinal().range(['#6bd155',  '#683ec2','#d34a76', '#ede737','#FF8C00','#00FFFF', '#EF767A', '#456990', '#48C0AA'])")
print(nodes)
print(sankey_data)

p <- sankeyNetwork(Links = sankey_data,
                   Nodes = nodes,
                   Source = "IDsource",  
                   Target = "IDtarget",  
                   Value = "Weight",   
                   NodeID = "name",    
                   LinkGroup = 'Source',
                   nodeWidth = 20,     
                   fontSize = 15,      
                   nodePadding = 30,   
                   sinksRight = FALSE,
                   NodeGroup = "name",
                   colourScale = my_colors
                   )

library(knitr)
library(rmarkdown)
rmarkdown::find_pandoc(dir='./anaconda/envs/R4.2.1/bin')
saveWidget(p, "sankey_domain_to_region.html", selfcontained = TRUE)

