### Pathway enrichment
library(reshape)
library(reshape2)
library(tidyverse)
library(dplyr)
setwd('./GO/')

#### Extract the Pvalue and genes in each pathway
enrichment_all <- read.csv("GO_AllLists.csv")
enrichment_all$names <- paste0(enrichment_all$Category, "_",enrichment_all$GO, "_", enrichment_all$Description)

enrichment_all_dcast <- dcast(enrichment_all, names ~ GeneList, value.var = "LogP", fill = 0)
enrichment_all_dcast <- column_to_rownames(enrichment_all_dcast, var = "names")

enrichment_all_dcast$sigPvalue_min <- apply(enrichment_all_dcast, 1, min)
enrichment_all_dcast$sigPvalue_median <- apply(enrichment_all_dcast, 1, median)
enrichment_all_dcast$sigPvalue_mean <- apply(enrichment_all_dcast, 1, mean)
enrichment_all_dcast$sigPvalue_diff <- enrichment_all_dcast$sigPvalue_mean - enrichment_all_dcast$sigPvalue_min

enrichment_all_Hits <- dcast(enrichment_all, names ~ GeneList, value.var = "Hits", fill = 0)
enrichment_all_Hits <- column_to_rownames(enrichment_all_Hits, var = "names")
enrichment_all_Hits <- enrichment_all_Hits[rownames(enrichment_all_dcast),]

enrichment_all_dcast_combin <- cbind(enrichment_all_dcast, enrichment_all_Hits[rownames(enrichment_all_dcast),])
enrichment_all_dcast_combin$category <- sapply(rownames(enrichment_all_dcast_combin), function(x) unlist(strsplit(x,"_", fixed = TRUE))[1])
enrichment_all_dcast_combin$GO <- sapply(rownames(enrichment_all_dcast_combin), function(x) unlist(strsplit(x,"_", fixed = TRUE))[2])
enrichment_all_dcast_combin$Description <- sapply(rownames(enrichment_all_dcast_combin), function(x) unlist(strsplit(x,"_", fixed = TRUE))[3])

write.csv(enrichment_all_dcast_combin, "GO_AllLists_dcast.csv")

### select the top 10 pathways with smallest P value
names(enrichment_all_dcast_combin)[duplicated(names(enrichment_all_dcast_combin))] <- paste0(names(enrichment_all_dcast_combin)[duplicated(names(enrichment_all_dcast_combin))], "_dup")
domain0_min <- enrichment_all_dcast_combin %>%
  arrange(domain0) %>%
  slice(1:10)

domain10_min <- enrichment_all_dcast_combin %>%
  arrange(domain10) %>%
  slice(1:10)

domain11_min <- enrichment_all_dcast_combin %>%
  arrange(domain11) %>%
  slice(1:10)

domain14_min <- enrichment_all_dcast_combin %>%
  arrange(domain14) %>%
  slice(1:10)

domain2_min <- enrichment_all_dcast_combin %>%
  arrange(domain2) %>%
  slice(1:10)

domain3_min <- enrichment_all_dcast_combin %>%
  arrange(domain3) %>%
  slice(1:10)


enrich_filt <- bind_rows(domain0_min,domain10_min,domain11_min,domain14_min,domain2_min,domain3_min) %>%
  distinct()  

## to plot the bubble plots
enrichment_all_filt <- enrichment_all[enrichment_all$GO%in%enrich_filt$GO,]
enrichment_all_filt$LogP <- -enrichment_all_filt$LogP
enrichment_all_filt$Description <- gsub("HALLMARK ", "", enrichment_all_filt$Description, fixed = TRUE)
enrichment_all_filt$Description <- str_to_title(enrichment_all_filt$Description)

new_order <- c('Signaling By Rho Gtpases, Miro Gtpases And Rhobtb3','Regulation Of Dna Metabolic Process','Regulation Of Cell Projection Organization','Programmed Cell Death','Neuron Projection Morphogenesis','Neuron Projection Development','Cell Projection Morphogenesis','Cell Morphogenesis','Apoptotic Execution Phase',
'Thermogenesis','Respiratory Electron Transport','Parkinson Disease','Oxidative Phosphorylation','Huntington Disease','Electron Transport Chain Oxphos System In Mitochondria','Diabetic Cardiomyopathy','Chemical Carcinogenesis - Reactive Oxygen Species','Cellular Response To Chemical Stress','Aerobic Respiration And Respiratory Electron Transport',
'Nuclear Receptors Meta Pathway','Nonalcoholic Fatty Liver Disease','Non-Alcoholic Fatty Liver Disease','Interferon Alpha/Beta Signaling','Alzheimer Disease',
'Signaling By Receptor Tyrosine Kinases','Pid Ap1 Pathway','Neutrophil Degranulation','Negative Regulation Of Catalytic Activity','Lysosome','Interferon Signaling',
'Energy Derivation By Oxidation Of Organic Compounds','Cellular Respiration',
'Vegfa Vegfr2 Signaling','Supramolecular Fiber Organization','Regulation Of Proteolysis','Naba Matrisome Associated','Muscle Structure Development','Keratinocyte Differentiation','Epithelial Cell Differentiation','Epithelial Cell Development','Epidermis Development','Cytokine Signaling In Immune System','Actin Filament-Based Process','Actin Cytoskeleton Organization')
new_order <- rev(new_order)
enrichment_all_filt$Description <- factor(enrichment_all_filt$Description, levels = new_order)

new_order <- c('domain11','domain2','domain10','domain0','domain3','domain14')
enrichment_all_filt$GeneList <- factor(enrichment_all_filt$GeneList, levels = new_order)

pdf("GO_AllLists_filter_bubble.pdf", width = 10, height = 14)
ggplot(enrichment_all_filt, aes(GeneList, Description)) +
  geom_point(aes(fill = LogP, size = Enrichment), alpha = 0.9, pch = 21, colour = "gray25") +
  scale_fill_gradient(low = '#3288bd', high = 'firebrick3') + 
  scale_size_continuous(range = c(1, 10)) + 
  labs(y = '', x = '', fill = '-log10pvalue', size = 'Enrichment') + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, vjust = 1, hjust = 1),
        axis.text = element_text(colour = "#000000", size = 17))
dev.off()