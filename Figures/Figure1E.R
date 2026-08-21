library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(qs)
library(future)
library(cowplot)
library(RColorBrewer)

Clusters_sub <- qread(paste0(path, "OPG/Clusters_sub.qs"))

DotPlot(Clusters_sub, idents = "FMC", features = features, group.by = "celltype_broad", dot.scale = 10) +
 theme(axis.title = element_blank()) +
 scale_y_discrete(limits = rev(c("TAM", "Oligodendrocyte", "T cell", "Tumor cell"))) +
 theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, vjust = 0.5, face = "italic")) +
 theme(axis.text.y = element_text(size = 20))
