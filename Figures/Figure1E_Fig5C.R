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

# Fig1E
DotPlot(Clusters_sub, idents = "FMC", features = "Tlr3", group.by = "celltype_broad", dot.scale = 10) +
 theme(axis.title = element_blank()) +
 scale_y_discrete(limits = rev(c("TAM", "Oligodendrocyte", "T cell", "Tumor cell"))) +
 theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 20, vjust = 0.5, face = "italic")) +
 theme(axis.text.y = element_text(size = 20))

# Fig5C
features = c("Ccl2", "Ccl3", "Ccl7", "Ccl12", "Pf4", "Cxcl9", "Cxcl10", "Cxcl13")
DotPlot(Clusters_sub, features = features, group.by = "celltype_broad", idents = "FMC") +
 theme(axis.title = element_blank()) +
 scale_y_discrete(limits = rev(c("TAM", "Oligodendrocyte", "T cell", "Tumor cell", "Astrocyte"))) +
 theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 20)) +
 theme(axis.text.y = element_text(size = 20)) +
 coord_flip()
