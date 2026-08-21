library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(qs)
library(future)
library(extrafont)
library(cowplot)
library(RColorBrewer)
library(nichenetr)
library(singleCellTK)
library(reticulate)
library(ggrepel)
library(tidyverse)

TAM_Clusters_sub <- qread(paste0(path, "TAM_Clusters_sub.qs"))

# supplementary Figure1C
features = c("TLR3")

DotPlot(TAM_Clusters_sub, features = features, group.by = "celltype_combine", dot.scale = 10, idents = c("PA")) +
    scale_y_discrete(limits = rev(levels(factor(TAM_Clusters_sub$celltype_combine))[1:3])) +
    theme(axis.text.y = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 90, face = "italic", vjust = 0.5, hjust = 1, size = 20)) +
    theme(axis.title = element_blank())


# Figure7G
DotPlot(TAM_Clusters_sub, features = "IFIH1", group.by = "celltype_combine", dot.scale = 10, idents = c("PA")) +
    scale_y_discrete(limits = rev(levels(factor(TAM_Clusters_sub$celltype_combine))[1:3])) +
    theme(axis.text.y = element_text(size = 20)) +
    theme(axis.text.x = element_text(angle = 90, face = "italic", vjust = 0.5, hjust = 1, size = 20)) +
    theme(axis.title = element_blank())
