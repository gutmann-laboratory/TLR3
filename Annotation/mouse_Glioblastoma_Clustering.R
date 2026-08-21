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

Clusters <- readRDS(paste0(mouse_data_path, "GBM_GSE235913/GSE235676_mouse.GBM.scRNA.integrated.rds"))

Clusters_WT <- subset(Clusters, genotype == "WT")

Clusters_WT$celltype_broad <- ifelse(Clusters_WT$anno_ident %in% c("pDCs", "cDC1", "cDC2", "Macrophages", "Microglial", "Monocytes"), "TAM", Clusters_WT$anno_ident)
Clusters_WT$tissue <- "GBM"
Clusters_WT$tumor <- "Tumor"
Clusters_WT$ID <- "Tumor_GBM"

qsave(Clusters_WT, paste0(path, "Clusters.qs"))
