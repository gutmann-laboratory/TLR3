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

Clusters <- readRDS(paste0(mouse_data_path, "Pancreas_GSE244992/GSE244992_08_processed_seurat.rds"))

Clusters$celltype_broad <- ifelse(Clusters$cell_type %in% c("B"), "B/Plasma cell",
                                  ifelse(Clusters$cell_type %in% c("cancer"), "Tumor",
                                        ifelse(Clusters$cell_type %in% c("cycling_CD8_T", "cytotoxic_CD8_T", "effector_CD4_T", "exhausted_CD8_T", "naive_CD8_T", "natural_killer", "regulatory_CD4_T"), "T/NK cell",
                                              ifelse(Clusters$cell_type %in% c("dendritic_cell", "monocyte", "myeloid_derived_supressor", "tumor_associated_macrophage"), "TAM",
                                                    ifelse(Clusters$cell_type %in% c("endothelial"), "Endothelial cell",
                                                          "Others")))))

Clusters$tissue <- "Pancreas"
Clusters$tumor <- "Tumor"
Clusters$ID <- "Tumor_Pancreas2"

qsave(Clusters, paste0(path, "Clusters.qs"))
