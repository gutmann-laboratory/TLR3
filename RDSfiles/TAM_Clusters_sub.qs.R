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

all_SO <- qread(paste0(path, "all_SO.qs"))

TAM <- subset(all_SO, celltype_broad1 == "TAM")
TAM <- DietSeurat(TAM)
TAM_list <- SplitObject(TAM, split.by = "Pathology")

TAM_list <- lapply(TAM_list, function(x){
    x <- FindVariableFeatures(x)
})

TAM_merge <- Reduce(merge, TAM_list)
VariableFeatures(TAM_merge) <- SelectIntegrationFeatures(TAM_list)

TAM_PCA <- ScaleData(TAM_merge) %>% RunPCA()

TAM_Clusters <- RunHarmony(TAM_PCA, group.by.vars = "Pathology", lambda = 0.01) %>%
RunUMAP(reduction = "harmony", min.dist = 0.8, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20L, dims = 1:10) %>%
FindNeighbors(reduction = "harmony", k.param = 20, dims = 1:10) %>% FindClusters(resolution = c(0.3, 0.5, 0.8, 1.0))

TAM_Clusters$seurat_clusters <- TAM_Clusters$RNA_snn_res.0.3

TAM_Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))

# remove 7, 9, 11, 13, 14 (nFeature < 1000), 12 (cellcycling)
TAM_Clusters_sub <- subset(TAM_Clusters, !seurat_clusters %in% c(7, 9, 11, 13, 14, 12))

TAM_Clusters_sub$celltype_new <- ifelse(TAM_Clusters_sub$seurat_clusters %in% c(6), "Monocytes",
                                       ifelse(TAM_Clusters_sub$seurat_clusters %in% c(5), "DCs",
                                             ifelse(TAM_Clusters_sub$seurat_clusters %in% c(3, 4), "Microglia",
                                                   ifelse(TAM_Clusters_sub$seurat_clusters %in% c(0, 1, 10 ,15), "Microglia",
                                                         ifelse(TAM_Clusters_sub$seurat_clusters %in% c(2), "Macrophages",
                                                               ifelse(TAM_Clusters_sub$seurat_clusters %in% c(8), "Monocytes", "Others"))))))

TAM_Clusters_sub$celltype_combine <- factor(TAM_Clusters_sub$celltype_combine, c("Microglia", "Macrophages", "Monocytes", "DCs"))

qread(paste0(path, "TAM_Clusters_sub.qs"))








