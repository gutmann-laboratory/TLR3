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

Idents(Clusters_sub) <- "celltype_broad"
TAM <- subset(Clusters_sub, idents = c("TAM")) %>% DietSeurat()

TAM_Clusters <- FindVariableFeatures(TAM, selection.method = "vst", nfeatures = 2000) %>% 
ScaleData() %>% RunPCA() %>% RunHarmony(group.by.vars = "orig.ident", lambda = 1) %>%
RunUMAP(reduction = "harmony", min.dist = 0.2, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 30L, dims = 1:10) %>%
FindNeighbors(reduction = "harmony", k.param = 30, dims = 1:10) %>% FindClusters(resolution = c(0.3, 0.35, 0.4, 0.45, 0.5, 0.8, 1.0))

TAM_Clusters$seurat_clusters <- TAM_Clusters$RNA_snn_res.1

TAM_Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))

# remove cluster 11, 14, 15 (nFeature < 1000, doublet with T cells)
TAM_Clusters_sub <- subset(TAM_Clusters, idents = c(0:10, 12:13))

TAM_Clusters_sub$celltype_new <- ifelse(TAM_Clusters_sub$seurat_clusters %in% c(0, 5, 1, 6, 12, 9, 14), "Microglia",
                                         ifelse(TAM_Clusters_sub$seurat_clusters %in% c(2, 3, 4, 10, 8), "Macrophages",
                                                ifelse(TAM_Clusters_sub$seurat_clusters %in% c(7, 13), "DCs",
                                                       ifelse(TAM_Clusters_sub$seurat_clusters %in% c("1"), "Monocytes", NA)))))))))))

TAM_Clusters_sub <- subset(TAM_Clusters_sub, !celltype_new %in% c("DCs"))

qsave(TAM_Clusters_sub, paste0(path, "OPG_TAM_Clusters_sub.qs"))
