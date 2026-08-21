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

SO_list <- Read10X_h5(paste0(mouse_data_path, "Pancreas_GSE203016/GSM6149735_1998-RM-1_filtered_feature_bc_matrix.h5")) %>%
    CreateSeuratObject(min.cells = 0, min.features = 0)

SO_list[["percent.mt"]] <- PercentageFeatureSet(SO_list, pattern = "^mt-")
SO_list[["percent.ribo"]] <- PercentageFeatureSet(SO_list, pattern = c('Rps','Rpl'))

filtered_list <- subset(SO_list, subset = nFeature_RNA > 500 & percent.mt < 10 & nCount_RNA < 25000) %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>%
    ScaleData()

sce_list <- as.SingleCellExperiment(filtered_list)

cell_number <- ncol(SO_list)
doublet_rate_list <- 0.008*cell_number/1000

scrublet_list <- runScrublet(sce_list, expectedDoubletRate = doublet_rate_list)

scrublet_SO_list <- as.Seurat(scrublet_list)

normalized_list <- subset(scrublet_SO_list, subset = scrublet_score < 0.30) %>%
    NormalizeData() %>%
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

merge_PCA <- ScaleData(normalized_list, verbose = FALSE) %>%
    RunPCA(npcs = 50, verbose = FALSE)

Clusters  <- RunUMAP(merge_PCA, reduction = "pca", dims = 1:20, min.dist = 0.2, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20) %>% 
FindNeighbors(reduction = "pca", dims = 1:20, k.param = 20) %>%
FindClusters(resolution = c(0.1, 0.2, 0.3, 0.5, 0.8))

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.3

Clusters$celltype_broad <- ifelse(Clusters$seurat_clusters %in% c(3), "T/NK cell",
                                     ifelse(Clusters$seurat_clusters %in% c(4), "TAM",
                                            ifelse(Clusters$seurat_clusters %in% c(8), "Mesothelial",
                                                       ifelse(Clusters$seurat_clusters %in% c(9), "Neutrophil",
                                                             ifelse(Clusters$seurat_clusters %in% c(0, 1, 2, 6), "Fibroblast",
                                                                    ifelse(Clusters$seurat_clusters %in% c(7), "Endothelial cell", 
                                                                           ifelse(Clusters$seurat_clusters %in% c(10), "Mast cell", "Acinar")))))))

Clusters$tissue <- "Pancreas"
Clusters$tumor <- "Tumor"
Clusters$ID <- "Tumor_Pancreas1"

qsave(Clusters, paste0(path, "Clusters.qs"))



                 
