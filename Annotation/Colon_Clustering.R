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

data <- read.csv(paste0(data_path, "/Colon_GSE200997/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv"), header = TRUE, row.names = 1)
metadata <- read.csv(paste0(data_path, "/Colon_GSE200997/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv"), header = TRUE, row.names = 1)

SO <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0, meta.data = metadata)

SO_list <- SplitObject(SO, split.by = "samples")

named_list <- list()

for (i in 1:23){
    result <- SO_list[[i]]
    result[["percent.mt"]] <- PercentageFeatureSet(result, pattern = "^MT-")
    result[["percent.ribo"]] <- PercentageFeatureSet(result, pattern = c('RPS','RPL'))
    named_list <- append(named_list, result)
}

filtered_list <- lapply(named_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 10 & nCount_RNA < 50000) %>% NormalizeData() %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% ScaleData()
    return(x)
}) 

sce_list <- lapply(filtered_list, function(x){
    x <- as.SingleCellExperiment(x)
    return(x)
})


doublet_rate_list <- c()

doublet_rate_list <- lapply(named_list, function(x){
    cell_number <- ncol(x)
    result <- 0.008*cell_number/1000
    doublet_rate_list <- c(doublet_rate_list, result)
})

scrublet_list <- mapply(function(x, y){
    result <- runScrublet(x, expectedDoubletRate = y)
    return(result)
}, sce_list, doublet_rate_list)

scrublet_SO_list <- lapply(scrublet_list, function(x){
    result <- as.Seurat(x)
    return(result)
})

scrublet_filter_list <- lapply(scrublet_SO_list, function(x){
    result <- subset(x, subset = scrublet_score < 0.20)
    return(result)
})

normalized_list <- lapply(scrublet_filter_list, function(x){
    result <- NormalizeData(x) %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
    return(result)
})

merge <- Reduce(merge, normalized_list)

VariableFeatures(merge) <- SelectIntegrationFeatures(object.list = normalized_list, nfeatures = 2000)

merge_scale <- ScaleData(merge, verbose = FALSE)

merge_PCA <- RunPCA(merge_scale, npcs = 50, verbose = FALSE)

Clusters  <- RunHarmony(merge_PCA, group.by.vars = c("samples"), lambda = 10, dims.use = 1:20) %>%
RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.2, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20) %>%
FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 20) %>%
FindClusters(resolution = c(0.1, 0.2, 0.3, 0.5, 0.8)) 

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.3

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))


# remove Cluster13 (nFeature < 750)
Clusters_sub <- subset(Clusters, seurat_clusters %in% c(0:12, 14))

Clusters_sub$celltype_broad <- ifelse(Clusters_sub$seurat_clusters %in% c(0, 1, 3, 4), "T/NK cell",
                                     ifelse(Clusters_sub$seurat_clusters %in% c(7), "TAM",
                                           ifelse(Clusters_sub$seurat_clusters %in% c(5, 12), "Epithelial cell", #Tumor
                                                 ifelse(Clusters_sub$seurat_clusters %in% c(9), "Endothelial cell",
                                                       ifelse(Clusters_sub$seurat_clusters %in% c(2, 6, 11), "B/Plasma cell",
                                                             ifelse(Clusters_sub$seurat_clusters %in% c(10), "Fibroblast",
                                                                    ifelse(Clusters_sub$seurat_clusters %in% c(8), "Proliferating Immune cell", 
                                                                           ifelse(Clusters_sub$seurat_clusters %in% c(14), "Mast cell", "Others"))))))))

Clusters_sub$tissue <- "Colon"
Clusters_sub$tumor <- ifelse(Clusters_sub$Condition == "Normal", "Normal", "Tumor")

Clusters_sub$ID <- ifelse(Clusters_sub$samples == "T_cac1", "Tumor_Colon1",
                          ifelse(Clusters_sub$samples == "T_cac2", "Tumor_Colon2",
                                 ifelse(Clusters_sub$samples == "T_cac3", "Tumor_Colon3",
                                        ifelse(Clusters_sub$samples == "T_cac4", "Tumor_Colon4",
                                               ifelse(Clusters_sub$samples == "T_cac5", "Tumor_Colon5",
                                                      ifelse(Clusters_sub$samples == "T_cac6", "Tumor_Colon6",
                                                             ifelse(Clusters_sub$samples == "T_cac7", "Tumor_Colon7",
                                                                    ifelse(Clusters_sub$samples == "T_cac8", "Tumor_Colon8",
                                                                           ifelse(Clusters_sub$samples == "T_cac9", "Tumor_Colon9",
                                                                                  ifelse(Clusters_sub$samples == "T_cac10", "Tumor_Colon10",
                                                                                         ifelse(Clusters_sub$samples == "T_cac11", "Tumor_Colon11",
                                                                                                ifelse(Clusters_sub$samples == "T_cac12", "Tumor_Colon12",
                                                                                                       ifelse(Clusters_sub$samples == "T_cac13", "Tumor_Colon13",
                                                                                                              ifelse(Clusters_sub$samples == "T_cac14", "Tumor_Colon14",
                                                                                                                     ifelse(Clusters_sub$samples == "T_cac15", "Tumor_Colon15",
                                                                                                                            ifelse(Clusters_sub$samples == "T_cac16", "Tumor_Colon16",
                                                                                                                                   ifelse(Clusters_sub$samples == "B_cac4", "Normal_Colon1",
                                                                                                                                          ifelse(Clusters_sub$samples == "B_cac6", "Normal_Colon2",
                                                                                                                                                 ifelse(Clusters_sub$samples == "B_cac7", "Normal_Colon3",
                                                                                                                                                        ifelse(Clusters_sub$samples == "B_cac10", "Normal_Colon4",
                                                                                                                                                               ifelse(Clusters_sub$samples == "B_cac11", "Normal_Colon5",
                                                                                                                                                                      ifelse(Clusters_sub$samples == "B_cac14", "Normal_Colon6",
                                                                                                                                                                             ifelse(Clusters_sub$samples == "B_cac15", "Normal_Colon7",
 "Others")))))))))))))))))))))))

qsave(Clusters_sub, paste0(path, "Clusters_sub.qs"))
