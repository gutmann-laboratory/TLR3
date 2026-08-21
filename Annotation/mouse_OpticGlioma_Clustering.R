library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(stringr)
library(qs)
library(future)
library(cowplot)
library(RColorBrewer)
library(singleCellTK)

name_list <- dir(mouse_data_path) 

SO_list <- lapply(name_list, function(x){
    data <- Read10X(data.dir = paste0(mouse_data_path, x))
    result <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
    result$orig.ident <- x
    return(result)
})

named_list <- list()

for (i in 1:n){
    result <- RenameCells(SO_list[[i]], add.cell.id = name_list[i])
    result[["percent.mt"]] <- PercentageFeatureSet(result, pattern = "^mt-")
    result[["percent.ribo"]] <- PercentageFeatureSet(result, pattern = c('rps','rpl'))
    named_list <- append(named_list, result)
}

filtered_list <- lapply(named_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 300 & percent.mt < 5 & nCount_RNA < 30000) %>% NormalizeData() %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% ScaleData()
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
    result <- subset(x, subset = scrublet_score < 0.25)
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

Clusters  <- RunHarmony(merge_PCA, group.by.vars = c("orig.ident"), lambda = 10, dims.use = 1:20) %>%
RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.1, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20) %>% 
FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 20) %>%
FindClusters(resolution = c(0.2, 0.3, 0.5))

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.5

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(mean_nFeature_RNA = median(nFeature_RNA))

# remove Cluster 7, 8, 19 (nFeature < 850)
Clusters_sub <- subset(Clusters, idents = c(0:6, 9:18, 20:24))

Clusters_sub$celltype_broad <- ifelse(Clusters_sub$seurat_clusters %in% c("3", "4", "11", "15", "24"), "TAM",
                                 ifelse(Clusters_sub$seurat_clusters %in% c("1", "18"), "T/NK cell",
                                       ifelse(Clusters_sub$seurat_clusters %in% c("20"), "proliferating cell",
                                             #ifelse(Clusters_sub$RNA_snn_res.0.5 %in% c("15"), "NK cell",
                                                   ifelse(Clusters_sub$seurat_clusters %in% c("0", "2", "5", "10"), "Oligodendrocyte",
                                                         #ifelse(Clusters_sub$RNA_snn_res.0.5 %in% c("8"), "Tumor cell",
                                                               ifelse(Clusters_sub$seurat_clusters %in% c("14", "16"), "Neuron",
                                                                     ifelse(Clusters_sub$seurat_clusters %in% c("22"), "Astrocyte",
                                                                           ifelse(Clusters_sub$seurat_clusters %in% c("9", "17"), "B cell",
                                                                                 ifelse(Clusters_sub$seurat_clusters %in% c("13"), "Stromal cell",
                                                                                       ifelse(Clusters_sub$seurat_clusters %in% c("21"), "Neutrophil",
                                                                                             ifelse(Clusters_sub$seurat_clusters %in% c("6"), "OPC",
                                                                                                    ifelse(Clusters_sub$seurat_clusters %in% c("12"), "Tumor cell",
                                                                                                   ifelse(Clusters_sub$seurat_clusters %in% c("23"), "Epithelial cell", NA))))))))))))

qsave(Clusters_sub, paste0(path, "Clusters_sub.qs"))



























