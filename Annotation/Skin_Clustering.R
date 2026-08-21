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

name_list <- dir(paste0(data_path, "Skin_GSE215120/")

data_list <- lapply(name_list, function(x){
    x <- dir(paste0(data_path, "Skin_GSE215120/", x))
})
                 
SO_list <- mapply(function(x, y){
    data <- Read10X_h5(file = paste0(data_path, "Skin_GSE215120/", x, "/", y))
    result <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
    result$orig.ident <- x
    return(result)
}, name_list, data_list)

named_list <- list()

for (i in 1:n){
    result <- RenameCells(SO_list[[i]], add.cell.id = name_list[i])
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

Clusters  <- RunHarmony(merge_PCA, group.by.vars = c("orig.ident"), lambda = 10, dims.use = 1:20) %>%
RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.2, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20) %>%
FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 20) %>%
FindClusters(resolution = c(0.1, 0.2, 0.3, 0.5, 0.8)) 

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.8

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))


Clusters$celltype_broad <- ifelse(Clusters$seurat_clusters %in% c(1, 4, 15, 18), "T/NK cell",
                                     ifelse(Clusters$seurat_clusters %in% c(17, 19), "TAM",
                                           ifelse(Clusters$seurat_clusters %in% c(0, 2, 3, 5, 6, 8, 9, 10, 12, 16), "Tumor", #Tumor
                                                 ifelse(Clusters$seurat_clusters %in% c(7), "Endothelial cell",
                                                       ifelse(Clusters$seurat_clusters %in% c(14), "B/Plasma cell",
                                                                    ifelse(Clusters$seurat_clusters %in% c(11, 13), "Fibroblast","Others"))))))

Clusters$tissue <- "Skin"
Clusters$tumor <- "Tumor"
Clusters$ID <- ifelse(Clusters$orig.ident == "Acral1", "Tumor_Skin1",
                      ifelse(Clusters$orig.ident == "Acral2", "Tumor_Skin2",
                             ifelse(Clusters$orig.ident == "Acral3", "Tumor_Skin3",
                                    ifelse(Clusters$orig.ident == "Acral4", "Tumor_Skin4",
                                           ifelse(Clusters$orig.ident == "Acral5", "Tumor_Skin5",
                                                  ifelse(Clusters$orig.ident == "Acral6", "Tumor_Skin6",
                                                         ifelse(Clusters$orig.ident == "Cutaneous1", "Tumor_Skin7",
                                                                ifelse(Clusters$orig.ident == "Cutaneous2", "Tumor_Skin8",
                                                                       ifelse(Clusters$orig.ident == "Cutaneous3", "Tumor_Skin9",
                                                                                             "Others")))))))))

qsave(Clusters, paste0(path, "Clusters.qs"))
