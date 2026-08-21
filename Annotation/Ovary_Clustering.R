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

name_list <- dir(paste0(data_path, "Ovary_GSE184880/"))  

SO_list <- lapply(name_list, function(x){
    data <- Read10X(data.dir = paste0(data_path, x))
    result <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
    result$orig.ident <- x
    return(result)
})

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

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.2

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))

Clusters$celltype_broad <- ifelse(Clusters$seurat_clusters %in% c(1, 2), "T/NK cell",
                                     ifelse(Clusters$seurat_clusters %in% c(3), "TAM",
                                           ifelse(Clusters$seurat_clusters %in% c(5), "Epithelial cell", #Tumor
                                                 ifelse(Clusters$seurat_clusters %in% c(6), "Endothelial cell",
                                                       ifelse(Clusters$seurat_clusters %in% c(8, 12), "B/Plasma cell",
                                                             ifelse(Clusters$seurat_clusters %in% c(0, 4, 9, 11), "Fibroblast",
                                                                    ifelse(Clusters$seurat_clusters %in% c(7), "Proliferating Immune cell",
                                                                   ifelse(Clusters$seurat_clusters %in% c(10), "SMC", "Others"))))))))
Clusters_sub$tissue <- "Ovary"
Clusters$tumor <- ifelse(str_detect(Clusters$orig.ident, "HGSOC"), "Tumor", "Normal")
Clusters$ID <- ifelse(Clusters$orig.ident == "HGSOC1", "Tumor_Ovary1", 
                     ifelse(Clusters$orig.ident == "HGSOC2", "Tumor_Ovary2",
                           ifelse(Clusters$orig.ident == "HGSOC3", "Tumor_Ovary3",
                                 ifelse(Clusters$orig.ident == "HGSOC4", "Tumor_Ovary4",
                                       ifelse(Clusters$orig.ident == "HGSOC5", "Tumor_Ovary5",
                                             ifelse(Clusters$orig.ident == "HGSOC6", "Tumor_Ovary6",
                                                   ifelse(Clusters$orig.ident == "HGSOC7", "Tumor_Ovary7",
                                                         ifelse(Clusters$orig.ident == "Normal1", "Normal_Ovary1",
                                                               ifelse(Clusters$orig.ident == "Normal2", "Normal_Ovary2",
                                                                     ifelse(Clusters$orig.ident == "Normal3", "Normal_Ovary3",
                                                                           ifelse(Clusters$orig.ident == "Normal4", "Normal_Ovary4",
                                                                                       ifelse(Clusters$orig.ident == "Normal5", "Normal_Ovary5",
                                                                                             "Others"))))))))))))

qsave(Clusters, paste0(path, "Clusters.qs"))
