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

name_list <- dir(paste0(mouse_data_path, "Medulloblastoma_GSE266048/"))

SO_list <- lapply(name_list, function(x){
    data <- Read10X(data.dir = paste0(mouse_data_path, "Medulloblastoma_GSE266048/", x))
    result <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
    result$orig.ident <- x
    return(result)
})

named_list <- list()

for (i in 1:3){
    result <- SO_list[[i]]
    result[["percent.mt"]] <- PercentageFeatureSet(result, pattern = "^mt-")
    result[["percent.ribo"]] <- PercentageFeatureSet(result, pattern = c('Rps','Rpl'))
    named_list <- append(named_list, result)
}

filtered_list <- lapply(named_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 10 & nCount_RNA < 25000) %>% NormalizeData() %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% ScaleData()
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
    result <- subset(x, subset = scrublet_score < 0.30)
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

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.3

Clusters$celltype_broad <- ifelse(Clusters$seurat_clusters %in% c(13), "T/NK cell",
                                     ifelse(Clusters$seurat_clusters %in% c(9), "TAM",
                                           ifelse(Clusters$seurat_clusters %in% c(0, 1, 4, 7, 8), "Tumor",
                                                  ifelse(Clusters$seurat_clusters %in% c(11, 12), "Astrocyte", 
                                                         ifelse(Clusters$seurat_clusters %in% c(2, 3), "Neuron", "Undefined")))))

Clusters$tissue <- "MB"
Clusters$tumor <- "Tumor"
Clusters$ID <- ifelse(Clusters$orig.ident == "Tumor_1", "Tumor1_MB",
                     ifelse(Clusters$orig.ident == "Tumor_2", "Tumor2_MB", "Tumor3_MB"))

qsave(Clusters, paste0(path, "Clusters.qs"))



                 
