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

data_list <- dir("PATH/Thyroid_GSE232237")

SO_list <- lapply(name_list, function(x){
    data <- read.table(paste0("PATH/Thyroid_GSE232237/", data_list[str_detect(data_list, x)]), sep = "\t", header = T, row.names = 1)
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


data <- read.table("PATH/Thyroid_N_GSE182416/GSE182416_Thyroid_normal_7samples_54726cells_raw_count.txt.gz", sep = "\t", header = T, row.names = 1)
metadata <- read.table("PATH/Thyroid_N_GSE182416/GSE182416_Thyroid_normal_7samples_metadata.txt.gz", sep = "\t", header = T, row.names = 1)
normal_SO <- CreateSeuratObject(counts = select(data, -Index), min.cells = 0, min.features = 0, meta.data = metadata)
normal_SO_sub <- subset(normal_SO, orig.ident %in% c("Thy05", "Thy06", "N3-GEX"))
normal_SO_list <- SplitObject(normal_SO_sub, split.by = "orig.ident")

normal_SO_list_new <- list()

for (i in 1:3){
    result <- normal_SO_list[[i]]
    result[["percent.mt"]] <- PercentageFeatureSet(result, pattern = "^MT-")
    result[["percent.ribo"]] <- PercentageFeatureSet(result, pattern = c('RPS','RPL'))
    normal_SO_list_new <- append(normal_SO_list_new, result)
}


merge_list <- c(named_list, normal_SO_list_new)

filtered_list <- lapply(merge_list, function(x){
    x <- subset(x, subset = nFeature_RNA > 500 & percent.mt < 10 & nCount_RNA < 50000) %>% NormalizeData() %>% FindVariableFeatures(selection.method = 'vst', nfeatures = 2000) %>% ScaleData()
    return(x)
}) 

sce_list <- lapply(filtered_list, function(x){
    x <- as.SingleCellExperiment(x)
    return(x)
})


doublet_rate_list <- c()

doublet_rate_list <- lapply(merge_list, function(x){
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

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.15

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))

Clusters$celltype_broad <- ifelse(Clusters$seurat_clusters %in% c(0, 1), "T/NK cell",
                                     ifelse(Clusters$seurat_clusters %in% c(2), "TAM",
                                           ifelse(Clusters$seurat_clusters %in% c(3, 11), "Thyroid cell",
                                                 ifelse(Clusters$seurat_clusters %in% c(4), "Endothelial cell",
                                                       ifelse(Clusters$seurat_clusters %in% c(5, 8, 9), "B/Plasma cell",
                                                             ifelse(Clusters$seurat_clusters %in% c(6, 7), "Fibroblast",
                                                                   ifelse(Clusters$seurat_clusters %in% c(10), "Mast cell", "Others")))))))

Clusters_sub$tissue <- "Thyroid"
Clusters$tumor <- ifelse(str_detect(Clusters$orig.ident, "AT|PT"), "Tumor", "Normal")

Clusters$ID <- ifelse(Clusters$orig.ident == "PT3", "Tumor_Thyroid1", 
                     ifelse(Clusters$orig.ident == "PT5", "Tumor_Thyroid2",
                           ifelse(Clusters$orig.ident == "PT7", "Tumor_Thyroid3",
                                 ifelse(Clusters$orig.ident == "PT8", "Tumor_Thyroid4",
                                       ifelse(Clusters$orig.ident == "PT9", "Tumor_Thyroid5",
                                             ifelse(Clusters$orig.ident == "PT10", "Tumor_Thyroid6",
                                                   ifelse(Clusters$orig.ident == "PT12", "Tumor_Thyroid7",
                                                         ifelse(Clusters$orig.ident == "AT9", "Tumor_Thyroid8",
                                                               ifelse(Clusters$orig.ident == "AT13", "Tumor_Thyroid9",
                                                                     ifelse(Clusters$orig.ident == "AT16", "Tumor_Thyroid10",
                                                                           ifelse(Clusters$orig.ident == "AT17", "Tumor_Thyroid11",
                                                                                 ifelse(Clusters$orig.ident == "AT20", "Tumor_Thyroid12",
                                                                                       ifelse(Clusters$orig.ident == "Thy05", "Normal_Thyroid1",
                                                                                             ifelse(Clusters$orig.ident == "Thy06", "Normal_Thyroid2",
                                                                                                   ifelse(Clusters$orig.ident == "N3-GEX", "Normal_Thyroid3", "Others")))))))))))))))

qsave(Clusters, paste0(path, "Clusters.qs"))
