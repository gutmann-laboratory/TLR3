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

name_list <- dir(data_path)

SO_list <- lapply(name_list, function(x){
    data <- Read10X(data.dir = paste0(data_path, x))
    result <- CreateSeuratObject(counts = data, min.cells = 0, min.features = 0)
    result$orig.ident <- x
    return(result)
})

data_list <- list()
data_list[[1]] <- read.table(paste0(data_path, "/P4_NI/GSM6047641_P4-1N_matrix.tsv.gz"), sep = "\t", header = T, row.names = 1)
data_list[[2]] <- read.table(paste0(data_path, "/P4_NS/GSM6047639_P4-2N_matrix.tsv.gz"), sep = "\t", header = T, row.names = 1)
data_list[[3]] <- read.table(paste0(data_path, "/P4_TI/GSM6047640_P4-1T_matrix.tsv"), sep = "\t", header = T, row.names = 1)
data_list[[4]] <- read.table(paste0(data_path, "/P4_TS1/GSM6047637_P4-2T1_matrix.tsv"), sep = "\t", header = T, row.names = 1)
data_list[[5]] <- read.table(paste0(data_path, "/P4_TS2/GSM6047638_P4-2T2_matrix.tsv"), sep = "\t", header = T, row.names = 1)

SO_list_2 <- lapply(data_list, function(x){
    result <- CreateSeuratObject(counts = x, min.cells = 0, min.features = 0)
    return(result)
})

SO_list_2[[1]]$orig.ident <- name_list[15]
SO_list_2[[2]]$orig.ident <- name_list[16]
SO_list_2[[3]]$orig.ident <- name_list[17]
SO_list_2[[4]]$orig.ident <- name_list[18]
SO_list_2[[5]]$orig.ident <- name_list[19]

merge_list <- c(SO_list, SO_list_2)

named_list <- list()

for (i in 1:n){
    result <- RenameCells(merge_list[[i]], add.cell.id = name_list[i])
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

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.2

Clusters@meta.data %>% group_by(seurat_clusters) %>% summarise(median_nFeature = median(nFeature_RNA))


# remove Cluster13, 14 (nFeature < 750)
Clusters_sub <- subset(Clusters, seurat_clusters %in% c(0:12, 15, 16))

Clusters_sub$celltype_broad <- ifelse(Clusters_sub$seurat_clusters %in% c(0, 2), "T/NK cell",
                                     ifelse(Clusters_sub$seurat_clusters %in% c(3, 4), "TAM",
                                           ifelse(Clusters_sub$seurat_clusters %in% c(1, 7, 9, 15, 16), "Epithelial cell", #Tumor
                                                 ifelse(Clusters_sub$seurat_clusters %in% c(8), "Endothelial cell",
                                                       ifelse(Clusters_sub$seurat_clusters %in% c(5, 11), "B/Plasma cell",
                                                             ifelse(Clusters_sub$seurat_clusters %in% c(12), "Fibroblast",
                                                                    ifelse(Clusters_sub$seurat_clusters %in% c(10), "Proliferating Immune cell",
                                                                   ifelse(Clusters_sub$seurat_clusters %in% c(6), "Mast cell", "Others"))))))))

Clusters_sub$tissue <- "Lung"
Clusters_sub$tumor <- ifelse(str_detect(Clusters_sub$orig.ident, "TI|TM|TS"), "Tumor", "Normal")

Clusters_sub$ID <- ifelse(Clusters_sub$orig.ident == "P1_TI", "Tumor_Lung1", 
                     ifelse(Clusters_sub$orig.ident == "P1_TM", "Tumor_Lung2",
                           ifelse(Clusters_sub$orig.ident == "P2_TI", "Tumor_Lung3",
                                 ifelse(Clusters_sub$orig.ident == "P2_TM", "Tumor_Lung4",
                                       ifelse(Clusters_sub$orig.ident == "P3_TI1", "Tumor_Lung5",
                                             ifelse(Clusters_sub$orig.ident == "P3_TI2", "Tumor_Lung6",
                                                   ifelse(Clusters_sub$orig.ident == "P3_TM1", "Tumor_Lung7",
                                                         ifelse(Clusters_sub$orig.ident == "P3_TM2", "Tumor_Lung8",
                                                               ifelse(Clusters_sub$orig.ident == "P4_TI", "Tumor_Lung9",
                                                                     ifelse(Clusters_sub$orig.ident == "P4_TS1", "Tumor_Lung10",
                                                                           ifelse(Clusters_sub$orig.ident == "P4_TS2", "Tumor_Lung11",
                                                                                       ifelse(Clusters_sub$orig.ident == "P1_NI", "Normal_Lung1",
                                                                                             ifelse(Clusters_sub$orig.ident == "P1_NM", "Normal_Lung2",
                                                                                                   ifelse(Clusters_sub$orig.ident == "P2_NI", "Normal_Lung3",
                                                                                                          ifelse(Clusters_sub$orig.ident == "P2_NM", "Normal_Lung4",
                                                                                                                 ifelse(Clusters_sub$orig.ident == "P3_NI", "Normal_Lung5",
                                                                                                                        ifelse(Clusters_sub$orig.ident == "P3_NM", "Normal_Lung6",
                                                                                                                               ifelse(Clusters_sub$orig.ident == "P4_NI", "Normal_Lung7",
                                                                                                                                      ifelse(Clusters_sub$orig.ident == "P4_NS", "Normal_Lung8","Others")))))))))))))))))))
qsave(Clusters_sub, paste0(path, "Clusters_sub.qs"))
