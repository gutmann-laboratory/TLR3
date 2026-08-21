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

tissue <- "OPG"
Clusters <- qread(paste0(path, "OPG/TAM_Clusters_sub.qs"))
TAM_list[[1]] <- Clusters
names(TAM_list)[1] <- tissue
TAM_list[[1]]

i = 2
tissue <- "Colon"
Clusters <- qread(paste0(path, "mouse_Colon/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 3
tissue <- "GBM"
Clusters <- qread(paste0(path, "mouse_GBM/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 4
tissue <- "Medulloblastoma"
Clusters <- qread(paste0(path, "mouse_MB/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 5
tissue <- "Pancreas1"
Clusters <- qread(paste0(path, "mouse_Pancreas1/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 6
tissue <- "Pancreas2"
Clusters <- qread(paste0(path, "mouse_Pancreas2/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 7
tissue <- "Pancreas3"
Clusters <- qread(paste0(path, "mouse_Pancreas3/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

i = 8
tissue <- "Skin"
Clusters <- qread(paste0(path, "mouse_Skin/Clusters.qs"))
TAM_list[[i]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[i] <- tissue

merge <- Reduce(merge, TAM_list)
VariableFeatures(merge) <- SelectIntegrationFeatures(object.list = TAM_list, nfeatures = 2000)
merge_scale <- ScaleData(merge, verbose = FALSE)
merge_PCA <- RunPCA(merge_scale, npcs = 50, verbose = FALSE)

Clusters  <- RunHarmony(merge_PCA, group.by.vars = c("ID"), lambda = 10, dims.use = 1:20) %>%
RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.2, umap.method = 'umap-learn', metric = "euclidean", n.neighbors = 20) %>% #min.dist小さければクラスターが離れる、各クラスターが凝集する感じ
FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 20) %>%
FindClusters(resolution = c(0.2, 0.3, 0.5, 0.8, 1.0, 1.2))

Clusters$seurat_clusters <- Clusters$RNA_snn_res.0.2

# remove Cluster10 (doublet with T)
Clusters_sub <- subset(Clusters, !(seurat_clusters %in% c(10)))

Clusters_sub$celltype_new <- ifelse(Clusters_sub$seurat_clusters %in% c(1), "Microglia",
                                   ifelse(Clusters_sub$seurat_clusters %in% c(2, 5, 6, 0, 7, 11, 12), "Macrophages",
                                                 ifelse(Clusters_sub$seurat_clusters %in% c(3), "Monocytes",
                                                        ifelse(Clusters_sub$seurat_clusters %in% c(4), "Monocytes",
                                                               ifelse(Clusters_sub$seurat_clusters %in% c(9, 8), "DCs",
                                                                   "Others"))))))

Clusters_sub <- subset(Clusters_sub, !celltype_new %in% c("DCs"))

Clusters_sub$tumor_tissue <- ifelse(Clusters_sub$ID %in% c("Het_PA", "Normal_PA"), "Normal_ON",
                                                     ifelse(Clusters_sub$ID %in% c("Tumor_Colon"), "Tumor_Colon",
                                                           ifelse(Clusters_sub$ID %in% c("Tumor_GBM"), "GBM",
                                                                 ifelse(Clusters_sub$ID %in% c("Tumor_PA"), "Tumor_ON",
                                                                       ifelse(Clusters_sub$ID %in% c("Tumor_Pancreas1", "Tumor_Pancreas2", "Tumor1_Pancreas3", "Tumor2_Pancreas3"), "Tumor_Pancreas",
                                                                             ifelse(Clusters_sub$ID %in% c("Tumor1_MB", "Tumor2_MB", "Tumor3_MB"), "MB",
                                                                                   ifelse(Clusters_sub$ID %in% c("Tumor1_Skin", "Tumor2_Skin"), "Tumor_Skin",
                                                                                         ifelse(Clusters_sub$ID %in% c("WT_Colon"), "Normal_Colon", "Others"))))))))


Clusters_sub$new_label <- ifelse(Clusters_sub$tumor_tissue %in% c("Tumor_ON", "GBM", "MB"), "Brain", 
                                 ifelse(Clusters_sub$tumor_tissue %in% c("Tumor_Pancreas"), "Pancreas",
                                       ifelse(Clusters_sub$tumor_tissue %in% c("Tumor_Colon"), "Colon",
                                             ifelse(Clusters_sub$tumor_tissue %in% c("Tumor_Skin"), "Skin", "Others"))))

Clusters_tumor <- subset(Clusters_sub, new_label != "Others")

# sFig1A
features = c("Tlr3")

DotPlot(Clusters_tumor, features = features, group.by = "new_label") +
    scale_y_discrete(limits = rev(levels(Clusters_tumor$new_label))) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) +
    theme(axis.text.y = element_text(size = 20))

# sFig1B

Clusters_brain <- subset(Clusters_tumor, subset = new_label == "Brain")
Clusters_brain$HGG <- ifelse(Clusters_brain$tumor_tissue %in% c("GBM", "MB"), "HGG", "LGG")

features = c("Tlr3")

DotPlot(Clusters_brain, features = features, group.by = "HGG", dot.scale = 10, scale.min = 0) +
    #scale_y_discrete(limits = c("")) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(size = 20, angle = 0, hjust = 0.5, vjust = 0.5, face = "italic")) +
    theme(axis.text.y = element_text(size = 20))


