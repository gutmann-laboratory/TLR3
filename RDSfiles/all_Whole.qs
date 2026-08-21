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

tissue <- "PA"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Whole_list <- list()
Whole_list[[1]] <- Clusters
names(Whole_list)[1] <- tissue

i = 2
tissue <- "Thyroid"
Clusters <- qread(paste0(path, tissue, "/Clusters.qs"))
Whole_list[[i]] <- subset(Clusters, tumor == "Tumor")
names(Whole_list)[i] <- tissue

i = 3
tissue <- "Lung"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Whole_list[[i]] <- subset(Clusters, tumor == "Tumor")
names(Whole_list)[i] <- tissue

i = 4
tissue <- "Pancreas"
Clusters <- qread(paste0(path, tissue, "/Clusters.qs"))
Whole_list[[i]] <- Clusters
names(Whole_list)[i] <- tissue

i = 5
tissue <- "Colon"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Whole_list[[i]] <- subset(Clusters, tumor == "Tumor")
names(Whole_list)[i] <- tissue

i = 6
tissue <- "Ovary"
Clusters <- qread(paste0(path, tissue, "/Clusters.qs"))
Whole_list[[i]] <- subset(Clusters, tumor == "Tumor")
names(Whole_list)[i] <- tissue

i = 7
tissue <- "MPNST"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Whole_list[[i]] <- Clusters
names(Whole_list)[i] <- tissue

i = 8
tissue <- "Skin"
Clusters <- qread(paste0(path, tissue, "/Clusters.qs"))
Whole_list[[i]] <- Clusters
names(Whole_list)[i] <- tissue

i = 9
tissue <- "HGG"
Clusters <- qread(paste0(path, tissue, "/HGG.qs"))
Clusters$ID <- ifelse(Clusters$tumor == "GBM", "GBM",
                     ifelse(Clusters$tumor == "DIPG", "DIPG", "Others"))
Clusters$tumor <- "Tumors"
Clusters$tissue <- "Brain"
Clusters$celltype_broad <- Clusters$celltype11.21
Whole_list[[i]] <- subset(Clusters, !ID %in% c("Others"))
names(Whole_list)[i] <- tissue

i = 10
tissue <- "GBM_Astrocytoma"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Whole_list[[i]] <- Clusters
names(Whole_list)[i] <- tissue

all_Whole <- Reduce(merge, Whole_list)

qsave(all_Whole, paste0(path, "all_Whole.qs"))
