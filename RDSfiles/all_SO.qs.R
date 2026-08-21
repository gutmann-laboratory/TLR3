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
library(tidyverse)
library(tidyr)

tissue <- "PA"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))
Clusters$ID <- "PA"

TAM_list <- list()
Olig_list <- list()
Tumor_list <- list()

TAM_list[[1]] <- subset(Clusters, celltype_broad == "TAM")
names(TAM_list)[1] <- tissue
TAM_list[[1]]

Tumor_list[[1]] <- subset(Clusters, celltype_broad == "Tumor")
names(Tumor_list)[1] <- tissue
Tumor_list[[1]]

tissue <- "HGG"
Clusters <- qread(paste0(path, tissue, "/HGG.qs"))

Clusters$ID <- ifelse(Clusters$tumor == "GBM", "GBM",
                     ifelse(Clusters$tumor == "DIPG", "DIPG", "Others"))

Clusters$celltype_broad <- Clusters$celltype11.21

i = 2
TAM_list[[2]] <- subset(Clusters, celltype11.21 %in% c("TAM"))
names(TAM_list)[2] <- tissue
TAM_list[[2]]

Tumor_list[[2]] <- subset(Clusters, celltype11.21 %in% c("Tumor1", "Tumor2", "Tumor3"))
names(Tumor_list)[2] <- tissue
Tumor_list[[2]]

Olig_list[[2]] <- subset(Clusters, celltype11.21 %in% c("Oligo"))
names(Olig_list)[2] <- tissue
Olig_list[[2]]

tissue <- "GBM_Astrocytoma"
Clusters <- qread(paste0(path, tissue, "/Clusters_sub.qs"))

Clusters$ID <- Clusters$orig.ident

TAM_list[[3]] <- subset(Clusters, celltype_broad %in% c("TAM"))
names(TAM_list)[3] <- tissue
TAM_list[[3]]

Tumor_list[[3]] <- subset(Clusters, celltype_broad %in% c("Tumor"))
names(Tumor_list)[3] <- tissue
Tumor_list[[3]]

Olig_list[[3]] <- subset(Clusters, celltype_broad %in% c("Oligodendrocyte"))
names(Olig_list)[3] <- tissue
Olig_list[[3]]

all_TAM <- Reduce(merge, TAM_list)
all_Oligo <- merge(Olig_list[[2]], Olig_list[[3]])
all_Tumor <- Reduce(merge, Tumor_list)

all_Tumor$celltype_broad <- "Tumor"

all_SO <- merge(x = all_TAM, y = c(all_Oligo, all_Tumor))

all_SO$Pathology <- ifelse(all_SO$ID == "Astrocytoma", "Grade II Astrocytoma", as.character(all_SO$ID))
all_SO$Pathology <- factor(all_SO$Pathology, c("PA", "DIPG", "Grade II Astrocytoma", "GBM"))

all_SO$celltype_broad1 <- ifelse(all_SO$celltype_broad == "Oligo", "Oligodendrocyte", all_SO$celltype_broad)

all_SO$Pathology_celltype <- paste(all_SO$Pathology, all_SO$celltype_broad1, sep = "_") #group.by and split.by
Idents(all_SO) <- "Pathology_celltype"

qsave(path, "all_SO.qs")
