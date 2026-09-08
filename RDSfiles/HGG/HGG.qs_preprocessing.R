library(Seurat)
library(DoubletFinder)
library(dplyr)
library(tidyr)
library(tibble)
library(stringr)
library(readr)
library(Matrix)

options(future.globals.maxSize = 1e10)

# Helper functions
source("run_doubletfinder_custom.R")

######
# GBM1
######
gbm1_data <- ReadMtx(mtx = paste0(path, "GBM1/E-MTAB-9435.aggregated_filtered_counts.mtx"),
                      features = paste0(path, "GBM1/E-MTAB-9435.aggregated_filtered_counts.mtx_rows"),
                      cells = paste0(path, "GBM1/E-MTAB-9435.aggregated_filtered_counts.mtx_cols"))
dim(gbm1_data)

gbm1 <- CreateSeuratObject(gbm1_data, project = "GBM1")
gbm1$samples <- sub("-.*", "", colnames(gbm1))
gbm1 <- gbm1 %>% subset(samples != "SAMEA10160202")
gbm1$samples <- factor(gbm1$samples, levels = c("SAMEA10160203",
                                                 "SAMEA10160204",
                                                 "SAMEA10160205"))

######
# GBM2
######
SF11644_data <- Read10X(paste0(path, "GBM2/SF11644"), gene.column = 1)
SF11644 <- CreateSeuratObject(SF11644_data)

SF11956_data <- Read10X(paste0(path, "GBM2/SF11956"), gene.column = 1)
SF11956 <- CreateSeuratObject(SF11956_data)

SF11977_data <- Read10X(paste0(path, "GBM2/SF11977"), gene.column = 1)
SF11977 <- CreateSeuratObject(SF11977_data)

SF11979_data <- Read10X(paste0(path, "GBM2/SF11979"), gene.column = 1)
SF11979 <- CreateSeuratObject(SF11979_data)

gbm2 <- merge(SF11644, y = c(SF11956, SF11977, SF11979),
              add.cell.ids = c("SF11644", "SF11956", "SF11977", "SF11979"),
              project = "GBM2")

gbm2$samples <- gsub("_.*", "", colnames(gbm2))
gbm2$samples <- factor(gbm2$samples, levels = c("SF11644", "SF11956",
                                                 "SF11977", "SF11979"))

gbm2[["RNA"]] <- JoinLayers(gbm2[["RNA"]])
gbm2_data <- GetAssayData(gbm2)

######
# Merge GBM1/GBM2
######
gbm_merged <- merge(gbm1, y = gbm2, project = "GBM")
gbm_merged

gbm_merged[["RNA"]] <- JoinLayers(gbm_merged[["RNA"]])

######
# DIPG
######
dipg_qumi <- read.delim(paste0(path, "dipg_qumi.txt")) # Output from DIPG_to_qumi.R quasi-UMI normalization

rownames(dipg_qumi) <- dipg_qumi$Gene
dipg_qumi <- dipg_qumi[, colnames(dipg_qumi) != "Gene"]

dipg <- CreateSeuratObject(dipg_qumi, project = "DIPG")
dipg$samples <- sub("\\..*", "", colnames(dipg))
dipg[["RNA"]]@misc$tpm <- dipg_tpm

########
# Merge GBM/DIPG
########
Merged <- merge(dipg, y = gbm_merged, add.cell.ids = c("DIPG", "GBM"),
                project = "HGG")
head(colnames(Merged))
Idents(Merged) <- "samples"
table(Idents(Merged))

Merged$tumor <- sub("_.*", "", colnames(Merged))
table(Merged$tumor)
study <- c(rep("DIPG", 6), rep("GBM1", 3), rep("GBM2", 4))
names(study) <- levels(Merged)
Merged <- RenameIdents(Merged, study)
Merged$study <- Idents(Merged)
table(Idents(Merged))

Merged[["RNA"]] <- JoinLayers(Merged[["RNA"]])

######################################
# QC
######################################
Idents(Merged) <- "samples"
Merged[["percent.mt"]] <- PercentageFeatureSet(Merged, pattern = "^MT-")
Merged[["percent.rib"]] <- PercentageFeatureSet(Merged, pattern = "RP[SL]")
VlnPlot(Merged, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"), ncol = 4)
median(Merged$nFeature_RNA)

QC_Merged <- subset(Merged, percent.mt < 10)

######################################
# DoubletFinder (run per sample)
######################################
# run_doubletfinder_custom() is defined in run_doubletfinder_custom.R
samp_split <- SplitObject(QC_Merged, split.by = "samples")
samp_split <- lapply(samp_split, run_doubletfinder_custom)

sglt_dblt_metadata <- data.frame(bind_rows(samp_split))
rownames(sglt_dblt_metadata) <- sglt_dblt_metadata$row_names
sglt_dblt_metadata$row_names <- NULL
head(sglt_dblt_metadata)
QC_Merged <- AddMetaData(QC_Merged, sglt_dblt_metadata, col.name = "doublet_finder")

# Compare QC metrics between doublets and singlets per sample
VlnPlot(QC_Merged, group.by = "samples", split.by = "doublet_finder",
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rib"),
        ncol = 3, pt.size = 0) + theme(legend.position = "right")

# Remove doublets
DF_Merged <- subset(QC_Merged, doublet_finder == "Singlet")
qsave(DF_Merged, paste0(path, "DF_Merged.qs")) # Input for HGG.qs.R
