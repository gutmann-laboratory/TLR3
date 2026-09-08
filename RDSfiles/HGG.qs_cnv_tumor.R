######################################
# Tumor identification using CNVs
######################################

library(CONICSmat)
library(Seurat)
library(qs)

Harmony_Merged <- qread(paste0(path, "HGG.qs"))

Harmony_Merged <- NormalizeData(Harmony_Merged, scale.factor = 1e6)
HGG_logCPM <- GetAssayData(Harmony_Merged)
HGG_log2CPCm <- log(expm1(HGG_logCPM) + 1, 2)
rm(HGG_logCPM)

HGG_log2CPCm <- as(HGG_log2CPCm, "CsparseMatrix")
dim(HGG_log2CPCm)

HGG_log2CPCm[1:5, 1:5]

patients <- Harmony_Merged$samples
unique(patients)

celltypes <- Harmony_Merged$celltype11.21

regions <- read.table("../chromosome_arm_positions_grch38.txt",  # Provided by CONICSmat authors
                      sep = "\t", row.names = 1, header = T)

gene_pos <- getGenePositions(rownames(HGG_log2CPCm), 
                             ensembl_version = "may2025.archive.ensembl.org")

HGG_log2CPCm <- filterMatrix(HGG_log2CPCm, gene_pos[, "hgnc_symbol"], 
                             minCells=5)

normFactor <- calcNormFactors(HGG_log2CPCm)

l <- plotAll(HGG_log2CPCm, normFactor, regions, gene_pos, "HGG")
lsub <- l[, c("7p", "7q", "10p", "10q", "22q")]

hisub <- plotHistogram(lsub2, HGG_log2CPCm, clusters = 4, 
                       zscoreThreshold = 4, patients, celltypes)

Harmony_Merged$cand_tumor <- ifelse(hisub2==1, "normal", "tumor")
qsave(Harmony_Merged, paste0(path, "HGG.qs"))
