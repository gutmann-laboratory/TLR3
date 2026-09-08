######################################
# Normalization of DIPG counts to quasi-UMIs
######################################

library(dplyr)
library(Matrix)
library(parallel)
library(Seurat)

## Obtaining and preprocessing read count data without UMIs
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(EDASeq))
library(scRNAseq)
require(quminorm)

dipg_tpm <- read.delim(paste0(path, "GSE102130_K27Mproject.RSEM.vh20170621.txt"))
rownames(dipg_tpm) <- dipg_tpm$Gene
dipg_tpm <- dipg_tpm[, colnames(dipg_tpm) != "Gene"]

dipg_seurat <- CreateSeuratObject(dipg_tpm, assay = "TPM", project = "DIPG")
dipg_sce <- as.SingleCellExperiment(dipg_seurat)
dipg_tpm <- GetAssayData(dipg_seurat)
assay(dipg_sce, "tpm") <- dipg_tpm
assay(dipg_sce, "counts") <- NULL

gbm_merged <- readRDS(paste0(path, "gbm_merged.RDS"))
m <- round(GetAssayData(gbm_merged))

## Fitting a Poisson-lognormal distribution to UMI count data
m<-m[rowSums(m)>0,]
keep<-which(colSums(m)>25500) # Keep top 1000
mkeep <- m[,keep]

ncores <- max(1, detectCores()-1)
chunks <- split(seq_len(ncol(mkeep)), cut(seq_len(ncol(mkeep)), ncores, labels = F))

cl <- makeCluster(ncores)
clusterExport(cl, varlist = c("mkeep", "chunks"))
clusterEvalQ(cl, {library(quminorm)})
fit_list <- parLapply(cl, chunks, function(cols) {
  poilog_mle(mkeep[,cols])
})

stopCluster(cl)

fit <- do.call(rbind, lapply(fit_list, function(x) x$sig))
summary(as.vector(fit)) # shape=2.2

plot(table(m[,1]),main="GBM UMI counts")
hist(log1p(m[,1]),main="GBM log(1+UMI counts)")


## Quantile normalization to obtain quasi-UMI counts
system.time(dipg_sce<-quminorm(dipg_sce,assayName="tpm",shape=2.2))
qumi<-assay(dipg_sce,"qumi_poilog_2.2")
plot(table(qumi[,1]),main="DIPG QUMI counts")
hist(log1p(qumi[,1]),main="DIPG log(1+QUMI counts)")

qumi.df <- as.data.frame(qumi)
qumi.df$Gene <- rownames(qumi.df)
write_delim(qumi.df, "dipg_qumi.txt", delim = "\t")

