library(Seurat)
library(harmony)
library(dplyr)
library(clustree)
library(cowplot)
library(RColorBrewer)
library(ggrepel)
library(qs)

DF_Merged <- qread(paste0(path, "DF_Merged.qs")) # Output from HGG.qs_preprocessing.R

DF_Merged[["RNA"]] <- JoinLayers(DF_Merged[["RNA"]])
DF_Merged[["RNA"]] <- split(DF_Merged[["RNA"]], f = DF_Merged$samples)

# log2 normalization
L2N_Merged <- DF_Merged %>% NormalizeData() %>% FindVariableFeatures(selection.method = "vst", nfeatures = 2000)
L2N_Merged[["RNA"]] <- JoinLayers(L2N_Merged[["RNA"]])

L2N_Merged <- L2N_Merged %>% ScaleData(verbose = FALSE) %>% RunPCA(npcs = 50, verbose = FALSE)

# Integration with Harmony
Harmony_Merged <- RunHarmony(L2N_Merged, group.by.vars = c("samples", "tumor"),
                              lambda = 10, dims.use = 1:20) %>%
  RunUMAP(reduction = "harmony", dims = 1:20, min.dist = 0.1,
          metric = "euclidean", n.neighbors = 20) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, k.param = 20) %>%
  FindClusters(resolution = c(0.2, 0.3, 0.5, 0.6, 0.7, 0.8))

clustree(Harmony_Merged, prefix = "RNA_snn_res.",
         node_colour = "PTPRC", node_colour_aggr = "mean")
clustree(Harmony_Merged, prefix = "RNA_snn_res.",
         node_colour = "AIF1", node_colour_aggr = "mean")

Harmony_Merged$seurat_clusters <- Harmony_Merged$RNA_snn_res.0.3
length(levels(Harmony_Merged$RNA_snn_res.0.3))
Harmony_Merged$seurat_clusters <-
  factor(Harmony_Merged$seurat_clusters,
         c(0:length(levels(Harmony_Merged$RNA_snn_res.0.3))))

qsave(Harmony_Merged, "HGG_no_cnv_tumor.qs") # Input for HGG.qs_cnv_tumor.R


Harmony_Merged_tumor <- qread("HGG_with_cnv_tumor.qs") # Output from HGG.qs_cnv_tumor.R

Idents(Harmony_Merged_tumor) <- "seurat_clusters"
table(Idents(Harmony_Merged_tumor), Harmony_Merged_tumor$cand_tumor) # Determine consensus celltype for each cluster
new.cluster.ids <- c("Tumor1", "TAM", "Tumor2", "Tumor3", "Tumor1", "Tumor2",
                      "TAM", "Tumor1", "Tumor2", "Tumor2", "TAM", "TAM", "TAM",
                      "Oligo", "EC", "Unknown")
names(new.cluster.ids) <- levels(Harmony_Merged_tumor)
Harmony_Merged_tumor <- RenameIdents(Harmony_Merged_tumor, new.cluster.ids)
Idents(Harmony_Merged_tumor) <- factor(Idents(Harmony_Merged_tumor),
                                  levels = c("Tumor1", "Tumor2", "Tumor3",
                                             "TAM", "Oligo", "EC", "Unknown"))
Harmony_Merged_tumor$celltype11.21 <- Idents(Harmony_Merged_tumor)

qsave(Harmony_Merged_tumor, paste0(path, "HGG.qs"))
