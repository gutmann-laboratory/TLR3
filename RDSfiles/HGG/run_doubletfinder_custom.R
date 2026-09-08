######################################
# run_doubletfinder_custom()
######################################
# Runs DoubletFinder on a single-sample Seurat object and returns a data
# frame of cell IDs with a 'Singlet'/'Doublet' classification.
#
# Args:
#   seu_sample_subset : Seurat object containing a single sample
#   multiplet_rate     : expected multiplet rate (optional). If NULL, the
#                        rate is estimated from the number of recovered
#                        cells using the 10X multiplet rate table below.
#
# Returns:
#   data.frame with columns 'row_names' (cell ID) and 'doublet_finder'
#   ('Singlet' / 'Doublet')
#
# Requires: Seurat, DoubletFinder, dplyr, tibble

run_doubletfinder_custom <- function(seu_sample_subset, multiplet_rate = NULL) {
  print(paste0("Sample ", unique(seu_sample_subset[["samples"]]), "..........."))

  if (is.null(multiplet_rate)) {
    print("multiplet_rate not provided....... estimating multiplet rate from cells in dataset")

    # 10X multiplet rate table (https://rpubs.com/kenneditodd/doublet_finder_example)
    multiplet_rates_10x <- data.frame(
      "Multiplet_rate" = c(0.004, 0.008, 0.0160, 0.023, 0.031, 0.039, 0.046, 0.054, 0.061, 0.069, 0.076),
      "Loaded_cells" = c(800, 1600, 3200, 4800, 6400, 8000, 9600, 11200, 12800, 14400, 16000),
      "Recovered_cells" = c(500, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000)
    )
    print(multiplet_rates_10x)

    if (nrow(seu_sample_subset@meta.data) < 500) {
      multiplet_rate <- nrow(seu_sample_subset@meta.data) / 125000
    } else {
      multiplet_rate <- multiplet_rates_10x %>%
        dplyr::filter(Recovered_cells < nrow(seu_sample_subset@meta.data)) %>%
        dplyr::slice(which.max(Recovered_cells)) %>% # select the min threshold for this many recovered cells
        dplyr::select(Multiplet_rate) %>% as.numeric(as.character()) # expected multiplet rate for that number of recovered cells
    }
    print(paste("Setting multiplet rate to", multiplet_rate))
  }

  # Standard preprocessing
  sample <- NormalizeData(seu_sample_subset)
  sample <- FindVariableFeatures(sample)
  sample <- ScaleData(sample)
  sample <- RunPCA(sample, nfeatures.print = 10)

  # Identify the number of informative PCs
  stdv <- sample[["pca"]]@stdev
  percent_stdv <- (stdv / sum(stdv)) * 100
  cumulative <- cumsum(percent_stdv)
  co1 <- which(cumulative > 90 & percent_stdv < 5)[1]
  co2 <- sort(which((percent_stdv[1:length(percent_stdv) - 1] -
                        percent_stdv[2:length(percent_stdv)]) > 0.1),
              decreasing = TRUE)[1] + 1
  min_pc <- min(co1, co2)

  sample <- RunUMAP(sample, dims = 1:min_pc)
  sample <- FindNeighbors(object = sample, dims = 1:min_pc)
  sample <- FindClusters(object = sample, resolution = 0.1)

  sweep_list <- paramSweep(sample, PCs = 1:min_pc, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_list)
  bcmvn <- find.pK(sweep_stats) # metric to find the optimal pK (max mean-variance normalized bimodality coefficient)
  optimal.pk <- bcmvn %>%
    dplyr::filter(BCmetric == max(BCmetric)) %>%
    dplyr::select(pK)
  optimal.pk <- as.numeric(as.character(optimal.pk[[1]]))

  # Homotypic doublet proportion estimate
  annotations <- sample@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations)

  nExp.poi <- round(multiplet_rate * nrow(sample@meta.data)) # expected number of multiplets
  nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop)) # expected number of heterotypic doublets

  sample <- doubletFinder(seu = sample,
                           PCs = 1:min_pc,
                           pK = optimal.pk, # neighborhood size used to compute artificial nearest neighbors
                           nExp = nExp.poi.adj) # number of expected real doublets
  colnames(sample@meta.data)[grepl("DF.classifications.*", colnames(sample@meta.data))] <- "doublet_finder"

  double_finder_res <- sample@meta.data["doublet_finder"]
  double_finder_res <- rownames_to_column(double_finder_res, "row_names")
  return(double_finder_res)
}
