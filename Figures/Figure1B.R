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

all_Whole <- qread(paste0(path, "all_Whole.qs"))

features = c("TLR3")

DotPlot(all_Whole, features = features, group.by = "tissue1", dot.scale = 10) +
    scale_y_discrete(limits = rev(levels(all_Whole$tissue1))) +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) +
    theme(axis.text.y = element_text(size = 20))
