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

all_SO <- qraed(path, "all_SO.qs")

# Figure 1C
features = "TLR3"

dot_data <- DotPlot(all_SO, features = features)
df <- dot_data$data

df <- df %>% separate(id, into = c("Pathology", "celltype"), sep = "_")

df$Pathology <- factor(df$Pathology, levels = c("PA", "DIPG", "Grade II Astrocytoma", "GBM"))
df$celltype <- factor(df$celltype, levels = c("Tumor", "Oligodendrocyte", "TAM"))

ggplot(df, aes(x = Pathology, y = celltype)) + 
  geom_point(aes(size = pct.exp, color = avg.exp.scaled)) + 
　scale_size(range = c(3, 10)) + 
  scale_color_gradient(low = "lightgrey", high = "blue") + 
  labs(title = "TLR3", x = "Tumor", y = "CellType") + 
  theme_cowplot() +
  theme(axis.title = element_blank()) +
  theme(axis.text.x = element_text(size = 20, vjust = 0.5)) +
  theme(axis.text.y = element_text(size = 20)) +
  scale_x_discrete(labels = c(
    "PA" = "PA",
    "DIPG" = "DIPG",
    "Grade II Astrocytoma" = "Grade II\nAstrocytoma",
    "GBM" = "GBM"
  ))





