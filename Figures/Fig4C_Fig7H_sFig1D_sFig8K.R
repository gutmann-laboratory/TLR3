TAM_Clusters_sub <- qread(paste0(path, "OPG_TAM_Clusters_sub.qs"))

# Fig4C
features = c("Tmem119", "Ms4a7", "Fgr")
DotPlot(TAM_Clusters_sub, features = rev(features), group.by = "celltype_new", dot.scale = 10, idents = "FMC") +
 theme(axis.title = element_blank()) +
 theme(axis.text.y = element_text(size = 25, face = "italic")) +
 theme(axis.text.x = element_text(size = 25, angle = 90, hjust = 1, vjust = 0.5)) +
 coord_flip()

# Fig7H
features = c("Ifih1")
DotPlot(TAM_Clusters_sub, features = features, group.by = "celltype_new", dot.scale = 10, idents = "FMC") +
 theme(axis.title = element_blank()) +
 theme(axis.text.y = element_text(size = 25, face = "italic")) +
 theme(axis.text.x = element_text(size = 25, angle = 90, hjust = 1, vjust = 0.5))


# sFig1D
features = c("Tlr3")
DotPlot(TAM_Clusters_sub, features = features, group.by = "celltype_new", dot.scale = 10, idents = "FMC") +
 theme(axis.title = element_blank()) +
 theme(axis.text.y = element_text(size = 25, face = "italic")) +
 theme(axis.text.x = element_text(size = 25, angle = 90, hjust = 1, vjust = 0.5))

# sFig8K
FMC <- subset(TAM_Clusters_sub, condition == "FMC")

features = c("Ifih1", "Ddx58") 

DotPlot(FMC, features = features, group.by = "condition") +
    theme(axis.title = element_blank()) +
    theme(axis.text.x = element_text(size = 20, angle = 45, hjust = 1)) +
    theme(axis.text.y = element_text(size = 20))
