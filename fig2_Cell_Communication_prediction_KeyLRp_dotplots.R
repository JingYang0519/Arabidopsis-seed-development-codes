library(Seurat)
library(ggplot2)
library(viridis)
library(patchwork)

rds <- readRDS("At_merge.seurat_clusters_renamed.rds")

unique(rds$Organ)
sub_rds <- rds
sub_rds <- subset(sub_rds, Organ %in% c("1_9_11h","2_28h"))
sub_rds$seurat_clusters_renamed <- gsub("1_9-11h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("2_28h_","",sub_rds$seurat_clusters_renamed)

### Import genelist ------------------------------------------------------------

## Test v2
p1 <- DotPlot(object = sub_rds, scale.by = "size", group.by = "CellClass", features = rev(unique(genelist$Ligand))) + #cols = c("#ffffff", "#448444")) +
  coord_flip() +
  scale_x_discrete(position = "top") +
  theme_bw() +
  scale_y_discrete(position = "right") +
    scale_color_viridis(option="rocket", begin = 1, end = 0.6)

p3 <- DotPlot(object = sub_rds, scale.by = "size", group.by = "CellClass", features = rev(unique(genelist$Receptor))) + #cols = c("#ffffff", "#448444")) +
  coord_flip() +
  theme_bw() +
  scale_y_discrete(position = "right")  + NoLegend() +
    scale_color_viridis(option="rocket", begin = 1, end = 0.6)

pdf("Comm.SEE.Genes.Ligand.pdf", 30, 10)
p1 + p3
dev.off()
