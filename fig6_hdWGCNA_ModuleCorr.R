library(DOSE)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(igraph)
library(org.At.tair.db)   
keytypes(org.At.tair.db)
library(Seurat)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(dplyr)


# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("f05.AnnoLabel.ModuleEigengenesConnectivity.rds")

pdf("p10.ModuleCorrelogram.pdf",6,5)
print(ModuleCorrelogram(seurat_obj))
dev.off()