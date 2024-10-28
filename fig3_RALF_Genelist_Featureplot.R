### library --------------------------------------------------------------------
library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(viridis)
library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(tidyverse)
library(reshape2)
library(ggdist)

### Load data ------------------------------------------------------------------

## Import genelist
genelist <- read.table("RALF.selected.txt", header = F)$V1

rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/At/STO_embryo_scRNA/all_time_merge/allChr_v1/rmR1802/all_sample_subCluster_CellType/At_merge.seurat_clusters_renamed.rds")

p1 <- FeaturePlot(rds, features = genelist, pt.size = 0.1, order = T, raster = F, ncol = 5) & scale_color_viridis(option="rocket", begin = 1, end = 0.3)
pdf("RALF.selected.rasterF.231010.v2.pdf", 3*5, 3*2)
print(p1)
dev.off()

p1 <- FeaturePlot(rds, features = genelist, pt.size = 0.1, order = T, raster = T, ncol = 5) & scale_color_viridis(option="rocket", begin = 1, end = 0.3)
pdf("RALF.selected.rasterT.231010.v2.pdf", 3*5, 3*2)
print(p1)
dev.off()
