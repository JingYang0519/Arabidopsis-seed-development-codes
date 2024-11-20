#------------------------------------------------#
#          Cell communication Analysis           #
#     Editor: Ruiying Chen   TEL:18995624500     #
#              Edited on: 2023.09.14             #
#------------------------------------------------#

### Reference ------------------------------------------------------------------
# library(CommPath)
# trace('findLRpath', edit = T, where = asNamespace("CommPath"))
# trace('scorePath', edit = T, where = asNamespace("CommPath"))
# trace('pathHeatmap', edit = T, where = asNamespace("CommPath"))

### Library --------------------------------------------------------------------
# BiocManager::install("scde")
# BiocManager::install("DEsingle")
# BiocManager::install("MAST")
# install.packages("edgebundleR")

library(PlantPhoneDB)
library(Seurat)
library(dplyr)
library(igraph)
library(RColorBrewer)
library(viridis)
library(ggplot2)
library(tidyverse)
library(edgebundleR)
library(ggplotify)
library(circlize)
library(RColorBrewer)
library(cowplot)
library(ggsci)
library(GSVA)
library(ComplexHeatmap)
library(viridis)
library(dendsort)
# library(org.At.tair.db)
# keytypes(org.At.tair.db) # org.At.tair ARACYC Mappings between TAIR identifiers and KEGG pathway identifiers

### 0. RDS data ---------------------------------------------------------

rds <- readRDS("At_merge.seurat_clusters_renamed.rds")
# sub_rds <- subset(rds, Organ == args[1])
sub_rds <- subset(rds, Organ == "1_9_11h")
sub_rds$seurat_clusters_renamed <- gsub("Reproductive_cell_1_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("Reproductive_cell_2_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("1_9-11h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("2_28h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("3_48h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("4_globular_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("5_heart_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("6_torpedo_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("7_bent_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("8_cotyledon_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("9_late_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("9_late_cotyledon_","",sub_rds$seurat_clusters_renamed)
sub_rds <- subset(sub_rds, seurat_clusters_renamed %in% c("Seed_coat","Central_cell","Sperm","Synergid_cell",
                                                          "Egg_cell","Apical_cell","Zygote","Endosperm","Embryo",
                                                          "Suspensor"))
sub_rds2 <- SCTransform(sub_rds, verbose = FALSE)
## Method 1: Average expression
data <-as.data.frame(AverageExpression(sub_rds2, group.by = "seurat_clusters_renamed", assays = "SCT"))
colnames(data) <- gsub("SCT.","",colnames(data))
data <- as.matrix(data)
## Method 2: Raw data
# data <- as.matrix(sub_rds2@assays$SCT@data)

### 1. Statistics -------------------------------------------------------
LRp.cmb <- read.delim("LRp_CellType.revised/LRp.cmb.p0.05.3154.revised.tsv", header = T)
head(LRp.cmb)

### Anno by database
# goannot <- select(org.At.tair.db, keys=keys(org.At.tair.db), columns="ARACYC")
# head(goannot)
# gs <- split(goannot$TAIR, goannot$ARACYC)
goannot <- read.delim("02.Araport11.Mar92021.geneid_GO_annotation.xls", header = T)
gs <- split(goannot$GeneID, goannot$Description)

### Test on 9-11h cells
unique(LRp.cmb$TimePoint)
LRp.cmb.test <- LRp.cmb[which(LRp.cmb$TimePoint == "01_9_11h"),]
LRp.cmb.test <- LRp.cmb.test[which(LRp.cmb.test$Ligands_cell == "Sperm"),]
lrlist <- unique(c(LRp.cmb.test$Ligands,LRp.cmb.test$Receptors))

which.overlap.list <- unlist(lapply(gs, function(x) {
  if (any(x %in% lrlist)) {
    return(TRUE)
  }
  else {
    return(FALSE)
  }
}))
path.lr.list <- gs[which(which.overlap.list)]
head(path.lr.list)
gsva.es <- gsva(data, path.lr.list, verbose=FALSE, method="zscore", parallel.sz = 4)
gsva.es <- as.data.frame(gsva.es)

ntop <- 3
top_rows <- vector("list", ncol(gsva.es))
for (col in 1:ncol(gsva.es)) {
  top_rows[[col]] <- rownames(gsva.es[order(-gsva.es[, col]), ][1:ntop, ])
}
idx <- unique(unlist(top_rows))
gsva.es.top <- gsva.es[which(row.names(gsva.es) %in% idx),]

# saveRDS(gsva.es.top, "Test.gsva.es.top.rds")
gsva.es.top <- readRDS("Test.gsva.es.top.rds")
set.seed(123)

row_dend = dendsort(hclust(dist(gsva.es.top)))
col_dend = dendsort(hclust(dist(t(gsva.es.top))))

h <- Heatmap(as.matrix(gsva.es.top), 
              col = colorRamp2(seq(min(gsva.es.top),max(gsva.es.top),(1/10)*(max(gsva.es.top)-min(gsva.es.top))),
                               magma(11, begin = 1,end = 0.3)),
              # heatmap_height = unit(35*ntop, "mm"),
             heatmap_height = unit(35*ntop, "mm"),
              heatmap_width = unit(100, "mm"),
             column_names_rot = 45,
             cluster_rows = row_dend, cluster_columns = col_dend, 
             name = "GSVA Zscore",
              show_row_names = T, show_column_names = T, 
             column_dend_reorder = T, row_dend_reorder = T,
             show_column_dend = FALSE, show_row_dend = FALSE,
             heatmap_legend_param = list(
               at = c(min(gsva.es.top), max(gsva.es.top)),
               labels = c("Low", "High"),
               title = "GSVA Zscore",
               legend_height = unit(4, "cm"),
               direction = "horizontal",
               title_position = "topcenter"
             ))
pdf("Test.pdf",20,10)
draw(h, heatmap_legend_side = "bottom")
dev.off()
