library(reshape2)
library(scales)
library(ggplot2)
library(dplyr)
library(gplots)
library(viridis)
library(Seurat)

celltrajectory <- readRDS("refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")

time_point <- unique(celltrajectory$Organ)
time_point

names(celltrajectory@meta.data)
unique(celltrajectory$celltype0515)

rds <- subset(celltrajectory, celltype0515 !="Sp")

rds
names(rds@meta.data)

deg_results <- list()

edges <- read.table("endo_edge.txt", header = FALSE, stringsAsFactors = FALSE)
if (nrow(edges) == 0) {
  stop("Error: endo_edge.txt is empty or not read correctly.")
}

for (i in unique(edges$V1)) {

  related_rows <- edges[edges$V1 == i, ]
  print(i)
  print(related_rows)

  cell_types <- unique(sapply(strsplit(as.character(related_rows$V2), ":"), `[`, 2))
  print(cell_types)

  time_points <- unique(sapply(strsplit(as.character(related_rows$V2), ":"), `[`, 1))
  print(time_points)

  time_point <- time_points[1]
  print(time_point)

  cells_to_keep <- subset(rds, subset = Organ == time_point & celltype0515 %in% cell_types)


  unique(Idents(cells_to_keep))

  cells_to_keep <- SetIdent(cells_to_keep, value = cells_to_keep@meta.data$celltype0515)

  unique(Idents(cells_to_keep))

  degs <- FindMarkers(cells_to_keep, ident.1 = cell_types[1], ident.2 = cell_types[2])

  deg_results[[i]] <- degs
}

all_degs <- data.frame()

for (group_name in names(deg_results)) {

  degs <- deg_results[[group_name]]
  
  if (!"gene" %in% colnames(degs)) {
    degs$gene <- rownames(degs)
  }

  degs$group <- group_name

  all_degs <- rbind(all_degs, degs)
}

write.csv(all_degs, "DEG_results.csv", row.names = FALSE)
