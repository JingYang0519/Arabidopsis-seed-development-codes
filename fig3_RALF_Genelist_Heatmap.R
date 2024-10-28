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
set.seed(123)

##### 10 RALFs
### Load data ------------------------------------------------------------------

genelist <- read.table("10.RALF.list.order.txt", header = T)
rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/At/STO_embryo_scRNA/all_time_merge/allChr_v1/rmR1802/all_sample_subCluster_CellType/At_merge.seurat_clusters_renamed.rds")

avgExp <- AverageExpression(rds, features = genelist$GeneID, group.by = "seurat_clusters_renamed")
avgExp <- as.data.frame(avgExp)
names(avgExp) <- gsub("RNA.","",names(avgExp))
order <- c("1_9.11h_Reproductive_cell_1_Egg_cell",
           "1_9.11h_Reproductive_cell_1_Sperm","2_28h_Reproductive_cell_1_Sperm",
           "1_9.11h_Reproductive_cell_2_Central_cell",
           "1_9.11h_Reproductive_cell_2_Synergid_cell", 
           "2_28h_Zygote.Basal_cell",
           "2_28h_Embryo","3_48h_Embryo","4_globular_Embryo","5_heart_Embryo","6_torpedo_Embryo","7_bent_Embryo","8_cotyledon_Embryo","9_late_cotyledon_Embryo",
           "2_28h_Endosperm","3_48h_Endosperm","4_globular_Endosperm","5_heart_Endosperm","6_torpedo_Endosperm","7_bent_Endosperm","8_cotyledon_Endosperm","9_late_cotyledon_Endosperm",
           "1_9.11h_Seed_coat","2_28h_Seed_coat","3_48h_Seed_coat","4_globular_Seed_coat","5_heart_Seed_coat","6_torpedo_Seed_coat","7_bent_Seed_coat","8_cotyledon_Seed_coat","9_late_cotyledon_Seed_coat",
           "3_48h_Suspensor","4_globular_Suspensor")
avgExp <- avgExp[,order]
order.tp <- unlist(lapply(order, function(x) substr(x, 1, 2)))
order.tp <- gsub("1_", "1.9-11h", order.tp)
order.tp <- gsub("2_", "2.28h", order.tp)
order.tp <- gsub("3_", "3.48h", order.tp)
order.tp <- gsub("4_", "4.Globular", order.tp)
order.tp <- gsub("5_", "5.Heart", order.tp)
order.tp <- gsub("6_", "6.Torpedo", order.tp)
order.tp <- gsub("7_", "7.Bent", order.tp)
order.tp <- gsub("8_", "8.Cotyledon", order.tp)
order.tp <- gsub("9_", "9.Late-cotyledon", order.tp)

### Define annotation information ----------------------------------------------
ha_column = HeatmapAnnotation(df = data.frame(TimePoints = order.tp,
                                              CellType = c(rep("1.ReproductiveCell_1_EggCell",1),
                                                           rep("2.ReproductiveCell_1_Sperm",2),
                                                           rep("3.ReproductiveCell_2_CentralCell",1),
                                                           rep("4.ReproductiveCell_2_SynergidCell",1),
                                                           rep("5.Zygote,basal",1),
                                                           # rep("6.ApicalCell",1),
                                                           rep("6.Embryo",8),
                                                           rep("7.Endosperm",8),
                                                           rep("8.SeedCoat",9),
                                                           rep("9.Suspensor",2))),
                              col = list(TimePoints = c("1.9-11h" = "#FDE725FF",
                                                        "2.28h" = "#B4DE2CFF",
                                                        "3.48h" = "#6DCD59FF",
                                                        "4.Globular" = "#35B779FF",
                                                        "5.Heart" = "#1F9E89FF",
                                                        "6.Torpedo" = "#26828EFF",
                                                        "7.Bent" = "#31688EFF",
                                                        "8.Cotyledon" = "#3E4A89FF",
                                                        "9.Late-cotyledon" = "#482878FF"),
                                         CellType = c("1.ReproductiveCell_1_EggCell" = "#7CC08C",
                                                      "2.ReproductiveCell_1_Sperm" = "#bce07e",
                                                      "3.ReproductiveCell_2_CentralCell" = "#a3b9f0",
                                                      "4.ReproductiveCell_2_SynergidCell" = "#FDE377",
                                                      "5.Zygote,basal" = "#3ca372",
                                                      "6.Embryo" = "#357B6E",
                                                      "7.Endosperm" = "#63a4cf",
                                                      "8.SeedCoat" = "#F4A370",
                                                      # "6.ApicalCell" = "#3fab98",
                                                      "9.Suspensor" = "#d962c2")))

#ha_row = HeatmapAnnotation(RALF.info = genelist[,3], col = list(RALF.info = c("Not_MP1" = "#444654", "MP1.2" = "#147a86", "MP1.1_unique" = "#ea595e")), which = "row")

ha_row = HeatmapAnnotation(RALF.info = genelist[,3], col = list(RALF.info = c("MP1.2" = "#4472C4", "MP1.1_unique" = "#C00000")), which = "row")
### Normalize avg expression ---------------------------------------------------
normalize <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

normalized_avgExp <- t(apply(avgExp, 1, normalize))
normalized_avgExp <- as.data.frame(normalized_avgExp)
row.names(normalized_avgExp) <- genelist$RALFs

h1 <- Heatmap(as.matrix(normalized_avgExp), 
              # col = colorRamp2(seq(0,1,0.1),c("#FFFFFF",rocket(10,begin = 1,end = 0.3))),
              col = colorRamp2(seq(0,1,0.1),c("#fcf7f2",rocket(10,begin = 1,end = 0.3))),
              heatmap_height = unit(100, "mm"),
              heatmap_width = unit(100, "mm"),
              #clustering_method_rows = "complete",
              #clustering_method_columns = "complete",
              top_annotation = ha_column,
              left_annotation = ha_row,
              cluster_rows = F, cluster_columns = F, name = "Normalized average expression",
              show_row_names = F, show_column_names = F)
pdf("AverageExpression.10.RALF.231010.v2.pdf",14,10)
print(h1)
dev.off()













# ##### Backup ----
# ### Load data ------------------------------------------------------------------
# 
# genelist <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/09.Expression.heatmaps/All.RALF.list.order.txt", header = T)
# rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/At/STO_embryo_scRNA/all_time_merge/allChr_v1/rmR1802/all_sample_subCluster_CellType/At_merge.seurat_clusters_renamed.rds")
# 
# avgExp <- AverageExpression(rds, features = genelist$GeneID, group.by = "seurat_clusters_renamed")
# avgExp <- as.data.frame(avgExp)
# names(avgExp) <- gsub("RNA.","",names(avgExp))
# order <- c("01_9.11h_Reproductive_cell_1_Egg_cell",
#            "01_9.11h_Reproductive_cell_1_Sperm","02_24h_Reproductive_cell_1_Sperm","03_28h_Reproductive_cell_1_Sperm",
#            "01_9.11h_Reproductive_cell_2_Central_cell",
#            "01_9.11h_Reproductive_cell_2_Synergid_cell", 
#            "02_24h_Zygote",
#            "02_24h_Apical_cell",
#            "03_28h_Embryo","04_48h_Embryo","05_globular_Embryo","06_heart_Embryo","07_torpedo_Embryo","08_bent_Embryo","09_cotyledon_Embryo","10_late_cotyledon_Embryo",
#            "02_24h_Endosperm","03_28h_Endosperm","04_48h_Endosperm","05_globular_Endosperm","06_heart_Endosperm","07_torpedo_Endosperm","08_bent_Endosperm","09_cotyledon_Endosperm","10_late_cotyledon_Endosperm",
#            "01_9.11h_Seed_coat","02_24h_Seed_coat","03_28h_Seed_coat","04_48h_Seed_coat","05_globular_Seed_coat","06_heart_Seed_coat","07_torpedo_Seed_coat","08_bent_Seed_coat","09_cotyledon_Seed_coat","10_late_cotyledon_Seed_coat",
#            "04_48h_Suspensor","05_globular_Suspensor")
# avgExp <- avgExp[,order]
# order.tp <- unlist(lapply(order, function(x) substr(x, 1, 2)))
# order.tp <- gsub("01", "1.9-11h", order.tp)
# order.tp <- gsub("02", "2.24h", order.tp)
# order.tp <- gsub("03", "3.28h", order.tp)
# order.tp <- gsub("04", "4.48h", order.tp)
# order.tp <- gsub("05", "5.Globular", order.tp)
# order.tp <- gsub("06", "6.Heart", order.tp)
# order.tp <- gsub("07", "7.Torpedo", order.tp)
# order.tp <- gsub("08", "8.Bent", order.tp)
# order.tp <- gsub("09", "9.Cotyledon", order.tp)
# order.tp <- gsub("10", "90.Late-cotyledon", order.tp)
# 
# ### Define annotation information ----------------------------------------------
# ha_column = HeatmapAnnotation(df = data.frame(TimePoints = order.tp,
#                                                CellType = c(rep("1.ReproductiveCell_1_EggCell",1),
#                                                             rep("2.ReproductiveCell_1_Sperm",3),
#                                                             rep("3.ReproductiveCell_2_CentralCell",1),
#                                                             rep("4.ReproductiveCell_2_SynergidCell",1),
#                                                             rep("5.Zygote",1),
#                                                             rep("6.ApicalCell",1),
#                                                             rep("7.Embryo",8),
#                                                             rep("8.Endosperm",9),
#                                                             rep("9.SeedCoat",10),
#                                                             rep("90.Suspensor",2))),
#                                col = list(TimePoints = c("1.9-11h" = "#FDE725FF",
#                                                          "2.24h" = "#B4DE2CFF",
#                                                          "3.28h" = "#6DCD59FF",
#                                                          "4.48h" = "#35B779FF",
#                                                          "5.Globular" = "#1F9E89FF",
#                                                          "6.Heart" = "#26828EFF",
#                                                          "7.Torpedo" = "#31688EFF",
#                                                          "8.Bent" = "#3E4A89FF",
#                                                          "9.Cotyledon" = "#482878FF",
#                                                          "90.Late-cotyledon" = "#440154FF"),
#                                           CellType = c("1.ReproductiveCell_1_EggCell" = "#7CC08C",
#                                                        "2.ReproductiveCell_1_Sperm" = "#bce07e",
#                                                        "3.ReproductiveCell_2_CentralCell" = "#a3b9f0",
#                                                        "4.ReproductiveCell_2_SynergidCell" = "#FDE377",
#                                                        "9.SeedCoat" = "#F4A370",
#                                                        "6.ApicalCell" = "#3fab98",
#                                                        "8.Endosperm" = "#63a4cf",
#                                                        "5.Zygote" = "#3ca372",
#                                                        "7.Embryo" = "#357B6E",
#                                                        "90.Suspensor" = "#d962c2")))
# 
# ha_row = HeatmapAnnotation(RALF.info = genelist[,3], col = list(RALF.info = c("Not_MP1" = "#444654", "MP1.2" = "#147a86", "MP1.1_unique" = "#ea595e")), which = "row")
# 
# ### Normalize avg expression ---------------------------------------------------
# normalize <- function(x) {
#   (x - min(x)) / (max(x) - min(x))
# }
# 
# normalized_avgExp <- t(apply(avgExp, 1, normalize))
# normalized_avgExp <- as.data.frame(normalized_avgExp)
# row.names(normalized_avgExp) <- genelist$RALFs
# 
# h1 <- Heatmap(as.matrix(normalized_avgExp), 
#               col = colorRamp2(c(0,1),c("#FFFFFF","#D44520")),
#               heatmap_height = unit(150, "mm"),
#               heatmap_width = unit(200, "mm"),
#               #clustering_method_rows = "complete",
#               #clustering_method_columns = "complete",
#               top_annotation = ha_column,
#               left_annotation = ha_row,
#               cluster_rows = F, cluster_columns = F, name = "Normalized average expression",
#               show_row_names = F, show_column_names = F)
# h2 <- Heatmap(as.matrix(normalized_avgExp), 
#               # col = colorRamp2(c(0,0.5,1),c("#3484BE","#FFFFFF","#f48d70")),
#               col = colorRamp2(c(0,1),c("#FFFFFF","#D44520")),
#               heatmap_height = unit(150, "mm"),
#               heatmap_width = unit(200, "mm"),
#               #clustering_method_rows = "complete",
#               #clustering_method_columns = "complete",
#               top_annotation = ha_column,
#               left_annotation = ha_row,
#               cluster_rows = F, cluster_columns = F, name = "Normalized average expression",
#               show_row_names = T, show_column_names = F)
# pdf("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/09.Expression.heatmaps/AverageExpression.all.RALF.v2.pdf",14,10)
# print(h1)
# print(h2)
# dev.off()