library(Seurat) 
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(viridis)

rds2 <- readRDS("AnnoLabel.merge.rds")
row.names(rds2@meta.data) <- gsub("embryo_9_11h", "embryo_9-11h", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("late_cotyledon", "late-cotyledon", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("07_torpedo_", "07_torpedo-", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("06_heart_", "06_heart-", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("08_walking_stick_", "08_walking_stick-", row.names(rds2@meta.data))

rds2$DetailAnno <- gsub("Seed.coat", "Seed.Coat", rds2$DetailAnno)

rds2$label2 <- paste0(rds2$Organ, ".", rds2$DetailAnno)
avgExp <- AverageExpression(rds2, group.by = "label2")
avgExp <- as.data.frame(avgExp)
colnames(avgExp) <- gsub("RNA.", "", colnames(avgExp))
df.s <- scale(avgExp, center = T, scale = T) 

df.spearman1 <- as.data.frame(cor(df.s, method="spearman"))
df <- as.data.frame(table(rds2$label2))
row.names(df) <- df$Var1
df.spearman1 <- as.data.frame(cbind(df.spearman1, df))
write.table(df.spearman1, "df.spearman1.beforelabel.tsv", quote = F, sep = "\t", row.names = T, col.names = T)

# h1 <- Heatmap(df.spearman1, 
#               col = rocket(10,begin = 1,end = 0.3),
#               # heatmap_height = unit(100, "mm"),
#               # heatmap_width = unit(100, "mm"),
#               cluster_rows = T, cluster_columns = T, name = "Gene expression\n Spearman correlation",
#               show_row_names = T, show_column_names = T)
# pdf("RNA.corrHeatmap.pdf", width = 15, height = 15)
# h1
# dev.off()

data <- read.delim("df.spearman1.afterlabel.tsv", header = T)
row.names(data) <- c(colnames(data[1:68]),"99_lab1","99_lab2")
col1 = colorRamp2(c(0.6, 0.8, 1), c("navy", "white", "firebrick3"))

level_Organ <- c(rep("1_9_11h",9), 
                rep("2_28h" ,8),
                rep("3_48h",7),
                rep("4_globular",7),
                rep("5_heart",8),
                rep("6_torpedo",8),
                rep("7_bent",9),
                rep("8_cotyledon",6),
                rep("9_late_cotyledon",6),
                rep("99_label",2)) %>% factor()

h <- circos.heatmap(as.matrix(data)[,1:68], 
               col = col1, 
               dend.side = "inside",
               dend.track.height = 0.07,
               cluster = T,
               rownames.side = "outside", 
               track.height = 0.5,
               cell.lwd = 0.1,
               split = level_Organ)

pdf("RNA.corrHeatmap.dendrogram.pdf", width = 50, height = 50)
circos.heatmap(as.matrix(data)[,1:68], 
               col = col1, 
               dend.side = "inside",
               dend.track.height = 0.07,
               cluster = T,
               rownames.side = "outside", 
               track.height = 0.5,
               cell.lwd = 0.1,
               split = level_Organ)
circos.clear()
dev.off()

pdf("RNA.corrHeatmap.label.ring.pdf", width = 50, height = 50)
circos.heatmap(as.matrix(data)[,1:3], 
               col = col1, 
               dend.side = "inside",
               dend.track.height = 0.07,
               cluster = T,
               rownames.side = "outside", 
               track.height = 0.5,
               cell.lwd = 0.1,
               split = level_Organ)
circos.clear()
dev.off()

CELL_META$order