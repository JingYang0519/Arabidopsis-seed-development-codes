library(Seurat) 
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(circlize)
library(ComplexHeatmap)
library(viridis)

rds <- readRDS("AnnoLabel.merge.rds")
umap <- rds[["umap"]]@cell.embeddings
row.names(umap) <- gsub("-","_",row.names(umap))
data <- cbind(rds@meta.data,umap)

data$DetailAnno <- gsub("Seed.coat", "Seed.Coat", data$DetailAnno)
print("Doing plot...")
### Palette
TimePoints.col = c("#FDE725FF","#B4DE2CFF","#6DCD59FF","#35B779FF","#1F9E89FF",
                   "#26828EFF","#31688EFF","#3E4A89FF","#482878FF")
CellType.col = c("#357B6E","#63a4cf","#934F3D",
                 "#7CC08C","#bce07e","#a3b9f0",
                 "#FDE377", "#F4A370",
                 "#d962c2","#AFAFAF","#3ca372")

CellType.sub = c("#a3b9f0",
                 "#7CC08C",
                 colorRampPalette(c("#3d8d7e","#26574e"))(4),
                  rev(colorRampPalette(c("#77afd5","#176da6"))(4)),
                 "#934F3D",
                 rev(colorRampPalette(c("#f8c19f","#fa945c"))(5)),
                 "#bce07e",
                 "#d962c2",
                 "#FDE377", 
                 "#AFAFAF",
                 "#3ca372")

pdf("Test.dimplot.DetailAnno.v2.pdf",15,10)
ggplot(data, aes(x = UMAP_1, y = UMAP_2, color = DetailAnno)) + geom_point(size = 0.5) +
  scale_color_manual(values = CellType.sub) +
  theme_classic()
dev.off()

### ----------------------------------------------------------------------------
library(Seurat) 
library(stringr)
library(RColorBrewer)
library(ggplot2)

rds2 <- readRDS("AnnoLabel.merge.rds")
row.names(rds2@meta.data) <- gsub("embryo_9_11h", "embryo_9-11h", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("late_cotyledon", "late-cotyledon", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("07_torpedo_", "07_torpedo-", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("06_heart_", "06_heart-", row.names(rds2@meta.data))
row.names(rds2@meta.data) <- gsub("08_walking_stick_", "08_walking_stick-", row.names(rds2@meta.data))

rds2$label2 <- paste0(rds2$Organ, ".", rds2$DetailAnno)
avgExp <- AverageExpression(rds2, group.by = "label2")
avgExp <- as.data.frame(avgExp)
colnames(avgExp) <- gsub("RNA.", "", colnames(avgExp))
df.s <- scale(avgExp, center = T, scale = T) 
df.spearman1 <- as.data.frame(cor(df.s, method="spearman"))
df.spearman1["99_label",] <- c(rep(1,71))
df.spearman1["99_label2",] <- c(rep(1,71))
df <- as.data.frame(table(rds2$label2))
tmp <- data.frame(Var1 = c("99_label","99_label2"),Freq = c(0,0))
df <- rbind(df, tmp)
row.names(df) <- df$Var1

# h1 <- Heatmap(df.spearman1, 
#               col = rocket(10,begin = 1,end = 0.3),
#               # heatmap_height = unit(100, "mm"),
#               # heatmap_width = unit(100, "mm"),
#               cluster_rows = T, cluster_columns = T, name = "Gene expression\n Spearman correlation",
#               show_row_names = T, show_column_names = T)
# pdf("RNA.corrHeatmap.pdf", width = 15, height = 15)
# h1
# dev.off()

level_Organ <- c(rep("1_9_11h",9), 
                rep("2_28h" ,8),
                rep("3_48h",7),
                rep("4_globular",7),
                rep("5_heart",9),
                rep("6_torpedo",9),
                rep("7_bent",9),
                rep("8_cotyledon",6),
                rep("9_late_cotyledon",7),
                rep("99_label",2)) %>% factor()

df.spearman1 <- rbind(df.spearman1, t(df))
col1 = colorRamp2(c(0.6, 0.8, 1), c("navy", "white", "firebrick3"))
pdf("RNA.corrHeatmap.v2.pdf", width = 50, height = 50)
circos.heatmap.initialize(as.matrix(df.spearman1), split = level_Organ)
circos.trackHist(df$Var1, df$Freq, bin.size = 0.2)
circos.heatmap(as.matrix(df.spearman1[1:72,]), 
               col = col1, 
               dend.side = "inside",
               cluster = T,
               rownames.side = "outside", 
               track.height = 0.5,
               cell.lwd = 0.8,
               split = level_Organ)
circos.clear()
dev.off()