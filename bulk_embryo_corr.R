library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/04.ZY/CellClass.rds")#读取单细胞RDS
rds1 <- rds
rds1$CellClass <- gsub("03_28h_Zygote","03_28h_Embryo",rds1$CellClass)#将28h的Zygote改为Embryo，即合并数据
DEG <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/zhangyi9/SCRNA_cor/bulk_cor/DEG_cor/DEG_wilcox_Embryo_Zygote.tsv", header = T)#读取单细胞DEG
DEG$cluster <- gsub("03_28h_Zygote","03_28h_Embryo",DEG$cluster)
DEG.list <- DEG$gene # head(DEG.list)
bulk <- read.table("/jdfsbjcas1/ST_BJ/P21Z28400N0234/zhangyi9/SCRNA_cor/bulk_cor/Bulk_counts_mean.txt", header = T)
row.names(bulk) <- bulk$TAIR10_id
inter <- intersect(row.names(bulk),DEG.list) # Intersection gene list of bulk genelist and scDEGs

bulk.select <- bulk[inter,]
bulk.select <- bulk.select[,-1]
rds.embryo <- subset(rds1,subset=CellClass %in% c("02_24h_Zygote","03_28h_Embryo","04_48h_Embryo","05_globular_Embryo","06_heart_Embryo","07_torpedo_Embryo","08_bent_Embryo","09_cotyledon_Embryo","10_late_cotyledon_Embryo"))
scRNA.select <- as.data.frame(AverageExpression(rds.embryo, assays="RNA", group.by="CellClass", features = inter)) # Calculate average expression for each "CellClass"

#method <- c("spearman","pearson"
corr <- cor(bulk.select, scRNA.select, method = "spearman")
colnames(corr) <- gsub("^[^_]*_", "", colnames(corr))
write.table(corr, file = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/zhangyi9/SCRNA_cor/bulk_cor/DEG_cor/sc_Bulk_DEG_cor_230626/Bulk_sc_DEG_cor_0628.tsv", quote = F, sep = "\t",row.names = T,col.names=T)
f1 <- colorRamp2(breaks = c(min(corr),0,0.7,max(corr)),c("#74ADD1","#FFFFBF","burlywood1","#D73027"))
corr <- corr[,c("24h_Zygote","28h_Embryo","48h_Embryo","globular_Embryo","heart_Embryo","torpedo_Embryo","bent_Embryo","cotyledon_Embryo","late_cotyledon_Embryo")]
# matrix.corr <- matrix(corr)
s2 <- Heatmap(as.matrix(corr), name = "spearman corr", 
                column_title = " ",row_title = " ",
                col = f1,
                cluster_rows = FALSE, cluster_columns = FALSE,
                show_column_names = TRUE,show_row_names = TRUE,
                row_names_side = "right",
                row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
                column_dend_side = "top",column_title_rot = 0,
                column_title_side = "top",
                column_names_rot = 45)
#pheatmap画法
s2 <- pheatmap(corr, border_color = "black",
               color = f1,
               cluster_row = FALSE,
               cluster_cols = FALSE) 
  pdf(paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/zhangyi9/SCRNA_cor/bulk_cor/DEG_cor/sc_Bulk_DEG_cor_230626/Bulk_scRNA.Embryo.Zygote.corr.spearman.pdf"),9,7)
  print(s2)
  dev.off()