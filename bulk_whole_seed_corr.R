library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
rds <- readRDS("CellClass.rds")
rds1 <- rds

rds1$CellClass <- gsub("^02_24h_.*","02_24h",rds1$CellClass)
rds1$CellClass <- gsub("^03_28h_.*","03_28h",rds1$CellClass)
rds1$CellClass <- gsub("^04_48h_.*","04_48h",rds1$CellClass)
rds1$CellClass <- gsub("^05_glo.*","05_globular",rds1$CellClass)
rds1$CellClass <- gsub("^06_h.*","06_heart",rds1$CellClass)
rds1$CellClass <- gsub("^07.*","07_torpedo",rds1$CellClass)
rds1$CellClass <- gsub("^08.*","08_bent",rds1$CellClass)
rds1$CellClass <- gsub("^09.*","09_cotyledon",rds1$CellClass)
rds1$CellClass <- gsub("^10.*","10_late_cotyledon",rds1$CellClass)

DEG <- read.table("SCRNA_cor/bulk_Whole_Seed_cor/sc_DEG_wilcox_WS.tsv", header = T)
DEG.list <- DEG$gene # head(DEG.list)
bulk <- read.table("SCRNA_cor/bulk_Whole_Seed_cor/bulk_whole_seed_mean.tsv", header = T)
subset_bulk <- subset(bulk, select = c("gene","preglobular","heart","linear_cotyledon","mature_green"))
row.names(subset_bulk) <- subset_bulk$gene
inter <- intersect(row.names(subset_bulk),DEG.list) # Intersection gene list of bulk genelist and scDEGs

bulk.select <- subset_bulk[inter,]
bulk.select <- bulk.select[,-1]
rds.embryo <- subset(rds1,subset=CellClass %in% c("02_24h","03_28h","04_48h","05_globular","06_heart","07_torpedo","08_bent","09_cotyledon","10_late_cotyledon"))
scRNA.select <- as.data.frame(AverageExpression(rds.embryo, assays="RNA", group.by="CellClass", features = inter)) # Calculate average expression for each "CellClass"

#method <- c("spearman","pearson"
corr <- cor(bulk.select, scRNA.select, method = "spearman")
write.table(corr, file = "Bulk_sc_WS_DEG_cor_0627.tsv", quote = F, sep = "\t",row.names = T,col.names=T)
colnames(corr) <- gsub("^[^_]*_", "", colnames(corr))
f1 <- colorRamp2(breaks = c(min(corr),0,0.8,max(corr)),c("#74ADD1","#FFFFBF","burlywood1","#D73027"))
corr <- corr[,c("24h","28h","48h","globular","heart","torpedo","bent","cotyledon","late_cotyledon")]
# matrix.corr <- matrix(corr)
#pheatmap
s2 <- pheatmap(corr, border_color = "black",
               color = f1,
               cluster_row = FALSE,
               cluster_cols = FALSE)

pdf(paste0("SCRNA_cor/bulk_Whole_Seed_cor/SC_bulk_WS_cor.pdf"),9,4)
print(s2)
dev.off()


#Heatmap
s2 <- Heatmap(as.matrix(corr), name = "spearman corr", 
              column_title = " ",row_title = " ",
              col = f1,
              #border = T,
              cluster_rows = FALSE, cluster_columns = FALSE,
              show_column_names = TRUE,show_row_names = TRUE,
              row_names_side = "right",
              row_dend_width = unit(10, "mm"),column_dend_height = unit(10,"mm"),
              column_dend_side = "top",column_title_rot = 0,
              column_title_side = "top",
              column_names_rot = 45)+pheatmap(corr, border_color = "black")
