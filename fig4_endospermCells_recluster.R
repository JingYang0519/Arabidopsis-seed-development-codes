library(Seurat)
library(Matrix)
library(dplyr)

rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202403.CellFlow/stage2-9.rm36LowQualityCells.VeryDetaiAnno.rds")

rds <- subset(rds, Organ != "9_late_cotyledon" & CellType == "Endosperm")
rds

# 重新计算高变等
rds <- NormalizeData(rds, verbose = FALSE)
rds <- FindVariableFeatures(rds, selection.method = "vst", nfeatures = 1500, verbose = FALSE)
rds <- ScaleData(rds)
rds <- RunPCA(rds, verbose = FALSE)
rds <- FindNeighbors(rds, dims = 1:10)
rds <- FindClusters(rds, resolution = 1.0)
rds <- RunUMAP(rds, dims = 1:20)

levels(rds)

scRNA = rds
names(scRNA@meta.data)

scRNA
names(scRNA@meta.data)
unique(scRNA$Organ)
unique(scRNA$DetailAnno)
unique(scRNA$Organ_Cluster.subtype)

### cellcoor
cluster_ID = as.data.frame(Idents(object = scRNA))
cluster_cor1 = as.data.frame(Embeddings(object = scRNA, reduction = "umap"))
res = cbind(cluster_ID, cluster_cor1)
write.table(res,"Endo_Cells_recluster.stage2-8.rm36LowQualityCells.cellcoor.tsv",sep="\t",quote = FALSE)

library(ggplot2)

pdf("EndoCells_recluster_UMAP.pdf")
p1 <- DimPlot(scRNA, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
print(p1)

p2 <- DimPlot(scRNA, reduction = "umap", group.by = "Organ", label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
print(p2)

p3 <- DimPlot(scRNA, reduction = "umap", group.by = "Batch", label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
print(p3)

p4 <- DimPlot(scRNA, reduction = "umap", group.by = "DetailAnno", label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
print(p4)

p5 <- DimPlot(scRNA, reduction = "umap", group.by = "Organ_Cluster.subtype", label = TRUE, repel = TRUE) + theme(aspect.ratio = 1)
print(p5)

dev.off()

### DEG
## cosg

DEG_cosg <- cosg(scRNA, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=200)

DEG_cosg_merge <- data.frame(cluster=rep(colnames(DEG_cosg$names), each=200), gene=as.vector(as.matrix(DEG_cosg$names)), score=as.vector(as.matrix(DEG_cosg$scores)))

write.table(DEG_cosg_merge, "Endo_Cells_recluster.stage2-8.rm36LowQualityCells.DEG_cosg.tsv", sep="\t", quote = FALSE, row.names=FALSE)


top_list <- c()
ntop <- as.numeric(3)

for (group in colnames(DEG_cosg$names)) {
    top_i <- DEG_cosg$names[group][1:ntop, 1] 
    top_list <- c(top_list, top_i)
}


pdf("Endo_Cells_recluster.stage2-8.rm36LowQualityCells.DEG_cosg.pdf", length(table(Idents(scRNA))) * 0.5, length(table(Idents(scRNA))) * 0.2 * ntop)

p <- DotPlot(scRNA, assay = 'RNA', features = unique(top_list)) + 
     RotatedAxis() + coord_flip() + 
     theme(axis.title.x=element_blank(), axis.title.y=element_blank())
print(p)
dev.off()

## wilcox

DEG_wilcox <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")

DEG_wilcox_order <- DEG_wilcox[order(DEG_wilcox$cluster, -DEG_wilcox$avg_log2FC),]

DEG_wilcox_order$avg_FC <- 2^DEG_wilcox_order$avg_log2FC

write.table(DEG_wilcox_order,"Endo_Cells_recluster.stage2-8.rm36LowQualityCells.DEG_wilcox.tsv", sep="\t", row.names=F, quote = FALSE)


DEG_wilcox_order %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC) -> top_list

pdf("Endo_Cells_recluster.stage2-8.rm36LowQualityCells.DEG_wilcox.pdf", length(table(Idents(scRNA))) * 0.5, length(table(Idents(scRNA))) * 0.2 * ntop)

p <- DotPlot(scRNA, assay = 'RNA', features = unique(top_list$gene)) + 
     RotatedAxis() + coord_flip() + 
     theme(axis.title.x=element_blank(), axis.title.y=element_blank())
print(p)
dev.off()

