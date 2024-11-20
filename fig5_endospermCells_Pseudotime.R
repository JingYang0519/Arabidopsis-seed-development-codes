library(monocle)
library(RColorBrewer)
library(ggplot2)
library(Seurat)
library(ggsci)
library(ggpubr)

rds <- readRDS("refined_anno0514.Endo_Cells_recluster.stage2-8.rm36LowQualityCells.rds")
rds

rds@active.ident <- factor(rds@meta.data$seurat_clusters_renamed)
names(rds@active.ident) <- rownames(rds@meta.data)

# Find diff genes betweeen cell types
dif_28h <- FindMarkers(rds, ident.1 = "2_28h_Endosperm", min.pct = 0.08)
dif_cotyledon <- FindMarkers(rds, ident.1 = "8_cotyledon_Endosperm", min.pct = 0.08)

#dif <- FindMarkers(rds, ident.1 = args[2], min.pct = 0.08)
write.table(dif_28h,file="2_28h_Endosperm_monocle_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")
#dif_28 <- rownames(dif_28h[which(dif_28h$p_val<0.05),])
dif_p<-dif_28h[which(dif_28h$p_val<0.05),]
dif_28 <- rownames(dif_p[which(dif_p$avg_log2FC>1),])
write.table(dif_cotyledon,file="8_cotyledon_Endosperm_monocle_DEG.txt",row.names=T,col.names=T,quote=F,sep="\t")
#dif_co <- rownames(dif_cotyledon[which(dif_cotyledon$p_val<0.05),])
dif_p<-dif_cotyledon[which(dif_cotyledon$p_val<0.05),]
dif_co <- rownames(dif_p[which(dif_p$avg_log2FC>1),])
dif_all<-union(dif_28,dif_co)
#dif<-setdiff(dif_all,stat2vsstat3[,7])
dif<-dif_all
length(dif)

umi <- data.frame(rds@assays$RNA@counts, check.names=FALSE)
metadata <- data.frame(rds@meta.data)

write.table(umi, file = "stage1-9.endo_umi.txt", quote = FALSE, sep = "\t")
write.table(metadata, file = "stage1-9.endo_metadata.txt", quote = FALSE, sep = "\t")

data <- umi

type <- metadata

expr_matrix <- data[,rownames(type)]   # expression matrix
gene_annotation <- data.frame(gene_id=rownames(expr_matrix),
                              gene_short_name=rownames(expr_matrix))
rownames(gene_annotation) <- gene_annotation[,1]   # gene annotation
sample_sheet <- type   # corresponding cell annotation information

pos <- which(rownames(gene_annotation) %in% rownames(expr_matrix))
gene_annotation <- gene_annotation[pos,]
expr_matrix <- expr_matrix[rownames(gene_annotation), rownames(sample_sheet)]
pd <- new("AnnotatedDataFrame", data = sample_sheet)   # cell annotation object
fd <- new("AnnotatedDataFrame", data = gene_annotation)   # gene annotation object

# UMI (default input): for negative binomal distribution data
cd <- newCellDataSet(as(as.matrix(expr_matrix), "sparseMatrix"), phenoData = pd,
                     featureData = fd, expressionFamily = negbinomial.size(),
                     lowerDetectionLimit = 0.5)

cd <- estimateSizeFactors(cd)   # help eliminate differences in mRNA capture between cells
cd <- estimateDispersions(cd)   # for subsequent differential expression analysis
cd <- detectGenes(cd, min_expr = 0.5)   # counts the number of genes expressed in each cell

# if ("monocle3" %in% loadedNamespaces()) {
#   detach("package:monocle3", unload = TRUE)
# }

library(monocle)
feature_data <- fData(cd)
expressed_genes <- row.names(subset(fData(cd), num_cells_expressed > nrow(sample_sheet) * 0.01))
length(expressed_genes)   # genes expressed in at least 1% cells with expression > 0.5

ordering_genes <- dif
ordering_genes <- intersect(ordering_genes, expressed_genes)
cd <- setOrderingFilter(cd, ordering_genes)

# Order Cells by Progress
cd <- reduceDimension(cd, max_components = 2, method = 'DDRTree')
cd <- orderCells(cd, reverse = T)   # reverse = F by default
pseudotime_7stagesDEG = cd
save(pseudotime_7stagesDEG,file="endo_allcells.pseudotime_7stagesDEG.rda")

names(rds@meta.data)
rds@meta.data$cellID <- rownames(rds@meta.data)
names(rds@meta.data)

meta <- pData(cd)

same_order <- rownames(rds@meta.data) == rownames(meta)

if (all(same_order)) {
  rds@meta.data$pseudotime <- meta$Pseudotime
} else {

  merged_meta <- merge(rds@meta.data, meta, by.x="cellID", by.y="rownames", all.x=TRUE)
  rds@meta.data <- merged_meta
}

colnames(rds@meta.data)

saveRDS(rds, file="stage1-9.endo_rmCZE_pseudotime.rds")

write.table(rds@meta.data, file = "stage1-9.endo_rmCZE_pseudotime_metadata.txt", quote = FALSE, sep = "\t")
