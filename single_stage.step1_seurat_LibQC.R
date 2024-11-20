#!/usr/bin/Rscript
args=commandArgs(T)

#.libPaths("./R/library")

### Seurat clustering
library(Seurat)
library(dplyr)
library(patchwork)
library(DoubletFinder)
library(ggplot2)
library(future)
library(data.table)
library(COSG)
ntop<-3
plan()
plan("multiprocess", workers = as.numeric(args[1]))
plan()
options(future.globals.maxSize= 429496729600)

SampleID <- args[2] #R1087_Root_AT_0316
path <- args[3] 
minUMIs <- as.numeric(args[4]) #1000
minGenes <- as.numeric(args[5]) #500
maxPercent.mt <- as.numeric(args[6]) #4
maxPercent.C <- as.numeric(args[7]) #10
dim.usage <- as.numeric(args[8]) #20
res.usage <- as.numeric(args[9]) #0.5
doublets.percentage <- as.numeric(args[10]) #0.05

count <- Read10X(paste0(path,"/",SampleID,"/04.Matrix/"),gene.column=1)
colnames(count) <- paste(SampleID,colnames(count),sep=":")

### 2. clustering
dir.create(args[11], recursive = TRUE)
setwd(args[11]) 
dir.create(SampleID, recursive = TRUE)
setwd(SampleID)
### Creat Seurat object
scRNA = CreateSeuratObject(counts = count, min.cells = 1, min.features=1)

cat("1.raw matrix:",dim(scRNA),"\n")

scRNA@meta.data$Batch <- SampleID
cat("Batch:\n")
table(scRNA@meta.data$Batch)

nReads_stat <- as.data.frame(fread(paste0(path,"/",SampleID,"/merge_cell.stat")),header=T,sep="\t")
nReads_stat$Cell <- paste(SampleID,nReads_stat$Cell,sep=":")
if(all(colnames(scRNA)==nReads_stat$Cell)){
				scRNA@meta.data$nReads <- nReads_stat$Raw
}else{
				rownames(nReads_stat) <- nReads_stat$Cell
				nReads_stat <- nReads_stat[colnames(scRNA),]
				scRNA@meta.data$nReads <- nReads_stat$Raw
}

### for Arabidopsis_thaliana
geneinfo <- as.data.frame(fread("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/script/Araport11.Mar92021.geneinfo_20220708.txt"),header=T,sep="\t")

Mtgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"]) %in% rownames(scRNA))]
cat("Mitochondrial gene:", length(Mtgene),"\n", Mtgene,"\n")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, features = Mtgene) 

Cgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"]) %in% rownames(scRNA))]
cat("Chloroplast gene:",length(Cgene),"\n", Cgene,"\n")
scRNA[["percent.C"]] <- PercentageFeatureSet(scRNA, features = Cgene)

### Plot QC plot
features <- c("nReads", "nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
pdf(paste0(SampleID,"_QC_VlnPlot_raw.pdf"), length(features)*4.5, 5)
VlnPlot(scRNA, group.by="Batch", pt.size=0, features = features, ncol = length(features))
dev.off()

pdf(paste0(SampleID,"_QC_Scatter_raw.pdf"), 7, 5)
for (i in features){
				if (i %in% c("nCount_RNA")){
								cat("This is nCount_RNA.\n")
				}else{
								p <- FeatureScatter(scRNA, group.by="Batch", feature1 = "nCount_RNA", feature2 = i)
								print(p)
						 }
}
dev.off()

qc_stat <- rbind(summary(scRNA@meta.data$nReads),summary(scRNA@meta.data$nCount_RNA),summary(scRNA@meta.data$nFeature_RNA),summary(scRNA@meta.data$percent.mt),summary(scRNA@meta.data$percent.C))
rownames(qc_stat) <- c("nReads","nUMIs","nGenes","percent.mt","percent.C")
write.table(qc_stat, file = paste0(SampleID,"_qc_stat_raw.tsv"), quote = F, sep = "\t",row.names = T,col.names=NA)


### filter low quality cell
scRNA <- subset(scRNA, subset = nCount_RNA>minUMIs & nFeature_RNA > minGenes & percent.mt < maxPercent.mt & percent.C < maxPercent.C)
cat("2.Filter Low Quality Cell matrix:",dim(scRNA),"\n")

scRNA <- NormalizeData(object = scRNA, verbose = FALSE)
scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = 3000, verbose = FALSE)
scRNA <- ScaleData(scRNA) 
scRNA <- RunPCA(scRNA, verbose=F)
scRNA <- RunUMAP(scRNA, dims = 1:dim.usage)

### Define Find_doublet function
Find_doublet <- function(data){
	sweep.res.list <- paramSweep_v3(data, PCs = 1:dim.usage, sct = FALSE)
	sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
	bcmvn <- find.pK(sweep.stats) ### output plot
	nExp_poi <- round(doublets.percentage*ncol(data))
	p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC==max(bcmvn$MeanBC),]$pK)) ### pK Selection
	data <- doubletFinder_v3(data, PCs = 1:dim.usage, pN = 0.25, pK = p, nExp = nExp_poi, reuse.pANN = FALSE, sct = FALSE) ### output plot
	colnames(data@meta.data)[ncol(data@meta.data)] = "doublet_info"
	return(data)
}

### Find and filter doublets
pdf(paste0(SampleID,"_find_doublet.pdf"),7,7)
scRNA <- Find_doublet(scRNA)
dev.off()
#write.table(scRNA@meta.data,paste0(SampleID,"_doublets_info.txt"),quote = F, sep = "\t",row.names = T,col.names=NA)

scRNA <- subset(scRNA,subset=doublet_info=="Singlet")

cat("3.FilterDoublet matrix:",dim(scRNA),"\n")
 
pdf(paste0(SampleID,"_QC_VlnPlot_filter_cell_v1.pdf"), length(features)*4.5, 5)
p <- VlnPlot(scRNA, group.by="Batch", pt.size=0, features = features, ncol = length(features))
print(p)
dev.off()

chr_gene_rmM<-as.vector(geneinfo[which(geneinfo$Chr !="ChrM") ,"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr !="ChrM") ,"GeneID"]) %in% rownames(scRNA))]
chr_gene_rmC<-as.vector(geneinfo[which(geneinfo$Chr !="ChrC") ,"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr !="ChrC") ,"GeneID"]) %in% rownames(scRNA))]
chr_gene<-intersect(chr_gene_rmM,chr_gene_rmC)

saveRDS(scRNA, file = paste0(SampleID,"_scRNA_allChr.rds"))
scRNA<-subset(scRNA,features = chr_gene)
saveRDS(scRNA, file = paste0(SampleID,"_scRNA_rmChrMC.rds"))

### Find cluster
scRNA <- FindNeighbors(scRNA, dims = 1:dim.usage)
scRNA <- FindClusters(scRNA, resolution = res.usage)

saveRDS(scRNA, file = paste0(SampleID,"_scRNA.rds"))



### Cluster_Batch
pdf(paste0(SampleID,"_umap_cluster.pdf"),8,7)
DimPlot(object = scRNA, reduction = "umap",pt.size = 0.1,label=T)+ggtitle(label = "umap_Cluster")
DimPlot(object = scRNA, reduction = "umap",group.by = "Batch",pt.size = 0.1,label=F)+ggtitle(label = "umap_Batch")
dev.off()

### umi_gene
features <- c("nReads","nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
pdf(paste0(SampleID,"_umap_gene_umi.pdf"),length(features)*4.5,5)
FeaturePlot(scRNA, order = T, features = features, ncol = length(features), reduction = "umap")
FeaturePlot(scRNA, order = T, features = features, ncol = length(features), reduction = "umap",max.cutoff = "q95")
dev.off()

### count matrix
write.table(as.matrix(scRNA@assays$RNA@counts),paste0(SampleID,"_count_mat_FilterDoublet.txt"),sep="\t",quote = FALSE,row.names = T,col.names=NA)

### DEG
#DEG <- FindAllMarkers(scRNA, only.pos = FALSE, verbose = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
#DEG_order <- DEG[order(DEG$cluster,-DEG$avg_log2FC),]
#DEG_order$avg_FC <- 2^DEG_order$avg_log2FC
#write.table(DEG_order,paste0(SampleID,"_DEG.txt"),sep="\t",quote = FALSE,row.names = T,col.names=NA)

## cosg
DEG_cosg <- cosg(scRNA, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=200)
DEG_cosg_merge <- data.frame(cluster=rep(colnames(DEG_cosg$names),each=200), gene=as.vector(as.matrix(DEG_cosg$names)), score=as.vector(as.matrix(DEG_cosg$scores)))
write.table(DEG_cosg_merge,paste0(SampleID,"_DEG_cosg.tsv"),sep="\t",quote = FALSE,row.names=FALSE)

top_list <- c()
for (group in colnames(DEG_cosg$names)){
	top_i <- DEG_cosg$names[group][1:ntop,1] #ntop
	top_list <- c(top_list,top_i)
}
      
pdf(paste0(SampleID,"_DEG_cosg.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
print(p)
dev.off()

## wilcox
DEG_wilcox <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
DEG_wilcox_order <- DEG_wilcox[order(DEG_wilcox$cluster,-DEG_wilcox$avg_log2FC),]
DEG_wilcox_order$avg_FC <- 2^DEG_wilcox_order$avg_log2FC
write.table(DEG_wilcox_order,paste0(SampleID,"_DEG_wilcox.tsv"),sep="\t",row.names=F,quote = FALSE)

DEG_wilcox_order %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC) -> top_list
pdf(paste0(SampleID,"_DEG_wilcox.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list$gene)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
print(p)
dev.off()
