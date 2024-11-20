#!/usr/bin/Rscript
args=commandArgs(T)

#.libPaths("./R/library")

### Seurat clustering
library(Seurat)
library(harmony)
library(COSG)
library(dplyr)
library(data.table)
library(patchwork)
library(ggplot2)
library(future)
library(PNWColors)

plan()
plan("multiprocess", workers = as.numeric(args[1])) #3
plan()
options(future.globals.maxSize= 429496729600)

Organ <- args[2] #Root
dim.usage <- as.numeric(args[3]) #20
res.usage <- as.numeric(args[4]) #0.8
minUMIs <- as.numeric(args[5]) #1000
minGenes <- as.numeric(args[6]) #500
maxPercent.mt <- as.numeric(args[7]) #4
maxPercent.C <- as.numeric(args[8]) #10
ToIntegrate <- args[9] #TRUE
ntop <- as.numeric(args[10]) #3 ntop: n top marker gene to plot
hvg <- as.numeric(args[11]) #3000
#RunTime <- args[12]
LibQC<-args[12]
OutDir<-args[13]
InDir<-args[14]
### 1. load data

setwd(LibQC)
file_list = list.files(pattern = "count_mat_FilterDoublet.txt",recursive=T,full.name=T) ## select Sample
cat(Organ,"has",length(file_list),"samples.\n")

if (length(file_list) > 1){
  cat(Organ,"has",length(file_list),"qualified samples to merge.\n")

  count = list()
  for(i in file_list){
		cat(i,"\n")
    count[[i]] <- as.data.frame(fread(i))
    colnames(count[[i]])[1] <- "ID"
  }

  counts <- Reduce(function(dtf1,dtf2) full_join(dtf1,dtf2,by="ID"),count)
  rownames(counts) = counts[,which(colnames(counts) %in% "ID")]
  counts[is.na(counts)] <- 0
  counts = counts[, -which(colnames(counts) %in% "ID")]
  cat(Organ,"raw count matrix:",dim(counts),"\n")

  ### 2. clustering
	dir.create(OutDir, recursive = TRUE)
	setwd(OutDir)
  ### Creat Seurat object
  scRNA = CreateSeuratObject(counts = counts, min.cells = 1, min.features=1)
  ## Warning: Feature names cannot have underscores ('_'), replacing with dashes ('-')
  scRNA@meta.data$Organ <- Organ
  scRNA@meta.data$Batch <- unlist(lapply(strsplit(as.character(colnames(scRNA)),":"),"[",1))
  cat("\n Organ")
  print(table(scRNA@meta.data$Organ))
  cat("\n Batch")
  print(table(scRNA@meta.data$Batch))

### input nReads_stat
nReads_stat <- list()
for(i in unique(scRNA@meta.data$Batch)){
#				nReads_stat[[i]] <- as.data.frame(fread(paste0("scRNA/embryo_72h/",i,"/merge_cell.stat")),header=T,sep="\t")
	nReads_stat[[i]] <- as.data.frame(fread(paste0(InDir,"/",i,"/merge_cell.stat")),header=T,sep="\t")
	nReads_stat[[i]]$Cell <- paste(i,nReads_stat[[i]]$Cell,sep=":")
}
### merge nReads_stat
nReads_stat_merge <- Reduce(function(dtf1,dtf2) rbind(dtf1,dtf2),nReads_stat)
### add nReads info to meta.data
if(all(colnames(scRNA)==nReads_stat_merge$Cell)){
				scRNA@meta.data$nReads <- nReads_stat_merge$Raw
}else{
				rownames(nReads_stat_merge) <- nReads_stat_merge$Cell
				nReads_stat_merge <- nReads_stat_merge[colnames(scRNA),]
				scRNA@meta.data$nReads <- nReads_stat_merge$Raw
}
 
 ### for Arabidopsis_thaliana
geneinfo <- as.data.frame(fread("Araport11.Mar92021.geneinfo_20220708.txt"),header=T,sep="\t")

Mtgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrM"),"GeneID"]) %in% rownames(scRNA))]
cat("Mitochondrial gene:", length(Mtgene),"\n", Mtgene,"\n")
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, features = Mtgene) 

Cgene <- as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr=="ChrC"),"GeneID"]) %in% rownames(scRNA))]
cat("Chloroplast gene:",length(Cgene),"\n", Cgene,"\n")
scRNA[["percent.C"]] <- PercentageFeatureSet(scRNA, features = Cgene)


  ### Plot QC plot
  features <- c("nReads", "nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
  pdf(paste0(Organ,"_QC_VlnPlot.pdf"), length(features)*4.5, 5)
  p <- VlnPlot(scRNA, group.by="Organ", pt.size=0, features = features, ncol = length(features))
  print(p)
	p <- VlnPlot(scRNA, group.by="Batch", pt.size=0, features = features, ncol = length(features))
  print(p)
	dev.off()

  pdf(paste0(Organ,"_QC_Scatter.pdf"), 7, 5)
	for (i in features){
					if (i %in% c("nCount_RNA")){
									cat("This is nCount_RNA.\n")
					}else{
									p <- FeatureScatter(scRNA, group.by="Organ", feature1 = "nCount_RNA", feature2 = i)
									print(p)
									p <- FeatureScatter(scRNA, group.by="Batch", feature1 = "nCount_RNA", feature2 = i)
									print(p)
							 }
	}
  dev.off()

  ### QC
  qc_stat <- rbind(summary(scRNA@meta.data$nReads),summary(scRNA@meta.data$nCount_RNA),summary(scRNA@meta.data$nFeature_RNA),summary(scRNA@meta.data$percent.mt),summary(scRNA@meta.data$percent.C))
  rownames(qc_stat) <- c("nReads","nUMIs","nGenes","percent.mt","percent.C")
  write.table(qc_stat, file = paste0(Organ,"_qc_stat.tsv"), quote = F, sep = "\t",row.names = T,col.names=NA)
  
  scRNA <- subset(scRNA, subset = nCount_RNA > minUMIs & nFeature_RNA > minGenes & percent.mt < maxPercent.mt & percent.C < maxPercent.C)

	saveRDS(scRNA, file = paste0(Organ,"_scRNA_allChr.rds"))
	
	chr_gene_rmM<-as.vector(geneinfo[which(geneinfo$Chr !="ChrM") ,"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr !="ChrM") ,"GeneID"]) %in% rownames(scRNA))]
	chr_gene_rmC<-as.vector(geneinfo[which(geneinfo$Chr !="ChrC") ,"GeneID"])[which(as.vector(geneinfo[which(geneinfo$Chr !="ChrC") ,"GeneID"]) %in% rownames(scRNA))]
	chr_gene<-intersect(chr_gene_rmM,chr_gene_rmC)

	scRNA<-subset(scRNA,features = chr_gene)
	saveRDS(scRNA, file = paste0(Organ,"_scRNA_rmChrMC.rds"))

  ### Plot filter_cell QC plot
  features <- c("nReads", "nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
  pdf(paste0(Organ,"_QC_VlnPlot_filter_cell.pdf"), length(features)*4.5, 5)
  p <- VlnPlot(scRNA, group.by="Organ", pt.size=0.01, features = features, ncol = length(features))
  print(p)
  p <- VlnPlot(scRNA, group.by="Batch", pt.size=0.01, features = features, ncol = length(features))
  print(p)
  dev.off()
	pdf(paste0(Organ,"_QC_VlnPlot_filter_cell_v1.pdf"), length(features)*4.5, 5)
	p <- VlnPlot(scRNA, group.by="Organ", pt.size=0, features = features, ncol = length(features))
	
print(p)
	p <- VlnPlot(scRNA, group.by="Batch", pt.size=0, features = features, ncol = length(features))
	print(p)
	dev.off()
  ### QC
  qc_stat <- rbind(summary(scRNA@meta.data$nReads),summary(scRNA@meta.data$nCount_RNA),summary(scRNA@meta.data$nFeature_RNA),summary(scRNA@meta.data$percent.mt),summary(scRNA@meta.data$percent.C))
  rownames(qc_stat) <- c("nReads","nUMIs","nGenes","percent.mt","percent.C")
  write.table(qc_stat, file = paste0(Organ,"_qc_stat_filter_cell.tsv"), quote = F, sep = "\t",row.names = T,col.names=NA)

################################################## 1. merge ##################################################

  scRNA <- NormalizeData(object = scRNA, verbose = FALSE)
  scRNA <- FindVariableFeatures(object = scRNA, selection.method = "vst", nfeatures = hvg, verbose = FALSE)
  scRNA <- ScaleData(scRNA) 
  scRNA <- RunPCA(scRNA,verbose=F)

  ### Find cluster
  scRNA <- FindNeighbors(scRNA, dims = 1:dim.usage)
  scRNA <- FindClusters(scRNA, resolution = res.usage)
  scRNA <- RunUMAP(scRNA, dims = 1:dim.usage)

  saveRDS(scRNA, file = paste0(Organ,"_scRNA.rds"))

  ### Cluster_Batch
  pdf(paste0(Organ,"_Cluster.pdf"),8,7)
  col <- pnw_palette("Sailboat",length(table(Idents(scRNA))),type="continuous")
	p <- DimPlot(object = scRNA, reduction = "umap",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_Cluster")
  print(p)
  #p <- DimPlot(object = scRNA, reduction = "umap",group.by = "Organ",pt.size = 0.1,label=F)+ggtitle(label = "umap_Organ")
  #print(p)
	#col <- pnw_palette("Bay",length(table(scRNA@meta.data$Batch)),type="discrete")
	col <- pnw_palette("Bay",length(table(scRNA@meta.data$Batch)),type="continuous")
  p <- DimPlot(object = scRNA, reduction = "umap",group.by = "Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_Batch")
  print(p)
  dev.off()

  pdf(paste0(Organ,"_Cluster_BatchSplit.pdf"),length(unique(scRNA@meta.data$Batch))*5+1,6)
	col <- pnw_palette("Sailboat",length(table(Idents(scRNA))),type="continuous")
  p <- DimPlot(scRNA, reduction = "umap", split.by = "Batch",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_BatchSplit")
  print(p)
	col <- pnw_palette("Sailboat",length(unique(scRNA@meta.data$Batch)),type="continuous")
  p <- DimPlot(scRNA, reduction = "umap", split.by = "Batch",group.by="Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_BatchSplit")
  print(p)
  dev.off()

  ### umi_gene_mt_C
  features <- c("nReads","nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
  pdf(paste0(Organ,"_UMAP_nReads_gene_umi_mt_C.pdf"),length(features)*5+1,5)
  p <- FeaturePlot(scRNA, order = T, features = features, reduction = "umap", ncol = length(features))
  print(p)
  p <- FeaturePlot(scRNA, order = T, features = features, reduction = "umap", max.cutoff = "q95", ncol = length(features))
  print(p)
  dev.off()

  ### marker feature plot
  #marker <- c('AT1G22710', 'AT3G25710', 'AT1G79580', 'AT5G49270')
  #pdf(paste0(Organ,"_UMAP_marker.pdf"),length(marker)*5+1,5)
  #p <- FeaturePlot(scRNA, features = marker, order = T, reduction = "umap", ncol = length(marker))
  #print(p)
  #p <- FeaturePlot(scRNA, features = marker, order = T, reduction = "umap", max.cutoff = "q95", ncol = length(marker))
  #print(p)
  #dev.off()

  ### meta data
  write.table(scRNA@meta.data,file = paste0(Organ,"_meta_data.tsv"),sep="\t",quote=FALSE)

  ### cellcoor
  cluster_ID = as.data.frame(Idents(object = scRNA))
  cluster_cor1 = as.data.frame(Embeddings(object = scRNA,reduction = "umap"))
  res = cbind(cluster_ID,cluster_cor1)
  write.table(res,paste0(Organ,"_cellcoor.tsv"),sep="\t",quote = FALSE)

  ### DEG
  ## cosg
  DEG_cosg <- cosg(scRNA, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=200)
  DEG_cosg_merge <- data.frame(cluster=rep(colnames(DEG_cosg$names),each=200), gene=as.vector(as.matrix(DEG_cosg$names)), score=as.vector(as.matrix(DEG_cosg$scores)))
  write.table(DEG_cosg_merge,paste0(Organ,"_DEG_cosg.tsv"),sep="\t",quote = FALSE,row.names=FALSE)

  top_list <- c()
  for (group in colnames(DEG_cosg$names)){
      top_i <- DEG_cosg$names[group][1:ntop,1] #ntop
      top_list <- c(top_list,top_i)
  }
      
  pdf(paste0(Organ,"_DEG_cosg.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
  p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

  ## wilcox
  DEG_wilcox <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
  DEG_wilcox_order <- DEG_wilcox[order(DEG_wilcox$cluster,-DEG_wilcox$avg_log2FC),]
  DEG_wilcox_order$avg_FC <- 2^DEG_wilcox_order$avg_log2FC
  write.table(DEG_wilcox_order,paste0(Organ,"_DEG_wilcox.tsv"),sep="\t",row.names=F,quote = FALSE)

  DEG_wilcox_order %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC) -> top_list
  pdf(paste0(Organ,"_DEG_wilcox.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
  p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list$gene)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

################################################## 2. integrated ##################################################
if (ToIntegrate){
  if (length(unique(scRNA@meta.data$Batch)) > 1){
    cat(Organ,"has",length(file_list),"(",file_list,")","qualified samples from",length(unique(scRNA@meta.data$Batch)),"(",unique(scRNA@meta.data$Batch),") Batchs to integrate.\n")
    
    ############################################## 2.1 Harmony
    cat("Run Harmony.\n")
    ##############################################
    p1 <- DimPlot(object = scRNA, reduction = "pca", pt.size = .1, group.by = "Batch")
    p2 <- VlnPlot(object = scRNA, features = "PC_1", group.by = "Batch", pt.size = .1)
    pdf(paste0(Organ,"_PCA.pdf"),12,5)
    print(p1+p2)
    dev.off()

    scRNA <- RunHarmony(scRNA, group.by.vars="Batch", verbose=F, plot_convergence = FALSE)

    p1 <- DimPlot(object = scRNA, reduction = "harmony", pt.size = .1, group.by = "Batch")
    p2 <- VlnPlot(object = scRNA, features = "harmony_1", group.by = "Batch", pt.size = .1)
    pdf(paste0(Organ,"_Harmony.pdf"),12,5)
    print(p1+p2)
    dev.off()

    scRNA <- scRNA %>%
    RunUMAP(reduction = "harmony", dims = 1:dim.usage) %>%
    FindNeighbors(reduction = "harmony", dims = 1:dim.usage) %>%
    FindClusters(resolution = res.usage) %>%
    identity()

    saveRDS(scRNA, file = paste0(Organ,"_scRNA.harmony.rds"))

  ### Cluster_Batch
	#col <- pnw_palette("Bay",length(table(Idents(scRNA))),type="continuous")
	#col <- pnw_palette("Sunset2",length(table(Idents(scRNA))),type="continuous")
	#col <- met.brewer(name="VanGogh1",n=length(table(Idents(scRNA))),type="continuous")
	pdf(paste0(Organ,"_Cluster_harmony.pdf"),8,7)
	col <- pnw_palette("Sailboat",length(table(Idents(scRNA))),type="continuous")
  p <- DimPlot(object = scRNA, reduction = "umap",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_Cluster")
  print(p)
  #p <- DimPlot(object = scRNA, reduction = "umap",group.by = "Organ",pt.size = 0.1,label=F)+ggtitle(label = "umap_Organ")
  #print(p)
	#col <- pnw_palette("Bay",length(table(scRNA@meta.data$Batch)),type="discrete")
	col <- pnw_palette("Bay",length(table(scRNA@meta.data$Batch)),type="continuous")
  p <- DimPlot(object = scRNA, reduction = "umap",group.by = "Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_Batch")
  print(p)
  dev.off()

  pdf(paste0(Organ,"_Cluster_BatchSplit_harmony.pdf"),length(unique(scRNA@meta.data$Batch))*5+1,6)
  col <- pnw_palette("Sailboat",length(table(Idents(scRNA))),type="continuous")
	p <- DimPlot(scRNA, reduction = "umap", split.by = "Batch",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_BatchSplit")
  print(p)
	col <- pnw_palette("Sailboat",length(unique(scRNA@meta.data$Batch)),type="continuous")
	p <- DimPlot(scRNA, reduction = "umap", split.by = "Batch",group.by="Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_BatchSplit")
  print(p)
  dev.off()

  ### umi_gene_mt_C
  features <- c("nReads","nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
  pdf(paste0(Organ,"_UMAP_nReads_gene_umi_mt_C_harmony.pdf"),length(features)*5+1,5)
  p <- FeaturePlot(scRNA, order = T, features = features, reduction = "umap", ncol = length(features))
  print(p)
  p <- FeaturePlot(scRNA, order = T, features = features, reduction = "umap", max.cutoff = "q95", ncol = length(features))
  print(p)
  dev.off()

  ### meta data
  write.table(scRNA@meta.data,file = paste0(Organ,"_meta_data_harmony.tsv"),sep="\t",quote=FALSE)

  ### cellcoor
  cluster_ID = as.data.frame(Idents(object = scRNA))
  cluster_cor1 = as.data.frame(Embeddings(object = scRNA,reduction = "umap"))
  res = cbind(cluster_ID,cluster_cor1)
  write.table(res,paste0(Organ,"_cellcoor_harmony.tsv"),sep="\t",quote = FALSE)

  ### DEG
  ## cosg
  DEG_cosg <- cosg(scRNA, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=200)
  DEG_cosg_merge <- data.frame(cluster=rep(colnames(DEG_cosg$names),each=200), gene=as.vector(as.matrix(DEG_cosg$names)), score=as.vector(as.matrix(DEG_cosg$scores)))
  write.table(DEG_cosg_merge,paste0(Organ,"_DEG_cosg_harmony.tsv"),sep="\t",quote = FALSE,row.names=FALSE)

  top_list <- c()
  for (group in colnames(DEG_cosg$names)){
      top_i <- DEG_cosg$names[group][1:ntop,1] #ntop
      top_list <- c(top_list,top_i)
  }
      
  pdf(paste0(Organ,"_DEG_cosg_harmony.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
  p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

  ## wilcox
  DEG_wilcox <- FindAllMarkers(scRNA, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
  DEG_wilcox_order <- DEG_wilcox[order(DEG_wilcox$cluster,-DEG_wilcox$avg_log2FC),]
  DEG_wilcox_order$avg_FC <- 2^DEG_wilcox_order$avg_log2FC
  write.table(DEG_wilcox_order,paste0(Organ,"_DEG_wilcox_harmony.tsv"),sep="\t",row.name=F,quote = FALSE)

  DEG_wilcox_order %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC) -> top_list
  pdf(paste0(Organ,"_DEG_wilcox_harmony.pdf"),length(table(Idents(scRNA)))*0.5,length(table(Idents(scRNA)))*0.2*ntop)
  p <- DotPlot(scRNA, assay = 'RNA', features =  unique(top_list$gene)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

    ############################################## 2.2 Seurat integrate
    cat("Run Seurat integrate.\n")
    ##############################################
    scRNA.list <- SplitObject(object = scRNA, split.by = "Batch")

    for (i in 1:length(x = scRNA.list)) {
        scRNA.list[[i]] <- NormalizeData(object = scRNA.list[[i]], verbose = FALSE)
        scRNA.list[[i]] <- FindVariableFeatures(object = scRNA.list[[i]], selection.method = "vst", nfeatures = hvg, verbose = FALSE)
    }

    reference.list <- scRNA.list
    scRNA.anchors <- FindIntegrationAnchors(object.list = reference.list, dims = 1:dim.usage)
    scRNA.integrated <- IntegrateData(anchorset = scRNA.anchors, dims = 1:dim.usage)

    DefaultAssay(object = scRNA.integrated) <- "integrated"
    scRNA.integrated <- ScaleData(object = scRNA.integrated, verbose = FALSE)
    scRNA.integrated <- RunPCA(object = scRNA.integrated, npcs = dim.usage, verbose = FALSE)
    scRNA.integrated <- RunUMAP(object = scRNA.integrated, reduction = "pca", dims = 1:dim.usage)

    scRNA.integrated <- FindNeighbors(object = scRNA.integrated,k.param=40,dims = 1:dim.usage)
    scRNA.integrated.Find <- FindClusters(object = scRNA.integrated,resolution = res.usage)
    scRNA.integrated.Find <- RunTSNE(object = scRNA.integrated.Find,dims = 1:dim.usage,seed.use = 2)
    DefaultAssay(scRNA.integrated.Find)="RNA"

    saveRDS(scRNA.integrated.Find, file = paste0(Organ,"_scRNA.integrated.Find.rds"))
      
    ### Cluster_Batch
    pdf(paste0(Organ,"_Cluster_integrated.pdf"),8,7)
		col <- pnw_palette("Sailboat",length(table(Idents(scRNA.integrated.Find))),type="continuous")
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_Cluster")
    print(p)
    #p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",group.by = "Organ",pt.size = 0.1,label=F)+ggtitle(label = "umap_Organ")
    #print(p)
		col <- pnw_palette("Bay",length(table(scRNA.integrated.Find@meta.data$Batch)),type="continuous")
		#col <- pnw_palette("Bay",length(table(scRNA.integrated.Find@meta.data$Batch)),type="discrete")
    p <- DimPlot(object = scRNA.integrated.Find, reduction = "umap",group.by = "Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_Batch")
    print(p)
    dev.off()

    pdf(paste0(Organ,"_Cluster_BatchSplit_integrated.pdf"),length(unique(scRNA@meta.data$Batch))*5+1,6)
		col <- pnw_palette("Sailboat",length(table(Idents(scRNA.integrated.Find))),type="continuous")
    p <- DimPlot(scRNA.integrated.Find, reduction = "umap", split.by = "Batch",pt.size = 0.1,label=T,cols=col)+ggtitle(label = "umap_BatchSplit")
    print(p)
    p <- DimPlot(scRNA.integrated.Find, reduction = "umap", split.by = "Batch",pt.size = 0.1,label=F,cols=col)+ggtitle(label = "umap_BatchSplit")
    print(p)
    dev.off()

    ### umi_gene_mt_C
    features <- c("nReads","nFeature_RNA", "nCount_RNA", "percent.mt","percent.C")
    pdf(paste0(Organ,"_UMAP_nReads_gene_umi_mt_C_integrated.pdf"),length(features)*5+1,5)
    p <- FeaturePlot(scRNA.integrated.Find, order = T, features = features, reduction = "umap", ncol = length(features))
    print(p)
    p <- FeaturePlot(scRNA.integrated.Find, order = T, features = features, reduction = "umap", max.cutoff = "q95", ncol = length(features))
    print(p)
    dev.off()

    ### meta data
    write.table(scRNA.integrated.Find@meta.data,file = paste0(Organ,"_meta_data_integrated.tsv"),sep="\t",quote=FALSE)

    ### cellcoor
    cluster_ID = as.data.frame(Idents(object = scRNA.integrated.Find))
    cluster_cor1 = as.data.frame(Embeddings(object = scRNA.integrated.Find, reduction = "umap"))
    res = cbind(cluster_ID, cluster_cor1)
    write.table(res,paste0(Organ,"_cellcoor_integrated.tsv"),sep="\t",quote = FALSE)

    ### DEG
    ## cosg
  DEG_cosg <- cosg(scRNA.integrated.Find, groups='all', assay='RNA', slot='data', mu=1, n_genes_user=200)
  DEG_cosg_merge <- data.frame(cluster=rep(colnames(DEG_cosg$names),each=200), gene=as.vector(as.matrix(DEG_cosg$names)), score=as.vector(as.matrix(DEG_cosg$scores)))
  write.table(DEG_cosg_merge,paste0(Organ,"_DEG_cosg_integrated.tsv"),sep="\t",quote = FALSE,row.names=FALSE)

  top_list <- c()
  for (group in colnames(DEG_cosg$names)){
      top_i <- DEG_cosg$names[group][1:ntop,1] #ntop
      top_list <- c(top_list,top_i)
  }
      
  pdf(paste0(Organ,"_DEG_cosg_integrated.pdf"),length(table(Idents(scRNA.integrated.Find)))*0.5,length(table(Idents(scRNA.integrated.Find)))*0.2*ntop)
  p <- DotPlot(scRNA.integrated.Find, assay = 'RNA', features =  unique(top_list)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

  ## wilcox
  DEG_wilcox <- FindAllMarkers(scRNA.integrated.Find, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25, test.use = "wilcox")
  DEG_wilcox_order <- DEG_wilcox[order(DEG_wilcox$cluster,-DEG_wilcox$avg_log2FC),]
  DEG_wilcox_order$avg_FC <- 2^DEG_wilcox_order$avg_log2FC
  write.table(DEG_wilcox_order,paste0(Organ,"_DEG_wilcox_integrated.tsv"),sep="\t",row.name=F,quote = FALSE)

  DEG_wilcox_order %>% group_by(cluster) %>% top_n(n = ntop, wt = avg_log2FC) -> top_list
  pdf(paste0(Organ,"_DEG_wilcox_integrated.pdf"),length(table(Idents(scRNA.integrated.Find)))*0.5,length(table(Idents(scRNA.integrated.Find)))*0.2*ntop)
  p <- DotPlot(scRNA.integrated.Find, assay = 'RNA', features =  unique(top_list$gene)) + RotatedAxis() + coord_flip() + theme(axis.title.x=element_blank(),axis.title.y=element_blank())
  print(p)
  dev.off()

  }else{
    cat(Organ,"has",length(file_list),"(",file_list,")","qualified samples from",length(unique(scRNA@meta.data$Batch)),"(",unique(scRNA@meta.data$Batch),") Batch.\n")
  }
} else {
  cat("Please think if to integrate again!\n")
}
} else {
  cat(Organ,"has",length(file_list),"qualified samples:",file_list,".\n")
}
