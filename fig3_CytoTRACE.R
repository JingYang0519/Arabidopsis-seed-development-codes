# Rscript 9-11h.28h.Syn.Cen.End.correct.rds 10 FALSE scRNA_plot/Pseudotime/2402 9-11h_28h



args <- commandArgs (T)
if (length(args) != 5) {
	stop ("Usage: Rscript CytoTRACE.R <rds> <ncpu> <enableFast: TRUE or FALSE> <outdir> <sample_name>\n")
}

library(CytoTRACE)
library(Seurat)
library(dplyr)
obj <- readRDS(args[1])
setwd(args[4])
#run CytoTRACE
expr_matrix <- as.matrix(obj@assays$RNA@data)
#result <- CytoTRACE(expr_matrix,enableFast = FALSE, ncores = 6) 
if (args[3] == "TRUE"){
	result <- CytoTRACE(expr_matrix,enableFast = TRUE, ncores = as.numeric(args[2])) 
}else{
	if (args[3] == "FALSE"){
		result <- CytoTRACE(expr_matrix,enableFast = FALSE, ncores = as.numeric(args[2]))
	}
}

save(result,file=paste0(args[4],"/",args[5],"_CytoTRACE_result.rda"))
#pdf (file=paste0(args[3],"/",args[4],"_CytoTRACE_15Genes.pdf"))
plotCytoGenes(result, numOfGenes = 15, outputDir = paste(args[4],"/",args[5],".",sep=""))
#dev.off()

# anno <- obj$subcelltype
anno <- obj$Cluster
anno <- as.character(anno)
names(anno) <- rownames(obj@meta.data)
plotCytoTRACE(result, emb = obj@reductions$umap@cell.embeddings,outputDir=paste(args[4],"/",args[5],".",sep=""),phenotype = anno)

	
#pdf (file=paste0(args[4],"/",args[5],"_CytoTRACE_15Genes.pdf"))
#p
#dev.off()
