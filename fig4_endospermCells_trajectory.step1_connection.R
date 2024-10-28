library(Seurat)
library(dplyr)
library(ggplot2)
library(Matrix)
library(htmlwidgets)
library(plotly)
library(monocle3)
library(FNN)
library(future)
library(future.apply)

source("Fig5_endospermCells_trajectory.help_code.R")

args = commandArgs(trailingOnly=TRUE)
kk = as.numeric(args[1])

work_path = "/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202404.Cell_Trajectory/endo0517"
rds <-readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202404.Cell_Trajectory/endo0514/refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")


time_point <- unique(rds$Organ)
cell_types <- unique(rds$celltype0515) 
time_point

time_1 = time_point[kk]
time_2 = time_point[kk+1]
print(time_1)
print(time_2)


subset_data = subset(rds, Organ == time_1 | Organ == time_2)

obj <- subset_data
obj.list <- SplitObject(obj, split.by = "Organ")
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = obj.list)
obj.list <- future_lapply(X = obj.list, FUN = function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
})

anchors <- FindIntegrationAnchors(object.list = obj.list, reduction = "rpca", dims = 1:50)
obj.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
obj.integrated <- ScaleData(obj.integrated, verbose = FALSE)
obj.integrated <- RunPCA(obj.integrated, npcs = 30, verbose = FALSE)

library(umap)

pca_results <- Embeddings(obj.integrated, "pca")

umap_results <- umap::umap(pca_results[, 1:30], n_components = 3, min_dist = 0.75)

if (is.list(umap_results) && "layout" %in% names(umap_results)) {

  colnames(umap_results$layout) <- paste0("UMAP_", seq_len(ncol(umap_results$layout)))
  
  umap_embedding <- CreateDimReducObject(embeddings = umap_results$layout, key = "UMAP_", assay = DefaultAssay(obj.integrated))
  
  obj.integrated[["umap"]] <- umap_embedding
} else {
  stop("The UMAP results do not have the expected format.")
}

obj.integrated[["umap"]] <- umap_embedding

saveRDS(obj.integrated, paste0(work_path, "/", time_1, "_", time_2, ".rds"))

emb = data.frame(Embeddings(object = obj.integrated, reduction = "umap"))
saveRDS(emb, file=paste0(work_path, "/", time_1, "_", time_2, "_umap3.rds"))


p <- DimPlot(obj.integrated, reduction = "umap", label = TRUE, split.by = "celltype0515",group.by="Organ") 
ggsave(paste0(work_path, "/", time_1, "_", time_2, "_UMAP.pdf"), plot = p, width = 40, height = 6)

kk = as.numeric(args[1])
time_i = time_point[kk]
time_j = time_point[kk+1]
print(time_i)
print(time_j)

subset_data = subset(rds, Organ == time_i | Organ == time_j)

if ("umap" %in% names(subset_data@reductions)) {
  umap_results <- subset_data@reductions$umap@cell.embeddings
  # Convert the UMAP results to a data.frame
  umap_df <- as.data.frame(umap_results)
} else {
  stop("UMAP results do not exist in the Seurat object.")
}


anno1 = subset(rds, Organ == time_i)
anno1$Anno = as.vector(anno1$celltype0515)
anno1@meta.data$day <- "pre"
anno1 = anno1[[]][,c("day", "Anno")]
anno1$stage = time_i

anno2 = subset(rds, Organ == time_j)
anno2$Anno = as.vector(anno2$celltype0515)
anno2@meta.data$day <- "nex"
anno2 = anno2[[]][,c("day", "Anno")]
anno2$stage = time_j

anno = rbind(anno1, anno2)

emb <- umap_df

if(nrow(emb) != nrow(anno)){
    print("Error!")
    print(xxx) 
}

anno = anno[rownames(emb),]

res = createLineage_Knn(emb, anno,  k_neigh = 5)

saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap.rds"))

for(k_i in c(8, 10, 15, 20)){
    res = createLineage_Knn(emb, anno,  k_neigh = k_i)
    saveRDS(res, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap_k_", k_i, ".rds"))
}

replication_times=500
dat = res
state_1 = row.names(dat[[1]])
state_2 = names(dat[[1]])
tmp_1 = matrix(NA, nrow(dat[[1]]), ncol(dat[[1]]))
for(i in 1:nrow(dat[[1]])){
    for(j in 1:ncol(dat[[1]])){
        xx = NULL
        for(k in 1:replication_times){
            xx = c(xx, dat[[k]][i,j])
        }
        tmp_1[i,j] = median(xx[!is.na(xx)])
    }
}
tmp_1 = data.frame(tmp_1)
row.names(tmp_1) = state_1
names(tmp_1) = state_2

write.csv(tmp_1, paste0(work_path, "/", time_i, "_", time_j, "_Knn_umap.csv"))


