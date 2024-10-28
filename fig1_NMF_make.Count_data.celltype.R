### Preparation for cNMF input raw count matrix
### Raimi Chen 2023/06/08

library(Seurat)

print("Loading RDS...")
out <- "/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/02.NMF/00.Count_data/"
rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/At/STO_embryo_scRNA/all_time_merge/allChr_v1/rmR1802/all_sample_subCluster_CellType/At_merge.seurat_clusters_renamed.rds")

if (!file.exists(out)) {
  dir.create(out)
}

print(unique(rds@meta.data$seurat_clusters_renamed))
print("Filtering cell types...")
rds@meta.data$seurat_clusters_renamed <- gsub("_c","C",rds@meta.data$seurat_clusters_renamed)
# rds@meta.data$seurat_clusters_renamed <- gsub("\\.","_",rds@meta.data$seurat_clusters_renamed)
rds@meta.data$seurat_clusters_renamed <- gsub("late_cotyledon","late-cotyledon",rds@meta.data$seurat_clusters_renamed)
ct.list <- unique(rds@meta.data$seurat_clusters_renamed)[!grepl(paste(c("Funiculus", "Unknown", "Zygote,Basal"), collapse = "|"), unique(rds@meta.data$seurat_clusters_renamed))]

print(paste0("Cell type number ",length(ct.list)))
print("Output raw count matrix...")
for (i in ct.list) {
#for (i in c("07.torpedo.Seed coat")) {
  print(i)
  sub.meta <- rds@meta.data[which(rds@meta.data$seurat_clusters_renamed == i),]
  sub.count <- as.data.frame(rds[["RNA"]]@counts[,row.names(sub.meta)])
  sub.count <- sub.count[rowSums(sub.count)>0,]
  sub.count <- data.frame(t(sub.count))
  write.table(sub.count, paste0(out,"/",i,".count.txt"),
              quote = F, sep = "\t", row.names = T, col.names = T)
}

print("Done!")
warnings()
