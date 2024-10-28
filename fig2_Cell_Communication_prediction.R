#------------------------------------------------#
#          Cell communication Analysis           #
#     Editor: Ruiying Chen   TEL:18995624500     #
#              Edited on: 2023.06.11             #
#------------------------------------------------#

# install.packages("/path/to/package.tar.gz", repos = NULL, type = "source")
args=commandArgs(T)

### Load Arabidopsis LR pair data
print("Preparing data...")
#load("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/03.CellComm/00.Data/LR_pair_ath.RDa") # LR pair
LR_pair <- read.delim("PlantPhoneDB.AT.revised.QBLRY.tsv", header = T)
LR_pair <- data.frame(Ligands = LR_pair$Ligands, Receptors = LR_pair$Receptors)

### Library
library(PlantPhoneDB)
library(Seurat)

### Perform cell communication
print("Loading Seurat object...")
rds <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/hanrui/At/STO_embryo_scRNA/all_time_merge/allChr_v1/rmR1802/all_sample_subCluster_CellType/At_merge.seurat_clusters_renamed.rds")

print("Subsetting Seurat object...")
sub_rds <- subset(rds, Organ == args[1])
sub_rds$seurat_clusters_renamed <- gsub("Reproductive_cell_1_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("Reproductive_cell_2_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("01_9-11h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("02_24h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("03_28h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("04_48h_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("05_globular_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("06_heart_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("07_torpedo_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("08_bent_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("09_cotyledon_","",sub_rds$seurat_clusters_renamed)
sub_rds$seurat_clusters_renamed <- gsub("10_late_cotyledon_","",sub_rds$seurat_clusters_renamed)

sub_rds <- subset(sub_rds, seurat_clusters_renamed %in% c("Seed_coat","Central_cell","Sperm","Synergid_cell","Egg_cell","Apical_cell","Zygote","Endosperm","Embryo","Suspensor"))

print("Performing SCTransform...")
sub_rds2 <- SCTransform(sub_rds, verbose = FALSE)

print("Calculating Ligand-receptor pairs...")
LRp <- LRscore(sub_rds2@assays$SCT@data, 
               LRdb=LR_pair, 
               cluster = sub_rds2$seurat_clusters_renamed, 
               #min.pct = 0.15,
               min.pct = 0.1,
               iterations = 100, 
               method='Average')

print("Saving Ligand-receptor pair info RDS...")
saveRDS(LRp, paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/03.CellComm/02.Results/LRp_CellType.revised/",args[1],".LRp_CellType.0.1.rds"))

LRp <- LRp[which(LRp$Pvalue<0.05),]
write.table(LRp, paste0("/jdfsbjcas1/ST_BJ/P21Z28400N0234/chenruiying/01.Programs/2023_Embryo/01.Scripts/03.CellComm/02.Results/LRp_CellType.revised/",args[1],".LRp_CellType.test.0.1.tsv"), quote = F, sep = "\t", row.names = F, col.names = T)
print("Please perform visualization in your disk!")
