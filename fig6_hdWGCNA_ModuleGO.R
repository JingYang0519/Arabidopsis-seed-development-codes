library(DOSE)
library(clusterProfiler)
library(pathview)
library(igraph)
library(org.At.tair.db)   
library(Seurat)
library(cowplot)
library(patchwork)
library(WGCNA)
library(hdWGCNA)
library(monocle3)
library(SeuratWrappers)
library(ggplot2)
library(viridis)
library(dplyr)
library(enrichplot)

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

# seurat_obj <- readRDS("f05.AnnoLabel.ModuleEigengenesConnectivity.rds")
# 
# MEs <- GetMEs(seurat_obj)
# modules <- GetModules(seurat_obj)
# mods <- levels(modules$module)
# mods <- mods[mods != "grey"]
# #mods <- mods[mods != "black"]
# #mods <- mods[mods != "magenta"]
# 
# list <- data.frame()
# for (cur_mod in mods) {
#   print(cur_mod)
#   cur <- subset(modules, module == cur_mod)
#   cur <- cur[, c("gene_name", paste0("kME_", cur_mod))]
#   colnames(cur)[2] <- "var"
#   list.tmp <- cur %>% arrange(desc(var)) %>% .$gene_name
#   list.tmp <- data.frame(Module = rep(cur_mod, length(list.tmp)), TopGenes = list.tmp, Rank = 1:length(list.tmp))
#   list <- rbind(list, list.tmp)
# }
# write.table(list, "Module.top.genelist.all.v2.tsv", sep = "\t", quote = F, row.names = F, col.names = T)


# Input data -------------------------------------------------------------------

ModuleGene <- read.delim("Module.top.genelist.all.v2.tsv", header = T)
type <- unique(ModuleGene$Module)
for (i in type) {
  gene <- unique(ModuleGene$TopGenes[which(ModuleGene$Module == i)])
  genelits = bitr(gene, 
                  fromType = "TAIR",
                  toType = "ENTREZID", 
                  OrgDb = "org.At.tair.db")
  id <- as.vector(genelits[,2])
  
  ### Perform GO
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.At.tair.db,
    ont           = "MF",    # options : "CC" "BP" "MF" "ALL"
    pAdjustMethod = "BH",
    readable      = FALSE,
    qvalueCutoff = 0.05, pvalueCutoff = 0.05) 
  
  write.table(ego.all, paste("Enrichment.MF.v2/GO.results.",i,".all.txt", sep = ""), 
              row.names = FALSE, sep = "\t")
  
  pl <- treeplot(pairwise_termsim(ego.all), hclust_method = "average", nCluster = 2)
  
  pdf("Test.GOSimilarity.pdf",20,5)
  print(pl)
  dev.off()
  
df1 <- as.data.frame(ego.all)
# df2 <- as.data.frame(ekk.all)
nn <- 15 # Set top nn

# GO top i
m <- dim(df1)[1]
if (m > nn) {
  top1 <- head(df1[,c(2,5)], nn)
  top1[,2] <- -log10(top1[,2])
  m <- nn
}else{
  top1 <- df1[,c(2,5)]
  top1[,2] <- -log10(top1[,2])
}

p1 <- ggplot(top1, aes(reorder(Description, pvalue), pvalue))+
  #geom_bar(stat = "identity", position = "dodge", fill = "#b50016", width = 0.5) + # Default bar color
  # geom_bar(stat = "identity", position = "dodge", fill = "#F16B00", width = 0.5) + # BP bar color
  geom_bar(stat = "identity", position = "dodge", fill = "#a912db", width = 0.5) + # BP bar color
  theme(axis.ticks.length = unit(0.2,'cm'),
        axis.text.x = element_text(size = 12, color = "black"),
        axis.text.y = element_text(size = 12, color = "black")) +
  guides(fill=guide_legend(title = NULL))+
  ggtitle(paste("Top ", m, " enriched GO terms", sep = "")) + coord_flip() +
  theme(axis.title = element_blank(), legend.position='none') +
  theme(panel.grid = element_blank()) +   
  theme(axis.line = element_line(size = 0.5, colour = "black")) +
  theme(strip.background = element_blank(), legend.position = "none", 
        panel.background = element_rect(fill = 'transparent')) + 
  theme(axis.text.x = element_text(size = 12, color = "black", face = "plain")) + 
  theme(axis.text.y = element_text(size = 12, color = "black", face = "plain"))

pdf(paste("Enrichment.MF.v2/", "GO_plot.",i,".pdf",sep = ""), width = 10, height = 5)
print(p1)
dev.off()

}

# Input data -------------------------------------------------------------------

ModuleGene <- read.delim("Module.top.genelist.all.v2.tsv", header = T)
type <- unique(ModuleGene$Module)
for (i in type) {
  gene <- unique(ModuleGene$TopGenes[which(ModuleGene$Module == i)])
  genelits = bitr(gene, 
                  fromType = "TAIR",
                  toType = "ENTREZID", 
                  OrgDb = "org.At.tair.db")
  write.table(genelits,paste("Enrichment.BP.v2/genelist.",i,".all.txt", sep = ""),
              sep = "\t",row.names = FALSE,quote = FALSE)
  id <- as.vector(genelits[,2])
  
  ### Perform GO
  ego.all <- enrichGO(
    gene          = id,
    keyType = "ENTREZID",
    OrgDb         = org.At.tair.db,
    ont           = "BP",    # options : "CC" "BP" "MF" "ALL"
    pAdjustMethod = "BH",
    readable      = FALSE,
    qvalueCutoff = 0.05, pvalueCutoff = 0.05) 
  
  write.table(ego.all, paste("Enrichment.BP.v2/GO.results.",i,".all.txt", sep = ""), 
              row.names = FALSE, sep = "\t")
  
  df1 <- as.data.frame(ego.all)
  # df2 <- as.data.frame(ekk.all)
  nn <- 15 # Set top nn
  
  # GO top i
  m <- dim(df1)[1]
  if (m > nn) {
    top1 <- head(df1[,c(2,5)], nn)
    top1[,2] <- -log10(top1[,2])
    m <- nn
  }else{
    top1 <- df1[,c(2,5)]
    top1[,2] <- -log10(top1[,2])
  }
  
  p1 <- ggplot(top1, aes(reorder(Description, pvalue), pvalue))+
    #geom_bar(stat = "identity", position = "dodge", fill = "#b50016", width = 0.5) + # Default bar color
    geom_bar(stat = "identity", position = "dodge", fill = "#F16B00", width = 0.5) + # BP bar color
    # geom_bar(stat = "identity", position = "dodge", fill = "#a912db", width = 0.5) + # BP bar color
    theme(axis.ticks.length = unit(0.2,'cm'),
          axis.text.x = element_text(size = 12, color = "black"),
          axis.text.y = element_text(size = 12, color = "black")) +
    guides(fill=guide_legend(title = NULL))+
    ggtitle(paste("Top ", m, " enriched GO terms", sep = "")) + coord_flip() +
    theme(axis.title = element_blank(), legend.position='none') +
    theme(panel.grid = element_blank()) +   
    theme(axis.line = element_line(size = 0.5, colour = "black")) +
    theme(strip.background = element_blank(), legend.position = "none", 
          panel.background = element_rect(fill = 'transparent')) + 
    theme(axis.text.x = element_text(size = 12, color = "black", face = "plain")) + 
    theme(axis.text.y = element_text(size = 12, color = "black", face = "plain"))
  
  pdf(paste("Enrichment.BP.v2/", "GO_plot.",i,".pdf",sep = ""), width = 10, height = 5)
  print(p1)
  dev.off()
  
}
