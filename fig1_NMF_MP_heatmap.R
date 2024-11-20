### Analysis on cNMF results
### Raimi Chen 2023/06/08

### 0. Library -----------------------------------------------------------------
print("Begin settings...")
library(tidyverse)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(factoextra)
library(viridis)

top_gene.num <- 100 # Filter program top gene threshold
usage_filter <- -1 # -1 for no filer, >0 for filter
out <- "02.Metaprogram.rm24h.v3" # Output directory

if (!file.exists(out)) {
  dir.create(out)
}

### 1. QC ----------------------------------------------------------------------
dirs <- dir("01.Process/") # Input cNMF results directory
print("Begin QC...")
for (i in dirs) {
  # i <- dirs[1] # Test
  usage.df <- read.table(paste0("01.Process/",i,"/",list.files(paste0("01.Process/",i))[grepl("usages",list.files(paste0("01.Process/",i)))]),
                         header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  colnames(usage.df) <- paste0(i,"_P",1:dim(usage.df)[2])
  usage.df <- usage.df / rowSums(usage.df) # Normalize usage (sum to 1)
  write.table(usage.df, paste("02.Metaprogram.rm24h.v3","/",i,"_program.usage.norm.txt",sep = ""),
              quote = F, sep = "\t",row.names = T, col.names = T)
    
  #QC1
  tmpdf1 <- gather(usage.df,"program","ratio")
  median_df1 <- tmpdf1 %>%
    group_by(program) %>%
    summarize(median_value = median(ratio))
  p1 <- ggplot(tmpdf1, aes(x=program,y=ratio))+geom_boxplot(outlier.shape = NA) + # +geom_jitter(color="#4E5C68",alpha=0.2,width = 0.01)+
    labs(title = i) + theme_classic() +
    geom_text(data = median_df1, aes(x = program, y = median_value, label = round(median_value, 4)),
              vjust = -1) + geom_hline(yintercept = usage_filter, linetype = "dashed", color = "red") +
    geom_text(aes(x = 0.8, y = usage_filter, label = paste0("Cutoff: ",usage_filter)), color = "red", vjust = -0.5)
  
  pdf(paste(out,"/",i,"_program.usage.norm.QC.pdf",sep = ""),14,10)
  print(p1)
  dev.off()

  #score
  score.df <- read.table(paste0("01.Process/",i,"/",list.files(paste0("01.Process/",i))[grepl("gene_spectra_score",list.files(paste0("01.Process/",i)))]),
                         header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  score.df <- as.data.frame(t(score.df))
  colnames(score.df)=paste0(i,"_P",1:dim(score.df)[2])
  topn.df <- as.data.frame(matrix(nrow = top_gene.num, ncol = ncol(score.df)))
  colnames(topn.df)=colnames(score.df)
    
    for (k in colnames(score.df)) {
      tmpv=score.df[,k]
      names(tmpv)=rownames(score.df)
      topn.df[,k]=names(rev(tail(sort(tmpv),top_gene.num)))
    }
    
    #save
    write.table(topn.df, file = paste(out,"/",i,"_program.Zscore.top",top_gene.num,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
    score.df$gene=rownames(score.df)
    write.table(score.df,file = paste(out,"/",i,"_program.Zscore.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
}

print("Please do heatmap maunually!")

q()
### 3. Combining programs ------------------------------------------------------
 
 print("Begin settings...")
 library(tidyverse)
 library(ggplot2)
 library(ComplexHeatmap)
 library(circlize)
 library(factoextra)
 library(viridis)
 
 top_gene.num <- 100 # Filter program top gene threshold
 usage_filter <- -1 # -1 for no filer, >0 for filter
 out <- "02.Metaprogram.rm24h.v3" # Output directory
 
dirs <- dir("01.Process/") # Input cNMF results directory
check.usage=data.frame()
for (i in dirs) {
  # i <- dirs[1] # Test
  usage.file <- paste(out,"/",i,"_program.usage.norm.txt",sep = "")
  usage.df <- read.table(usage.file,header = T,row.names = 1,sep = "\t",stringsAsFactors = F)
  check.usage <- rbind(check.usage,as.data.frame(apply(usage.df, 2, median)))
}
  
colnames(check.usage)=c("median_ratio")
check.usage$sample_programs <- rownames(check.usage)
check.usage <- check.usage%>%arrange(median_ratio)
check.usage$sample_programs=factor(check.usage$sample_programs,levels = check.usage$sample_programs)
  
  linex=sum(check.usage$median_ratio < usage_filter)
  p2 <- ggplot(check.usage, aes(x=sample_programs, y=median_ratio))+geom_point()+
    geom_hline(yintercept = usage_filter,color="red")+
    geom_vline(xintercept = linex+0.5,color="red")+
    geom_text(aes(x = 3, y = usage_filter, label = paste0("Cutoff: ",usage_filter)), color = "red", vjust = -0.5) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  pdf("02.Metaprogram.rm24h.v3/00.check.usage.pdf",14,6)
  print(p2)
  dev.off()
  
  maybe.bg <- as.character(check.usage$sample_programs[check.usage$median_ratio < usage_filter])
  
print("Do heatmap...")
  
  all.score.df=data.frame()
  all.score.topn.df=data.frame()
  
  for (i in dirs) {
    score.file=paste(out,"/",i,"_program.Zscore.txt",sep = "")
    score.df=read.table(score.file,header = T,sep = "\t",stringsAsFactors = F)
    if (i==dirs[1]) {all.score.df=score.df}
    if (i!=dirs[1]) {
      all.score.df=all.score.df%>%inner_join(score.df,by="gene")
    }
    
    score.topn.file=paste(out,"/",i,"_program.Zscore.top",top_gene.num,"gene.txt",sep = "")
    score.topn.df=read.table(score.topn.file,header = T,sep = "\t",stringsAsFactors = F)
    if (i==dirs[1]) {all.score.topn.df=score.topn.df}
    if (i!=dirs[1]) {
      all.score.topn.df=cbind(all.score.topn.df,score.topn.df)
    }
  }
  
  rownames(all.score.df)=all.score.df$gene
  all.score.df$gene=NULL
  all.score.df=all.score.df[rowSums(is.na(all.score.df)) == 0,] 
  all.score.rm.df=all.score.df[,setdiff(colnames(all.score.df),maybe.bg)] 
  all.score.rm.df.cor=cor(all.score.rm.df,method = "spearman")
  
  # all.score.rm.df.cor[all.score.rm.df.cor < cor_min]=cor_min
  # all.score.rm.df.cor[all.score.rm.df.cor > cor_max]=cor_max
  
  # magma(100)
  # 
  
  ### Add 2 row annotation
  # colnames(all.score.rm.df.cor)
ha_column2 = HeatmapAnnotation(df = data.frame(TimePoints = c(rep("1.9-11h", 13),
                                                                rep("3.28h", 26),
                                                                rep("4.48h", 16),
                                                                rep("5.Globular", 11),
                                                                rep("6.Heart", 8),
                                                                rep("7.Torpedo", 10),
                                                                rep("8.Bent",11),
                                                                rep("9.Cotyledon",8),
                                                                rep("90.Late-cotyledon",11)),
                                                 CellType = c(rep("1.ReproductiveCell_1_EggCell",3),rep("2.ReproductiveCell_1_Sperm",3),rep("3.ReproductiveCell_2_CentralCell",2),rep("4.ReproductiveCell_2_SynergidCell",3),rep("9.SeedCoat",2),
                                                              rep("7.Embryo",4),rep("8.Endosperm",6),rep("2.ReproductiveCell_1_Sperm",7),rep("9.SeedCoat",3),rep("5.ZygoteB",6),
                                                              rep("7.Embryo",2),rep("8.Endosperm",4),rep("9.SeedCoat",7),rep("90.Suspensor",3),
                                                              rep("7.Embryo",4),rep("8.Endosperm",2),rep("9.SeedCoat",2),rep("90.Suspensor",3),
                                                              rep("7.Embryo",3),rep("8.Endosperm",2),rep("9.SeedCoat",3),
                                                              rep("7.Embryo",6),rep("8.Endosperm",2),rep("9.SeedCoat",2),
                                                              rep("7.Embryo",4),rep("8.Endosperm",3),rep("9.SeedCoat",4),
                                                              rep("7.Embryo",3),rep("8.Endosperm",3),rep("9.SeedCoat",2),
                                                              rep("7.Embryo",6),rep("8.Endosperm",3),rep("9.SeedCoat",2))),
                                col = list(TimePoints = c("1.9-11h" = "#FDE725FF",
                                                          "3.28h" = "#B4DE2CFF",
                                                          "4.48h" = "#6DCD59FF",
                                                          "5.Globular" = "#35B779FF",
                                                          "6.Heart" = "#1F9E89FF",
                                                          "7.Torpedo" = "#26828EFF",
                                                          "8.Bent" = "#31688EFF",
                                                          "9.Cotyledon" = "#3E4A89FF",
                                                          "90.Late-cotyledon" = "#482878FF"),
                                           CellType = c("1.ReproductiveCell_1_EggCell" = "#7CC08C",
                                                        "2.ReproductiveCell_1_Sperm" = "#bce07e",
                                                        "3.ReproductiveCell_2_CentralCell" = "#a3b9f0",
                                                        "4.ReproductiveCell_2_SynergidCell" = "#FDE377",
                                                        "9.SeedCoat" = "#F4A370",
                                                        "8.Endosperm" = "#63a4cf",
                                                        "5.ZygoteB" = "#3ca372",
                                                        "7.Embryo" = "#357B6E",
                                                        "90.Suspensor" = "#d962c2")))
  # Correlation heatmap
  set.seed(123)
  h2 <- Heatmap(as.matrix(all.score.rm.df.cor), 
                #col = colorRamp2(c(-1,0.2,0.4,0.6,0.7,1.0), c("#FFFFFFFF","#FD9A6AFF","#E85362FF","#D6456CFF","#C03A76FF","#160F3BFF")),
                # col=rev(magma(20)),
                col = colorRamp2(c(0,0.1,0.2,0.3,0.5,0.7,0.8,0.9,1),c("#FFFFFF",rev(mako(8, alpha = 1, begin = 0.3, end = 1)))),
                # col = colorRampPalette(colors = c("#46466A", "#5083A4", "#6ECACD","white","#e65a5d","red","black"))(100),
                clustering_method_rows = "complete",
                clustering_method_columns = "complete",
                top_annotation = ha_column2,
                cluster_rows = T, cluster_columns = T, name = "Spearman corr.",
                show_row_names = F, show_column_names = F)
  h2 <- draw(h2)
  col.list <- column_order(h2) 
  col.order <- colnames(as.matrix(all.score.rm.df.cor)[,col.list])
  write.table(col.order, "02.Metaprogram.rm24h.v3/01.program_spearman_cor.col.order.tsv", quote = F, sep = "\t", row.names = F, col.names = F)
  
  pdf("02.Metaprogram.rm24h.v3/01.program_spearman_cor.heatmap.pdf",14,10)
  print(h2)
  dev.off()

  write.table(all.score.rm.df.cor,file = paste("02.Metaprogram.rm24h.v3/01.cor_heatmap_data.txt",sep = ""),quote = F,sep = "\t",row.names = T,col.names = T)

  all.score.topn.rm.df=all.score.topn.df[,setdiff(colnames(all.score.topn.df),maybe.bg)] # Filter background noise
  write.table(all.score.topn.rm.df,file = paste("02.Metaprogram.rm24h.v3/01.program_top",top_gene.num,"gene.txt",sep = ""),quote = F,sep = "\t",row.names = F,col.names = T)
