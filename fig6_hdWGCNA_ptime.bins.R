library(DOSE)
library(clusterProfiler)
library(pathview)
library(ggplot2)
library(igraph)
library(org.At.tair.db)   
keytypes(org.At.tair.db)
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


# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("f05.AnnoLabel.ModuleEigengenesConnectivity.rds")

### 1. Load data ---------------------------------------------------------------
mydata <- data.frame(Cell = row.names(seurat_obj@meta.data), Stage = seurat_obj$Organ, Pseudotime = seurat_obj$pseudotime)
# mydata <- read.delim("ptime.dt.cell_type.tsv", header = T)
newdata <- data.frame(Pseudotime = seq(min(mydata$Pseudotime),max(mydata$Pseudotime),length.out = 21))
mydata$Column <- cut(mydata$Pseudotime, breaks = newdata$Pseudotime, include.lowest = TRUE, right = T)

cell_counts <- mydata %>%
  group_by(Stage,Column) %>%
  summarise(cell_count = n())

### 2. Stacked histogram -------------------------------------------------------
## Calculate percentage
data <- cell_counts
data.sum <- aggregate(data$cell_count, by = list(data$Column), sum) # View(data.sum)
data.merge <- merge(data, data.sum, all.x = TRUE, all.y = TRUE, by.x = "Column", by.y = "Group.1") # View(data.merge)
data.merge$percentage <- 100*(data.merge$cell_count/data.merge$x) # View(data.merge)

## Plot
pal <- c("#35B779FF","#1F9E89FF","#26828EFF","#31688EFF","#3E4A89FF","#482878FF")
p <- ggplot(data.merge,aes(x = Column, 
                           y = percentage, 
                           fill = Stage)) + 
  scale_y_continuous(expand = c(0, 0)) + # plot start from the origin
  geom_bar(stat='identity', width = 1) + # No space between bars
  scale_fill_manual(values = pal) +
  #theme_minimal()
  #theme_classic() #+ RotatedAxis()
  theme_gray()

pdf("p09.Stacked.ptime.pdf",9,5)
print(p)
dev.off()