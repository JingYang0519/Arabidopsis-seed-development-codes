## load data -------------------------------------------------------------------
cyto <- read.delim("9-11h_28h.CytoTRACE_plot_table.txt", header = T)

cells1 <- read.delim("1_9_11h_Synergid.cell.Common.genelist.tsv", header = F)
cells2 <- read.delim("2-28h_Endosperm.cell.Common.genelist.tsv", header = F)

cells2$V1 <- gsub("03_embryo_28h", "02_embryo_28h", cells2$V1)

## Add typeII label -------------------------------------------------------
for (i in 1:nrow(cyto)) {
  if (row.names(cyto)[i] %in% c(cells1$V1,cells2$V1)) {
    cyto$Phenotype[i] <- paste0(cyto$Phenotype[i], "_typeII")
  }
}

cyto <- cyto[which(cyto$Phenotype != "04.28h_C11.1_Endosperm_typeII"),]

## Plot data -------------------------------------------------------------------
library(dplyr)
cyto.median <- cyto %>%
  group_by(Phenotype) %>%
  summarize(CytoTRACE_median = median(CytoTRACE))

arrange(cyto.median, desc(CytoTRACE_median))


cyto$Phenotype <- factor(cyto$Phenotype, levels = c("04.28h_C11.1_Endosperm","03.9-11h_C22_Central_cell","01.9-11h_C10.1_Central_cell",
                                                    "02.9-11h_C21_Syndergid_cell","05.28h_C13_Endosperm",
                                                    "02.9-11h_C21_Syndergid_cell_typeII","05.28h_C13_Endosperm_typeII"))
library(ggplot2)
p <- ggplot(cyto, aes(x=Phenotype, y=CytoTRACE)) + 
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf("CytoTRACE.typeII.pdf",10,8)
p
dev.off()

## t.test significance ---------------------------------------------------------
celltype1 <- "02.9-11h_C21_Syndergid_cell"
celltype2 <- "02.9-11h_C21_Syndergid_cell_typeII"
paste0(celltype1," (n = ",length(which(cyto$Phenotype == celltype1)),")")
paste0(celltype2," (n = ",length(which(cyto$Phenotype == celltype2)),")")
wilcox.test(cyto$CytoTRACE[which(cyto$Phenotype == celltype1)],cyto$CytoTRACE[which(cyto$Phenotype == celltype2)])

celltype1 <- "05.28h_C13_Endosperm"
celltype2 <- "05.28h_C13_Endosperm_typeII"
paste0(celltype1," (n = ",length(which(cyto$Phenotype == celltype1)),")")
paste0(celltype2," (n = ",length(which(cyto$Phenotype == celltype2)),")")
wilcox.test(cyto$CytoTRACE[which(cyto$Phenotype == celltype1)],cyto$CytoTRACE[which(cyto$Phenotype == celltype2)])

celltype1 <- "03.9-11h_C22_Central_cell"
celltype2 <- "01.9-11h_C10.1_Central_cell"
paste0(celltype1," (n = ",length(which(cyto$Phenotype == celltype1)),")")
paste0(celltype2," (n = ",length(which(cyto$Phenotype == celltype2)),")")
wilcox.test(cyto$CytoTRACE[which(cyto$Phenotype == celltype1)],cyto$CytoTRACE[which(cyto$Phenotype == celltype2)])
