library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyr)
library(dplyr)

# Load the RDS file that contains the Seurat object with Pseudotime information
rds <- readRDS("pseudotime/0517.stage2-8.endo_rmCZE_pseudotime.rds")

# Normalize the data
rds <- NormalizeData(rds)


# Specify the genes of interest
genes <- c("AT1G55600", "AT5G60440", "AT1G65330", "AT1G49770", "AT3G26744", "AT1G75080")

# Extract expression data for the genes of interest
exp <- as.data.frame(t(as.data.frame(rds@assays$RNA@data[genes,])))

# Extract metadata including Pseudotime
meta <- rds@meta.data

# Combine metadata with expression data
merge.meta <- cbind(meta, exp)

# Rename the columns with gene names
colnames(merge.meta)[(ncol(meta) + 1):(ncol(meta) + length(genes))] <- c("MINI3", "AGL62", "PHE1", "ZOU", "ICE1", "BZR1")

# Normalize expression data by scaling between 0 and 1
dat2 <- merge.meta
for (gene in c("MINI3", "AGL62", "PHE1", "ZOU", "ICE1", "BZR1")) {
  dat2[[gene]] <- (dat2[[gene]] - min(dat2[[gene]], na.rm = TRUE)) / (max(dat2[[gene]], na.rm = TRUE) - min(dat2[[gene]], na.rm = TRUE))
}

# Convert data from wide to long format
dat_long <- gather(dat2, key = "gene", value = "Expression", c("MINI3", "AGL62", "PHE1", "ZOU", "ICE1", "BZR1"))


# Plot the normalized expression data
p <- ggplot(dat_long, aes(x = Pseudotime, y = Expression, color = gene, fill = gene)) + 
  #geom_point(alpha = 0.5) +
  geom_smooth() +
  theme_minimal() +
  ggtitle("Normalized Expression Over Pseudotime")
p

# Save the plot as a PDF file
pdf("Normalized_Expression_Pseudotime.pdf", width = 15, height = 9)
print(p)
dev.off()