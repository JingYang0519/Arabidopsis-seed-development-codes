library(Seurat)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyr)
library(dplyr)

rds <- readRDS("../pseudotime/0517.stage2-8.endo_rmCZE_pseudotime.rds")

# Normalize the data
rds <- NormalizeData(rds)

names(rds@meta.data)
unique(rds$endotype0514)
unique(rds$endoDetail0514)

data <- read.table("gene_info0619_for_bubblePlot.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

genes <- data$TG
gene_names <- data$Symbol

gene_id_name_map <- setNames(gene_names, genes)

exp <- as.data.frame(t(as.data.frame(rds@assays$RNA@data[genes,])))
meta <- rds@meta.data

merge.meta <- cbind(meta, exp)
colnames(merge.meta)[(ncol(meta) + 1):(ncol(meta) + length(genes))] <- gene_names

avg_exp <- merge.meta %>%
  group_by(endotype0514) %>%
  summarise(across(all_of(gene_names), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-endotype0514, names_to = "gene_name", values_to = "value")


anno_levels <- c("SE_CC", "SE_SynC", "SE", "MCE", "MCE-like", "PEN", 
                    "CZE","CZE_PEN",
                    "ESR", "ESR_PEN",
                    "ESR_MCE", "PMZ_type_I", "PMZ_type_II", "AL_type_I", "AL_type_II")

# Convert endotype0514 to a factor with the specified levels
merge.meta$endotype0514 <- factor(merge.meta$endotype0514, levels = anno_levels)

# Now proceed with the calculation of average expression and cell percentage as before
avg_exp <- merge.meta %>%
  group_by(endotype0514) %>%
  summarise(across(all_of(gene_names), ~mean(.x, na.rm = TRUE))) %>%
  pivot_longer(-endotype0514, names_to = "gene_name", values_to = "value")

cell_percentage <- merge.meta %>%
  group_by(endotype0514) %>%
  summarise(across(all_of(gene_names), ~mean(.x > 0, na.rm = TRUE))) %>%
  pivot_longer(-endotype0514, names_to = "gene_name", values_to = "percentage")

# Join the cell percentage to the average expression data
avg_exp <- left_join(avg_exp, cell_percentage, by = c("endotype0514", "gene_name"))

# Add gene ID and gene label columns
avg_exp <- avg_exp %>%
  mutate(gene_id = names(gene_id_name_map[match(gene_name, gene_names)]),
         gene_label = paste(gene_id, gene_name, sep = " | "))

# Convert gene_label column to a factor with the specified order
avg_exp$gene_label <- factor(avg_exp$gene_label, levels = paste(genes, gene_id_name_map[genes], sep = " | "))

# Define color gradient and create the bubble plot with smaller points
color_gradient <- scale_color_gradient(low = "#E3EEEF", high = "#C85D4D")
bubble_plot1 <- ggplot(avg_exp, aes(x = endotype0514, y = gene_label, size = percentage, color = value)) +
  geom_point() +
  color_gradient +
  scale_size_continuous(range = c(1, 6)) +  # Adjusted range for smaller points
  theme_minimal() +
  labs(x = "Endotype", y = "Gene ID | Gene Name", size = "Cell Proportion", color = "Expression Level") +
  #coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.major = element_blank(),  # Remove major grid lines
        panel.grid.minor = element_blank()   # Remove minor grid lines
  )


bubble_plot1

ggsave("bubble_plot.yGene_xCellType_0619.pdf", plot = bubble_plot1, width = 8, height = 9, units = "in")