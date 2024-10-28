library(Seurat)
library(ComplexHeatmap)
library(circlize)
library(dplyr)

key_genes <- read.table("cellTrajectory.gene_expPlot.txt", header = TRUE, stringsAsFactors = FALSE, sep = "\t")

celltrajectory <- readRDS("/jdfsbjcas1/ST_BJ/P21Z28400N0234/yangjing7/01.Proj/202404.Cell_Trajectory/endo0514/refined_anno0514.endosperm.stage1-9.for_celltrajectory.rds")
celltrajectory@meta.data$celltype0515<-gsub("Endo_CZE","CZE",celltrajectory@meta.data$celltype0515)

rds <- subset(celltrajectory, celltype0515 != "Sp")

gene_to_symbol <- setNames(key_genes$Symbol, key_genes$gene)

key_genes_unique <- key_genes %>% distinct(gene, .keep_all = TRUE)
genes <- key_genes_unique$gene
exp_data <- GetAssayData(rds, assay = "RNA", slot = "data")[genes, ]

exp_data_matrix <- as.matrix(exp_data)


rownames(exp_data_matrix) <- genes

meta_data <- rds@meta.data

exp_data_df <- data.frame(t(exp_data_matrix), check.names = FALSE)

exp_data_df$cell <- rownames(exp_data_df)
meta_data$cell <- rownames(meta_data)

exp_data_df <- merge(exp_data_df, meta_data, by = "cell")

avg_expression <- exp_data_df %>%
  group_by(Organ, celltype0515) %>%
  summarise(across(all_of(genes), mean, na.rm = TRUE), .groups = 'drop')

numeric_columns <- sapply(avg_expression, is.numeric)

avg_expression_z <- as.data.frame(scale(avg_expression[, numeric_columns]))
rownames(avg_expression_z) <- paste(avg_expression$Organ, avg_expression$celltype0515, sep = ":")

avg_expression_z_matrix =avg_expression_z
# Transpose avg_expression_z_matrix if genes are currently columns
# Ensure that the row names of avg_expression_z_matrix are gene IDs
avg_expression_z_matrix <- t(avg_expression_z_matrix)

# Map the gene IDs to Symbols
gene_symbols <- gene_to_symbol[rownames(avg_expression_z_matrix)]

# Handle missing gene IDs by setting their symbols to the gene ID itself
gene_symbols[is.na(gene_symbols)] <- rownames(avg_expression_z_matrix)[is.na(gene_symbols)]

# Ensure that the length of gene_symbols matches the number of rows in avg_expression_z_matrix
if (length(gene_symbols) != nrow(avg_expression_z_matrix)) {
  stop("The number of gene symbols does not match the number of rows in the expression matrix.")
}

# Create the row annotation
row_anno <- rowAnnotation(Symbols = anno_text(gene_symbols))

trajectory_colors <- c(
  "#FBE5B8",   #CC
  "#c2dcbf",   #Syn
  "#ef6547",   #SE_CC
  "#E11A1D",   #SE_SynC
  "#FF8C00",   #SE
  "#C2A289",   #MCE
  "#A46D2C",   #MCE-like
  "#9A9AF8",   #PEN 
  "#dbdcde",   #CZE
  "#A4A5A7",   #CZE_PEN
  "#FAD1E0",   #ESR
  "#F683D8",   #ESR_PEN
  "#EF5D8C",   #ESR_MCE
  "#65C2A4",   #PMZ_type_I
  "#1E803D",   #PMZ_type_II
  "#9FB2DA",   #AL_type_I
  "#1B79AF"    #AL_type_II
)


anno_levels1 <- c("CC","Syn","SE_CC", "SE_SynC", "SE", "MCE", "MCE-like", "PEN", 
                  "CZE","CZE_PEN","ESR", "ESR_PEN",
                  "ESR_MCE", "PMZ_type_I", "PMZ_type_II", "AL_type_I", "AL_type_II")

rds@meta.data$celltype0515 <- factor(rds@meta.data$celltype0515, levels = anno_levels1)

color_map1 <- setNames(trajectory_colors, anno_levels1)

devStage <- c("#caeac2","#ccd6bc", "#A7A9BC","#ECB884","#E4E45F","#4758A2","#E08D8B","#AF8CBB","#7A5A86")

anno_levels2 <- c("1_9_11h","2_28h", "3_48h", "4_globular",
             "5_heart", "6_torpedo", "7_bent", "8_cotyledon", "9_late_cotyledon")

rds@meta.data$Organ <- factor(rds@meta.data$Organ, levels = anno_levels2)

color_map2 <- setNames(devStage, anno_levels2)

anno_info <- strsplit(rownames(avg_expression_z), ":", fixed = TRUE)
anno_df <- data.frame(
  Organ = sapply(anno_info, `[`, 1),
  celltype0515 = sapply(anno_info, `[`, 2),
  stringsAsFactors = FALSE
)


ha <- HeatmapAnnotation(df = anno_df,
                        col = list(Organ = color_map2, celltype0515 = color_map1))

# Create a row annotation for gene symbols
row_anno <- rowAnnotation(Symbols = anno_text(gene_symbols))

# This assumes that 'anno_df' has the same order as the columns in 'avg_expression_z_matrix'
ordered_avg_expression_z_matrix <- avg_expression_z_matrix[, order(anno_df$Organ)]

# Draw the heatmap without clustering and with gene symbols as row annotation
pdf("Gene_Expression_Heatmap_organOrder1.pdf", width = 20, height = 15)
p <- Heatmap(ordered_avg_expression_z_matrix, name = "Z-score", top_annotation = ha,
        show_row_names = FALSE, show_column_names = FALSE,
        cluster_rows = FALSE, cluster_columns = FALSE) + row_anno
print(p)
dev.off()
p