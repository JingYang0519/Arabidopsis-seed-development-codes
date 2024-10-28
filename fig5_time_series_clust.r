require(Seurat);  require(tidyr); require(dplyr); require(ggplot2)
require(stringr); require(data.table); require(RColorBrewer);
require(patchwork); require(ComplexHeatmap); set.seed(123)

source('plotting_theme.R')

recluster <- readRDS("pseudotime/0517.stage2-8.endo_rmCZE_pseudotime.rds")

scale_this <- function(x) { (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE) }

# Calculate HVGs
recluster <- ScaleData(FindVariableFeatures(NormalizeData(recluster),
    selection.method = 'vst', nfeatures=3000))
high_var_genes <- VariableFeatures(recluster)
length(high_var_genes)

# Write HVGs to a text file
write.table(high_var_genes, file = "0517.HVG3000_info.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Extract Pseudotime information from recluster
recluster$cellID <- colnames(recluster)
meta_data <- data.frame(cell = recluster$cellID, Pseudotime = recluster$Pseudotime)

exp_data <- GetAssayData(recluster, assay = "RNA", slot = "data")[high_var_genes, ]

exp_data <- as.matrix(exp_data)

exp_data <- t(exp_data)

# Merge meta and expression data
merge.meta <- cbind(meta_data, exp_data) %>%
  pivot_longer(cols = all_of(high_var_genes)) %>%
  mutate(new_bin = cut_number(Pseudotime, 100)) %>%
  data.frame

# Save RDS
#saveRDS(merge.meta, file = "0517.Endosperm_rmCZE_cluster_snRNA.pseudotime_100bins.rds")

# Group data by gene and pseudotime bin
merge.meta2 <- merge.meta %>%
  group_by(name, new_bin) %>%
  summarize(mean_value = mean(value, na.rm = TRUE),
            mean_ptime = mean(Pseudotime, na.rm = TRUE), .groups = 'drop') %>%
  data.frame
  
rna <- merge.meta2
# Perform smoothing on already averaged data
n_adjacent_windows <- 15 # 15 looks good
total_windows <- sort(unique(rna$new_bin))
n_windows <- length(total_windows)
start_edge_idx <- 1:n_adjacent_windows
end_edge_idx <- (length(total_windows)-(n_adjacent_windows-1)):length(total_windows)
start_edge_bins <- total_windows[start_edge_idx]
end_edge_bins <- total_windows[end_edge_idx]
used_srnas <- rbindlist(lapply(1:n_windows, function(i) {
    print(paste0('working on bin: ', total_windows[i]))

    # decide which windows to average over
    if (i %in% start_edge_idx) {# beginning edge
        place <- 'first'
        start_idx <- 1; end_idx <- i+n_adjacent_windows
    }

    if (i %in% end_edge_idx) { # end edge
        place <- 'end'
        start_idx <- i-n_adjacent_windows; end_idx <- n_windows
    }
    if (!(i %in% start_edge_idx) & !(i %in% end_edge_idx)) { # in between

        place <- 'middle'
        start_idx <- i-n_adjacent_windows; end_idx <- i+n_adjacent_windows
    }
    print(paste0('averaging over: ', start_idx, ' to ', end_idx))
    central_bin_ptime <- mean((rna %>% filter(new_bin == total_windows[i]))$mean_ptime)
    return(rna %>% filter(new_bin %in% total_windows[start_idx:end_idx]) %>% 
        group_by(name) %>%
		summarize(mean_value=mean(mean_value, trim=0.1),
            mean_ptime=central_bin_ptime, # sets ptime to central bin
            .groups='drop') %>% 
        mutate(new_bin=total_windows[i]) %>%
        select(name, new_bin, mean_value, mean_ptime))
})) 

rna_wide <- used_srnas %>%
    pivot_wider(id_cols=c('new_bin', 'mean_ptime'),
		names_from=name, values_from=mean_value)
df_mat <- rna_wide %>% select(-new_bin, -mean_ptime) %>% as.matrix
rownames(df_mat) <- rna_wide$new_bin

# maybe only keep highly variable genes
new_vars <- apply(df_mat, 2, sd)
df_mat <- scale(df_mat[,rev(order(new_vars))[1:2000]])

top2000_var_genes <- names(rev(sort(new_vars)))[1:2000]

write.table(top2000_var_genes, file = "0517.HVG1500_2nd.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)


# redo above but with var genes
rna_wide <- used_srnas %>%
    filter(name %in% names(rev(sort(new_vars)))[1:2e3]) %>%
    group_by(name) %>% mutate(mean_value=scale_this(mean_value)) %>% ungroup %>%
    pivot_wider(id_cols=c('new_bin', 'mean_ptime'),
		names_from=name, values_from=mean_value)
df_mat <- rna_wide %>% select(-new_bin, -mean_ptime) %>% as.matrix
rownames(df_mat) <- rna_wide$new_bin

# # # tyring slightly fancier clustering strategy
require(dtw); require(dtwclust)
clusts0 <- tsclust(t(df_mat), type='partitional',
	k=2:35, # to check different values
	distance='dtw_basic', centroid='pam', trace=T)	
saveRDS(clusts0, '0517.HVG1500.Endosperm_rm24h_rmCZE_cluster_scRNA.pseudotime_dtw_clust.rds')

k_clusters <- 5 #from the previous calculation

clusts <- data.frame(name=colnames(df_mat), gene_cluster=clusts0[[4]]@cluster) %>%
	# for now using same previous order, but might need to reorder
	merge(data.frame(gene_cluster=1:k_clusters, new_clusters=c(2,1,5,3,4))) %>%
	identity()

rnasc <- used_srnas %>% # sliding window values
    filter(name %in% colnames(df_mat)) %>% group_by(name) %>%
	mutate(mean_value=scale_this(mean_value)) %>% ungroup %>%
    merge(clusts) %>% identity()

# to define the new order of clusters
rnasc %>% group_by(gene_cluster, new_bin) %>%
    summarize(mean_value=mean(mean_value), .groups='drop') %>%
		group_by(gene_cluster) %>%
    # ok ordering clusters by the first time it reaches 80% of max value
    summarize(max_value=max(mean_value),
        ugh=new_bin[which(mean_value > (.8*max_value))[1]],
        .groups='drop') %>%
    arrange(ugh, -max_value)

# # Let's make a version with an averged pattern
some_colors <- RColorBrewer::brewer.pal(n=k_clusters, name='Dark2')
names(some_colors) <- as.character(1:k_clusters)

rnasc <- used_srnas %>%

    filter(name %in% colnames(df_mat)) %>%

    group_by(name) %>%
	
	mutate(mean_value=scale_this(mean_value)) %>%
	
    ungroup %>%

    merge(clusts) %>%

    identity()
	

rnasc %>%

    group_by(gene_cluster, new_bin) %>%

    summarize(mean_value=mean(mean_value), .groups='drop') %>%

    group_by(gene_cluster) %>%

    summarize(max_value=max(mean_value),
        ugh=new_bin[which(mean_value > (.8*max_value))[1]],
        .groups='drop') %>%

    arrange(ugh, -max_value)
	

new_cluster_order <- rnasc %>%
    group_by(gene_cluster, new_bin) %>%
    summarize(mean_value = mean(mean_value), .groups = 'drop') %>%
    group_by(gene_cluster) %>%
    summarize(max_value = max(mean_value),
              ugh = new_bin[which(mean_value >= (0.8 * max_value))[1]],
              .groups = 'drop') %>%
    arrange(ugh, -max_value) %>%
    mutate(new_order = 1:n())


clusts <- clusts %>%
  left_join(new_cluster_order %>% select(gene_cluster, new_order), by = "gene_cluster") %>%
  arrange(new_order) %>%
  mutate(new_clusters = new_order)


write.table(clusts, '0517.HVG1500_C5.pseudotime_clusters_ordered.txt', quote = F, row.names = F, sep = '\t')


some_colors <- RColorBrewer::brewer.pal(n=k_clusters, name='Dark2')

names(some_colors) <- as.character(1:k_clusters)

head(rna_wide)
rnasc %>% write.table('0517.HVG1500_C5.pseudotime_clusters.txt',
    quote=F, row.names=F, sep='\t')
	

p1 <- rnasc %>% 
    ggplot(aes(
        x=mean_ptime,  
        y=mean_value,  
        color=factor(new_clusters), 
        group=factor(new_clusters)  
    )) +

    geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"),
        size=0.5, se=F) + 

    scale_color_manual("Cluster", values=some_colors) +

    guides(color=guide_legend(override.aes=list(size=1.5), keywidth=unit(5, 'mm'))) +

    theme(legend.position='top', panel.grid.major.x=element_line(size=0.1,
            color=scales::alpha('black', 0.2)))  +

    xlab('Pseudotime') + ylab("Scaled expression") + 
    theme(axis.text.x = element_blank(),  
            axis.ticks.x = element_blank())  

ggsave('0517.HVG1500_C5.pseudotime_cluster_fit.pdf', p1, height = 37,
    width = 83, units='mm', useDingbats=FALSE)
	
# make a heatmap
some_colors2 <- some_colors
names(some_colors2) <- as.character(1:k_clusters)
ha <- rowAnnotation(Clusters=as.character(clusts$new_clusters),
    col = list(Clusters=some_colors2), show_annotation_name = F, show_legend=F)


data <- read.table("106Endosperm_markers.xls", header = TRUE, sep = "\t")


InterestGene <- data$Gene
GeneLabel <- paste(data$Symbol,data$Notes, sep = "|")


some_colors2 <- some_colors
names(some_colors2) <- as.character(1:k_clusters)

ha <- rowAnnotation(Clusters=as.character(clusts$new_clusters),
    col = list(Clusters=some_colors2), show_annotation_name = F, show_legend=F)
	

retainedInterestGene <- InterestGene[InterestGene %in% clusts$name]
retainedGeneLabel <- GeneLabel[InterestGene %in% clusts$name]
length(retainedInterestGene)


print(retainedInterestGene)
print(retainedGeneLabel)


retainedInterestGene_idx <- match(retainedInterestGene, clusts$name)


ha2 <- rowAnnotation(
    foo = anno_mark(at = retainedInterestGene_idx, labels = retainedGeneLabel,
        link_width = unit(4, "mm"),
        labels_gp = grid::gpar(fontfamily='Helvetica', fontsize=6))
)


p <- Heatmap(t(df_mat)[clusts$name,], cluster_columns=F,
    show_column_names=F,
    right_annotation = ha2,
    column_names_gp = grid::gpar(fontsize = 1, fontfamily="Helvetica"),
    name='Expression',
    heatmap_legend_param = list(legend_gp = gpar(fontsize = 6),
        legend_title_gp = gpar(fontsize = 7, fontface = "bold")),
    row_split=as.character(clusts$new_clusters), 
    row_gap = unit(1, "mm"), left_annotation=ha,
    use_raster = F, 
    show_row_names=F, cluster_rows=F
)
p 


pdf("0530markers_refine.HVG1500_C5.pseudotime_heatmap.pdf")
draw(p)
dev.off()
